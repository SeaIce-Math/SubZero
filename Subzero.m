close all; clear all;

paths
%% Set Flags

RIDGING=false; 

FRACTURES=false;

PERIODIC=false;

PACKING = false;

WELDING = false;

CORNERS = true;

COLLISION = true;

AVERAGE = false; %average certain values in time

RAFTING = false;

KEEP_MIN = false; %Retain floes in simulation that get smaller set threshold

ifPlot = true; %Plot floe figures or not?

ifPlotStress = false; %Plot floe figures with floes shaded by stress values

ifPlotStressStrain = false; %Create plot with homogenized stress and strain values

justfrac = false;

%% Initialize model vars

dt=10; %Time step in sec
height.mean = 0.25; %mean value of thickness for initial floes
height.delta = 0; %maximum deviation about the mean thickness if a distribution is desired


nDTpack = 500; %Time streps between creation of new floes
rho_ice=920; % kg/m3
[ocean, HFo, h0]=initialize_ocean(dt,nDTpack); %Define ocean currents

%Define Atmospheric winds
winds.x = ocean.Xocn; winds.y = ocean.Yocn;
U0 = 0; V0 = 0;
winds.u=U0*ones(size(ocean.Xocn));
winds.v=V0*ones(size(ocean.Xocn));

%Define boundaries
c2_boundary=initialize_boundaries();
f=2*pi/(24*3600);
Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
min_floe_size = 4*Lx*Ly/10000; %set minimum floe size for initialization

%Initialize Floe state
target_concentration = 1; %Set target concentration for initial conditions
[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,10,min_floe_size);

%create polyshape of any boundary floes
if Nb >0
    polyboundary = union([Floe(1:Nb).poly]); 
else
    polyboundary = [];
end
c2_boundary_poly = polyshape(c2_boundary');
c2_border = polyshape(2*[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly]'); c2_border = subtract(c2_border, c2_boundary_poly);
floebound = initialize_floe_values(c2_border, height);
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
min_floe_size = (4*Lx*Ly-sum(cat(1,Floe(1:Nb).area)))/20000; %set minimum floe size for simulation

%calculate elastic modulus used in floe interactions
global Modulus
Modulus = 1.5e3*(mean(sqrt(cat(1,Floe.area)))+min(sqrt(cat(1,Floe.area))));
save('Modulus.mat','Modulus');

%% Set output parameters 

dhdt = 1; %atmospheric heat flux

nDTOut=150; %Output frequency (in number of time steps)

nSnapshots=50; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nSimp = 20; %Timesteps between simplification of floe boundaries

% Set parallel computing inputs
%nPar = 18; %Number of workers for parfor
%poolobj = gcp('nocreate'); % If no pool, do not create new one.
%if isempty(poolobj)
%    parpool(nPar);
%else
%    delete(poolobj);
%    parpool(nPar);
%end

target_concentration=1; %Set target concentration for when creation of new elements will be run
tStart = tic; 

%ocean forces need updated on less frequent time scale
doInt.flag = true; 
doInt.step = 10; %Timesteps between update of ocean stresses

% specify coarse grid size
LxO= 2*max(ocean.Xo);LyO= 2*max(ocean.Yo);
Nx=10; Ny=10;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);

%Initiailize Eulearian Data
[eulerian_data] = calc_eulerian_data(Floe,Nx,Ny,Nb,c2_boundary,PERIODIC);
Vd = zeros(Ny,Nx,2);
Vdnew=zeros(Ny, Nx);
SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
U = zeros(Ny, Nx); V = zeros(Ny, Nx);
dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
Sig = zeros(Ny, Nx); mass = zeros(Ny,Nx);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
[Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);

%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end
if isempty(dir('Floes')); disp('Creating folder: Floes'); mkdir('Floes'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=figure('Position',[100 100 1000 500],'visible','on');
end

%% Solving for floe trajectories
tic; 
while im_num<nSnapshots

    if mod(i_step,10)==0
        disp(' ');
        toc
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        disp(['number of collisions: ' num2str(numCollisions)]);
        disp(' ');
        tic
        doInt.flag=true;
    else
        doInt.flag=false;
    end

    %Dynamic simplification of floe shape
    if mod(i_step,nSimp)==0
        FloeOld = Floe;
        parfor j=1:length(FloeOld)
            FloeOld(j).poly = polyshape(FloeOld(j).c_alpha'+[FloeOld(j).Xi FloeOld(j).Yi]);
        end
        floenew = [];
        for ii = 1:length(Floe)
            ParFloes(ii).floenew = [];
            ParFloes(ii).kill = [];
            ParFloes(ii).verts = 0;
        end
        if isfield(Floe,'poly')
            Floe=rmfield(Floe,{'poly'});
        end
        parfor ii = 1:length(Floe)
            floe = Floe(ii);
            if length(Floe(ii).c0) > 30
                [floe2,kill] = FloeSimplify(Floe(ii),0,FloeOld,polyboundary);
                if isfield(floe2,'poly')
                    floe2=rmfield(floe2,{'poly'});
                end
                if isempty(kill)
                    ParFloes(ii).kill = [ParFloes(ii).kill kill];
                end
                for jj = 1:length(floe2)
                    if jj == 1
                        Floe(ii) = floe2(jj);
                        ParFloes(ii).verts=length(floe2(jj).c0(1,:));
                    else
                        ParFloes(ii).floenew = [ParFloes(ii).floenew floe2(jj)];
                    end
                end
            end
        end
        clear FloeOld
        if isfield(Floe,'poly')
            Floe=rmfield(Floe,{'poly'});
        end
        floenew =[]; kill = [];
        for ii = 1:length(ParFloes)
            floenew = [floenew ParFloes(ii).floenew];
            kill = [kill ParFloes(ii).kill];
        end
        kill = unique(kill(kill>Nb));
        live = cat(1,Floe.alive);
        live(kill) = 0;
        Floe(live==0)=[]; 
        Floe =[Floe floenew];
    end

    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        

        [eulerian_data] = calc_eulerian_data(Floe,Nx,Ny,Nb,c2_boundary,PERIODIC);
        if ifPlot
            [fig] = plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nb,PERIODIC);
%              saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        end
        
        %Average coarse eulerian data
        if AVERAGE
            SigXXa = SigXX/fix(nDTOut); SigYXa = SigYX/fix(nDTOut);
            SigXYa = SigXY/fix(nDTOut); SigYYa = SigYY/fix(nDTOut);
            Eux = Eux/fix(nDTOut); Evx = Evx/fix(nDTOut);
            Euy = Euy/fix(nDTOut); Evy = Evy/fix(nDTOut);
            U = U/fix(nDTOut); V = V/fix(nDTOut);
            dU = dU/fix(nDTOut); dV = dV/fix(nDTOut);
            Sig = Sig/fix(nDTOut); 
            mass = mass/fix(nDTOut);
        else
            SigXXa = squeeze(eulerian_data.stressxx); SigYXa = squeeze(eulerian_data.stressyx);
            SigXYa = squeeze(eulerian_data.stressxy); SigYYa = squeeze(eulerian_data.stressyy);
            Eux = squeeze(eulerian_data.strainux); Evx = squeeze(eulerian_data.strainvx);
            Euy = squeeze(eulerian_data.strainuy); Evy = squeeze(eulerian_data.strainvy);
            U = U+squeeze(eulerian_data.u);V = V+squeeze(eulerian_data.v);
            dU = dU+squeeze(eulerian_data.du);dV = dV+squeeze(eulerian_data.dv);
            Sig = Sig+squeeze(eulerian_data.stress);
        end
        
                
        if ifPlotStress
            [fig] =plot_basic_stress(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
            figure(2)
            plot(Xc,SigXXa,'kx','linewidth',2); title(['Time = ' num2str(Time/3600) ' hours'],'fontsize',24);
            drawnow
        end
        
        if ifPlotStressStrain
            fig2 = figure(fig2);
            SigO = Sig;
            subplot(2,4,1); imagesc(Xc,Yc,SigXXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{xx}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,2); imagesc(Xc,Yc,SigYXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{yx}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,5); imagesc(Xc,Yc,SigXYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{xy}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,6); imagesc(Xc,Yc,SigYYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{yy}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,3); imagesc(Xc,Yc,Eux); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{11}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,4); imagesc(Xc,Yc,Evx); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{21}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,7); imagesc(Xc,Yc,Euy); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{12}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,8); imagesc(Xc,Yc,Evy); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{22}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            saveas(fig2,['./figs/' num2str(im_num,'Stress%03.f') '.jpg'],'jpg');
        end
        
    end
    
    %Periodic creation of new elements in open ocean
    if PACKING && h0 > 0
        if mod(i_step,nDTpack)==0
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = create_new_ice(Floe,c2_boundary,dhdt,Vd,target_concentration, height, min_floe_size, PERIODIC,3,3,Nb);
        end
    end
    
    %output and reset data values
    if mod(i_step,nDTOut)==0
        save(['./Floes/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','eulerian_data','SigXXa','SigXYa', 'SigYYa','U','dU','mass');
        SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
        SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
        Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
        Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
        U = zeros(Ny, Nx); V = zeros(Ny, Nx);
        dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
        Sig = zeros(Ny, Nx); mass = zeros(Ny, Nx);
        
        M = cat(1,Floe.mass);
        Mtot(im_num) = sum(M)+sum(Vdnew(:));
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING);

    %Average coarse eulerian data
    if AVERAGE
        [eulerian_data] = calc_eulerian_data(Floe,Nx,Ny,Nb,c2_boundary,PERIODIC);
        SigXX = SigXX+squeeze(eulerian_data.stressxx); SigYX = SigYX+squeeze(eulerian_data.stressyx);
        SigXY = SigXY+squeeze(eulerian_data.stressxy); SigYY = SigYY+squeeze(eulerian_data.stressyy);
        Eux = Eux+squeeze(eulerian_data.strainux); Evx = Evx+squeeze(eulerian_data.strainvx);
        Euy = Euy+squeeze(eulerian_data.strainuy); Evy = Evy+squeeze(eulerian_data.strainvy);
        U = U+squeeze(eulerian_data.u);V = V+squeeze(eulerian_data.v);
        dU = dU+squeeze(eulerian_data.du);dV = dV+squeeze(eulerian_data.dv); 
        Sig = Sig+squeeze(eulerian_data.stress);
        mass = mass+squeeze(eulerian_data.Mtot);
    end
    
    % Weld floes
    if WELDING && mod(i_step,25)==0 && dhdt > 0
        Fweld = 150;%Welding probability coefficient
        A=cat(1,Floe.area);
        if max(A) > Amax
           Amax = max(A); 
        end
        if mod(i_step,5000)==0
            Floe = weld(Floe,Nb,Fweld,c2_boundary,Amax/2,1,1);
        elseif mod(i_step,500)==0
            Floe = weld(Floe,Nb,Fweld,c2_boundary,Amax/3,2,2);
        else
            Floe = weld(Floe,Nb,Fweld,c2_boundary,Amax/3,3,3);
        end
    end
    
    % Fracture floes based upon stress
    if FRACTURES && mod(i_step,75)==0 
        compactness = 1; 
        [Floe,Princ] = fracture(Floe,Nb,min_floe_size,compactness);
    end
    
    % Remove sharp corners of interacting floes
    if CORNERS && mod(i_step,10)==0
        
        keep = rand(length(Floe),1)>0.7;
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep(overlapArea>0.15) = 0;
        keep(1:Nb) = ones(Nb,1); %Do not fracture boundaries
        if KEEP_MIN %Do not allow floes below threshold to fracture
            keep(cat(1,Floe.area)<min_floe_size)=1;
        end
        fracturedFloes=corners(Floe(~keep),Nb,Floe,c2_boundary);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
    end
    
    %Advect the dissolved mass
    Area=cat(1,Floe.area);
    if ~KEEP_MIN
        dissolvedNEW = calc_dissolved_mass(Floe(Area<min_floe_size),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
    end
%     Vdnew = Advect_Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt);
    Vdnew = Vd(:,:,1)+dissolvedNEW;
    dissolvedNEW=zeros(Ny,Nx);
    Vd(:,:,2) = Vd(:,:,1);
    Vd(:,:,1) = Vdnew;
    
    % Remove small floes from simulation
    Area=cat(1,Floe.area);
    if ~KEEP_MIN
        if sum(Area<min_floe_size)>0, display(['num of small floes killed:' num2str(sum(Area<min_floe_size))]); end
        Floe=Floe(Area> min_floe_size);
    end
    live = cat(1,Floe.alive);
    Floe(live == 0) = [];
    
    
    Time=Time+dt; i_step=i_step+1; %update time index

    
end

tEnd = toc(tStart)
