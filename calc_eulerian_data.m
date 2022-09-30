function [eulerian_data] = calc_eulerian_data(Floe,Nx,Ny,Nb,c2_boundary,PERIODIC)
%% Function to take information of all floes and average them over a corase grained area
id = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id)

%Identify only the live floes
live = cat(1,Floe.alive);
Floe(live==0)=[];

%Remove any boundary floes so they do not count towards floe averages
if Nb > 0
    for ii =1:length(Floe)
        poly(ii) = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
    end
    boundaries = poly(1);
    if Nb>1
        for ii = 2:Nb
            boundaries = union(boundaries,poly(ii));
        end
    end
    clear poly
else
    boundaries = [];
end
Floe(1:Nb) = []; 

%Create ghost floes for periodic floe states
Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));%c2 must be symmetric around x=0 for channel boundary conditions.
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            
        end
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%Create coarse grid and coarse floe variables
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
[xx,yy] = meshgrid(0.5*(x(1:end-1)+x(2:end)),0.5*(y(1:end-1)+y(2:end)));
xf = cat(1,Floe.Xi);
yf = cat(1,Floe.Yi);
rmax = cat(1,Floe.rmax);
dx = abs(x(2)-x(1));

dy = abs(y(2)-y(1));
r_max = sqrt((dx/2)^2+(dy/2)^2);
eulerian_data.u = zeros(Ny,Nx);
eulerian_data.v = zeros(Ny,Nx);
eulerian_data.du = zeros(Ny,Nx);
eulerian_data.dv = zeros(Ny,Nx);
eulerian_data.stress = zeros(Ny,Nx);
eulerian_data.stressxx = zeros(Ny,Nx);
eulerian_data.stressyx = zeros(Ny,Nx);
eulerian_data.stressxy = zeros(Ny,Nx);
eulerian_data.stressyy = zeros(Ny,Nx);
eulerian_data.strainux = zeros(Ny,Nx);
eulerian_data.strainvx = zeros(Ny,Nx);
eulerian_data.strainuy = zeros(Ny,Nx);
eulerian_data.strainvy = zeros(Ny,Nx);
eulerian_data.c = zeros(Ny,Nx);
eulerian_data.Over = zeros(Ny,Nx);
eulerian_data.Mtot = zeros(Ny,Nx);
eulerian_data.area = zeros(Ny,Nx);
eulerian_data.h = zeros(Ny,Nx);


mass = cat(1,Floe.mass);
mass(isnan(mass)==1)=0;
A = cat(1,Floe.area);
Overlap = cat(1,Floe.OverlapArea);
A(isnan(A)==1)=0;
U = cat(1,Floe.Ui);
U(isnan(U)==1)=0;
V = cat(1,Floe.Vi);
V(isnan(V)==1)=0;
H = cat(1,Floe.h);
H(isnan(H)==1)=0;
dU = cat(1,Floe.dUi_p);
dU(isnan(dU)==1)=0;
dV = cat(1,Floe.dVi_p);
dV(isnan(dV)==1)=0;
Stress = zeros(2,2,length(Floe));
Strain = zeros(2,2,length(Floe));
for ii = 1:length(Floe)
    Stress(:,:,ii) = Floe(ii).Stress;
    Strain(:,:,ii) = Floe(ii).strain;
    poly(ii) = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

%Find potential interactions
potentialInteractions = zeros(Ny,Nx,length(Floe));
for ii = 1:length(Floe)
    pint = sqrt((xx-xf(ii)).^2+(yy-yf(ii)).^2)-(rmax(ii)+r_max);
    pint(pint>0) = 0;
    pint(pint<0) = 1;
    potentialInteractions(:,:,ii) = pint;
end


%% Loop to find coarse averages

for ii = 1:Nx
    for jj = 1:Ny
        live = logical(squeeze(potentialInteractions(jj,ii,:)));
        Mtot = sum(cat(1,Floe(live == 1).mass));
        if Mtot>0 
            
            %create the box and identify and subtract off any boundaries
            bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
            box = polyshape(bound(1,:), bound(2,:));
            if ~isempty(boundaries)
                box = subtract(box,boundaries);
            end
            
            %Find all floes from the potentially interacting ones that have
            %a piece in this area. Values weighted by mass in cell
            FloeNums = 1:length(Floe);
            FloeNums(~live) = [];
            overlap = intersect(box,poly(FloeNums));
            Aover = area(overlap)';
            FloeNums(Aover==0)=[];
            Aover(Aover==0) = [];
            A2 = A(FloeNums);
            Mtot = sum(mass(FloeNums).*Aover./A2);
            Atot = sum(Aover);
            eulerian_data.c(jj,ii) = Atot/area(box);
            if Mtot>0
                eulerian_data.Over(jj,ii) = sum(Overlap(FloeNums))/length(Aover);
                eulerian_data.Mtot(jj,ii) = Mtot;
                eulerian_data.area(jj,ii) = Atot;
                eulerian_data.h(jj,ii) = sum(H(FloeNums).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.u(jj,ii) = sum(U(FloeNums).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.v(jj,ii) = sum(V(FloeNums).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.du(jj,ii) = sum(dU(FloeNums).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.dv(jj,ii) = sum(dV(FloeNums).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.stressxx(jj,ii) = sum(squeeze(Stress(1,1,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.stressyx(jj,ii) = sum(squeeze(Stress(1,2,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.stressxy(jj,ii) = sum(squeeze(Stress(2,1,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.stressyy(jj,ii) = sum(squeeze(Stress(2,2,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.strainux(jj,ii) = sum(squeeze(Strain(1,1,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.strainvx(jj,ii) = sum(squeeze(Strain(1,2,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.strainuy(jj,ii) = sum(squeeze(Strain(2,1,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.strainvy(jj,ii) = sum(squeeze(Strain(2,2,(FloeNums))).*mass(FloeNums).*Aover./A2)./Mtot;
                eulerian_data.stress(jj,ii) = max(eig([eulerian_data.stressxx(jj,ii) eulerian_data.stressyx(jj,ii); eulerian_data.stressxy(jj,ii) eulerian_data.stressyy(jj,ii)]));
                if abs(eulerian_data.stress(jj,ii))> 1e8
                    eulerian_data.stress(jj,ii) = 0;
                end
            end
        end
    end
end

warning('on',id)

end

