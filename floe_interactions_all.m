function [Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, ocean, winds,c2_boundary, dt, HFo, min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING)

id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
global Modulus

Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
live = cat(1,Floe.alive);
Floe(live==0)=[];
FloeNums = 1:length(Floe);

%% Find ghost floes if periodicity is being used
N0=length(Floe);
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    translation = [];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            FloeNums(end+1) = -abs(FloeNums(i));
            parent=[parent  i];
            translation = [translation; -2*Lx*sign(x(i)) 0];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            FloeNums(end+1) = -abs(FloeNums(i));
            parent=[parent  i];
            translation = [translation; 0 -2*Ly*sign(y(i))];

        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%% Find potential interactions for calculating floe interactions
%Find length of new Floe variable including the ghost floes
N=length(Floe);
x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
rmax=cat(1,Floe.rmax);
alive=cat(1,Floe.alive);

for i=1+Nb:N  %do interactions with boundary in a separate parfor loop
    
    Floe(i).interactions=[];
    
    Floe(i).OverlapArea = 0;
    
    Floe(i).potentialInteractions=[];
    
    Floe(i).collision_force=[0 0];
    
    Floe(i).Stress=zeros(2);
    
    Floe(i).collision_torque=0;
    
    k=1;
    
    %check if parents are already interacting
    mems = [];
    if FloeNums(i) < 0
        num = abs(FloeNums(i));
        if ~isempty(Floe(num).potentialInteractions)
            mems = cat(1,Floe(num).potentialInteractions.floeNum);
        end
    end
    
    if ( alive(i) && ~isnan(x(i)) ) && COLLISION
        for j=1:N
            if j>i && alive(j) && sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j)) && (~ismember(abs(FloeNums(j)),mems) || 2*(rmax(i)+rmax(j))>min(2*Lx,2*Ly))% if floes are potentially overlapping
                Floe(i).potentialInteractions(k).floeNum=j;
                Floe(i).potentialInteractions(k).c=[Floe(j).c_alpha(1,:)+x(j); Floe(j).c_alpha(2,:)+y(j)];
                Floe(i).potentialInteractions(k).Ui=Floe(j).Ui;
                Floe(i).potentialInteractions(k).Vi=Floe(j).Vi;
                Floe(i).potentialInteractions(k).h=Floe(j).h;
                Floe(i).potentialInteractions(k).area=Floe(j).area;
                Floe(i).potentialInteractions(k).Xi=x(j);
                Floe(i).potentialInteractions(k).Yi=y(j);
                Floe(i).potentialInteractions(k).ksi_ice = Floe(j).ksi_ice;
                mems = [mems; FloeNums(j)];
                k=k+1;
            end
            
        end
        
    end
end

%% Calculate floe-floe interaction forces
kill = zeros(1,N0); transfer = kill;
%for i=1+Nb:N  %now the interactions could be calculated in a parfor loop!
parfor i=1+Nb:N  %now the interactions could be calculated in a parfor loop!
        
    if ~isempty(Floe(i).potentialInteractions)
        
        for k=1:length(Floe(i).potentialInteractions)
            
            floeNum=Floe(i).potentialInteractions(k).floeNum;
            
            [force_j,P_j, overlap] = floe_interactions(Floe(i),Floe(i).potentialInteractions(k),c2_boundary,PERIODIC,Modulus,dt);
            
            if sum(abs(force_j(:)))~=0
                Floe(i).interactions=[Floe(i).interactions ; floeNum*ones(size(force_j,1),1) force_j P_j zeros(size(force_j,1),1) overlap'];
                Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
            elseif isinf(overlap)
                if i <= N0 && sign(overlap)>0
                    kill(i) = i;
                    transfer(i)=floeNum;
                elseif floeNum <= N0
                    kill(i) = floeNum;
                end
            end
            
        end
        
    end
    if ~PERIODIC
        [force_b, P_j, overlap] = floe_interactions(Floe(i), floebound,c2_boundary,PERIODIC,Modulus,dt);
        in = inpolygon(x(i),y(i),c2_boundary(1,:)',c2_boundary(2,:)');
        if ~in
            Floe(i).alive = 0;
        end
        
        if sum(abs(force_b(:)))~=0
            [mm,~] = size(P_j);
            for ii =1:mm
                if abs(P_j(ii,2)) == Ly
                    force_b(ii,1) = 0;
                end
                if abs(P_j(ii,1)) == Lx
                    force_b(ii,2) = 0;
                end
            end
            Floe(i).interactions=[Floe(i).interactions ; Inf*ones(size(force_b,1),1) force_b P_j zeros(size(force_b,1),1) overlap'];
            Floe(i).OverlapArea = sum(overlap)+Floe(i).OverlapArea;
            Floe(i).potentialInteractions(end+1).floeNum = Inf;
            Floe(i).potentialInteractions(end).c = c2_boundary;
        end
    end
    
end
for i = 1:length(kill)
    if abs(kill(i)-i)>0 && kill(i)>0
        transfer(kill(i))=i;
    end
end


if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end

%% Fill the lower part of the interacton matrix (floe_i,floe_j) for floes with j<i
for i=1:N %this has to be done sequentially
      
    if ~isempty(Floe(i).interactions)
        
        a=Floe(i).interactions;
        
        indx=a(:,1);
        
        for j=1:length(indx)
            
            if indx(j)<=N && indx(j)>i
                Floe(indx(j)).interactions=[Floe(indx(j)).interactions; i -a(j,2:3) a(j,4:5) 0 a(j,7)];   % 0 is torque here that is to be calculated below
                Floe(indx(j)).OverlapArea = Floe(indx(j)).OverlapArea + a(j,7);
                m = size(Floe(indx(j)).potentialInteractions,2);
                Floe(indx(j)).potentialInteractions(m+1).floeNum=i;
                Floe(indx(j)).potentialInteractions(m+1).c=[Floe(i).c_alpha(1,:)+x(i); Floe(i).c_alpha(2,:)+y(i)];
                Floe(indx(j)).potentialInteractions(m+1).Ui=Floe(i).Ui;
                Floe(indx(j)).potentialInteractions(m+1).Vi=Floe(i).Vi;
                Floe(indx(j)).potentialInteractions(m+1).Xi=x(i);
                Floe(indx(j)).potentialInteractions(m+1).Yi=y(i);
                Floe(indx(j)).potentialInteractions(m+1).ksi_ice = Floe(i).ksi_ice;
            end
            
        end

    end

end


%% calculate all torques from forces
if PERIODIC
    
%     for i=N0+1:N %do this in parfor
   parfor i=N0+1:N %do this in parfor
        
        if ~isempty(Floe(i).interactions)
            
            a=Floe(i).interactions;
            r=[x(i) y(i)];
            for k=1:size(a,1)
                floe_Rforce=a(k,4:5);
                floe_force=a(k,2:3);
                floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
                Floe(i).interactions(k,6)=floe_torque(3);
            end
            
            Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1);
            Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1);
            
        end
        
    end
     %add forces and torques from ghost floes to their parents; ghost floes
    %begin with the index N0+1
    for i=1:length(parent)
        Floe(parent(i)).collision_force =Floe(parent(i)).collision_force +Floe(N0+i).collision_force;
        Floe(parent(i)).collision_torque=Floe(parent(i)).collision_torque+Floe(N0+i).collision_torque;
    end
end

%% Calculate updates to floe trajectories
parfor i=1+Nb:N0
    
    if ~isempty(Floe(i).interactions)
        
       a=Floe(i).interactions;
       r=[x(i) y(i)];
        for k=1:size(a,1)
            floe_Rforce=a(k,4:5);
            floe_force=a(k,2:3);
            floe_torque=cross([floe_Rforce-r 0], [floe_force 0]);
            Floe(i).interactions(k,6)=floe_torque(3);
        end
        
       Floe(i).collision_force=sum(Floe(i).interactions(:,2:3),1)+Floe(i).collision_force;
       Floe(i).collision_torque=sum(Floe(i).interactions(:,6),1)+Floe(i).collision_torque;
        
    end
    
    if PERIODIC
            
        if abs(Floe(i).Xi)>Lx %if floe got out of periodic bounds, put it on the other end
            Floe(i).Xi=Floe(i).Xi-2*Lx*sign(Floe(i).Xi);
        end
        
        if abs(Floe(i).Yi)>Ly %if floe got out of periodic bounds, put it on the other end
            Floe(i).Yi=Floe(i).Yi-2*Ly*sign(Floe(i).Yi);
        end
        
    end
    
   %Do the timestepping now that forces and torques are known.
    if Floe(i).alive
        [tmp,Fx,Fy] =calc_trajectory(dt,ocean,winds,Floe(i),HFo,doInt);
        if (isempty(tmp) || isnan(x(i)) ), kill(i)=i; else; Floe(i)=tmp; end
    end

end


%%  Ridging 
floenew = [];
Ridged = zeros(1,length(Floe));
if RIDGING && doInt.flag
    %Create a function to control probability that ridging will occur
    h = cat(1,Floe.h);
     keepR=rand(length(Floe),1)<0.05;
    for ii=1+Nb:N0
        
        if Floe(ii).alive && ~isempty(Floe(ii).interactions)
            a = Floe(ii).interactions;
            c1 = Floe(ii).c_alpha+[Floe(ii).Xi; Floe(ii).Yi];
            abound = zeros(1+Nb,1);
            if ~isempty(a)
                if ~isempty(InterX(c1,c2_boundary))
                    abound(1+Nb) = 1;
                end
                a(isinf(a(:,1)),:)=[];
            end

            if  ~keepR(ii) && h(ii)<5  && ~isempty(a)
                clear overlap;
                for jj = 1:size(a,1)
                    if a(jj,1) < length(Floe)+1
                        overlap(jj) = a(jj,7)/min([Floe(ii).area Floe(a(jj,1)).area]);
                    else
                        overlap(jj) = a(jj,7)/Floe(ii).area;
                    end
                end
                overlap(overlap<1e-6) = nan; overlap(overlap>0.95) = nan;
                overlappingFloes = a(~isnan(overlap),1);
                overlappingFloes = unique(overlappingFloes);
                abound(overlappingFloes<Nb+1) = 1;
                for jj = length(overlappingFloes):-1:1
                    if Ridged(overlappingFloes(jj))
                        overlappingFloes(jj)=[];
                    end
                end
                for jj = 1:length(overlappingFloes)
                    if Floe(overlappingFloes(jj)).h < 5
                        [Floe1, Floe2] = ridge(Floe(ii),Floe(overlappingFloes(jj)),c2_boundary_poly,PERIODIC,min_floe_size);
                        if length(Floe1) > 1
                            Floe(ii) = Floe1(1);
                            Ridged(ii) = 1;
                            floenew = [floenew Floe1(2:end)];
                        else
                            Floe(ii) = Floe1;
                            if Floe1.alive == 0
                                kill(ii) = ii;
                            end
                        end
                        if length(Floe2) > 1
                            Floe(overlappingFloes(jj)) = Floe2(1);
                            Ridged(overlappingFloes(jj)) = 1;
                            floenew = [floenew Floe2(2:end)];
                        else
                            Floe(overlappingFloes(jj)) = Floe2;
                            if Floe2.alive == 0 && overlappingFloes(jj) <= N0
                                kill(overlappingFloes(jj)) = overlappingFloes(jj);
                            end
                        end
                    end
                end

            end
            if sum(abound)>0 && h(ii)<1.25 && Floe(ii).area > min_floe_size
                for jj = 1:length(abound)
                    if abound(jj) == 1 && jj == Nb+1
                        [Floe1, ~] = ridge(Floe(ii),floebound,c2_boundary_poly,PERIODIC,min_floe_size);
                        poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
                    elseif abound(jj)  == 1
                        poly = polyshape([Floe(abound(jj)).c_alpha(1,:)+Floe(abound(jj)).Xi; Floeabound(jj).c_alpha(2,:)+Floe(abound(jj)).Yi]');
                        [Floe1, ~] = ridge(Floe(ii),Floe(abound(jj)),poly,PERIODIC,min_floe_size);
                    end
                    if length(Floe1) > 1
                        Floe(ii) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(ii) = Floe1;
                        if Floe1.alive == 0
                            kill(ii) = ii;
                        end
                    end
                end
            end
        end
    end
end

%% Rafting
Rafted = zeros(1,length(Floe));
if RAFTING && doInt.flag
    %Create a function to control probability that ridging will occur
    h = cat(1,Floe.h);
    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
    keepR=rand(length(Floe),1)>0.5*overlapArea;
    for ii=1+Nb:N0
        
        if Floe(ii).alive && ~isempty(Floe(ii).interactions)
            a = Floe(ii).interactions;
            c1 = Floe(ii).c_alpha+[Floe(ii).Xi; Floe(ii).Yi];
            abound = zeros(1+Nb,1);
            if ~isempty(a)
                if ~isempty(InterX(c1,c2_boundary))
                    abound(1+Nb) = 1;
                end
                a(isinf(a(:,1)),:)=[];
            end

            if  ~keepR(ii) && h(ii)<0.25  && ~isempty(a)
                clear overlap;
                for jj = 1:size(a,1)
                    if a(jj,1) < length(Floe)+1
                        overlap(jj) = a(jj,7)/min([Floe(ii).area Floe(a(jj,1)).area]);
                    else
                        overlap(jj) = a(jj,7)/Floe(ii).area;
                    end
                end
                overlap(overlap<1e-6) = nan; overlap(overlap>0.95) = nan;
                overlappingFloes = a(~isnan(overlap),1);
                overlappingFloes = unique(overlappingFloes);
                abound(overlappingFloes<Nb+1) = 1;
                for jj = length(overlappingFloes):-1:1
                    if Rafted(overlappingFloes(jj))
                        overlappingFloes(jj)=[];
                    end
                end
                for jj = 1:length(overlappingFloes)
                    if Floe(overlappingFloes(jj)).h < 0.25
                        [Floe1, Floe2] = raft(Floe(ii),Floe(overlappingFloes(jj)),c2_boundary_poly,PERIODIC,min_floe_size);
                        
                        if length(Floe1) > 1
                            Floe(ii) = Floe1(1);
                            Rafted(ii) = 1;
                            floenew = [floenew Floe1(2:end)];
                        else
                            Floe(ii) = Floe1;
                            if Floe1.alive == 0
                                kill(ii) = ii;
                            end
                        end
                        if length(Floe2) > 1
                            Floe(overlappingFloes(jj)) = Floe2(1);
                            Rafted(overlappingFloes(jj)) = 1;
                            floenew = [floenew Floe2(2:end)];
                        else
                            Floe(overlappingFloes(jj)) = Floe2;
                            if Floe2.alive == 0 && overlappingFloes(jj) <= N0
                                kill(overlappingFloes(jj)) = overlappingFloes(jj);
                            end
                        end
                    end
                end

            end
            if sum(abound)>0 && h(ii)<0.25 && Floe(ii).area > min_floe_size
                for jj = 1:length(abound)
                    if abound(jj) == 1 && jj == Nb+1
                        [Floe1, ~] = raft(Floe(ii),floebound,c2_boundary_poly,PERIODIC,min_floe_size);
                        poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
                    elseif abound(jj)  == 1
                        poly = polyshape([Floe(abound(jj)).c_alpha(1,:)+Floe(abound(jj)).Xi; Floeabound(jj).c_alpha(2,:)+Floe(abound(jj)).Yi]');
                        [Floe1, ~] = raft(Floe(ii),Floe(abound(jj)),poly,PERIODIC,min_floe_size);
                    end
                    if length(Floe1) > 1
                        Floe(ii) = Floe1(1);
                        floenew = [floenew Floe1(2:end)];
                    else
                        Floe(ii) = Floe1;
                        if Floe1.alive == 0
                            kill(ii) = ii;
                        end
                    end
                end
            end
        end
    end
end

%% Remove any floes that were created for computations or lost from interactions
Floe=Floe(1:N0); % ditch the ghost floes.

if ~isempty(kill(kill>0)) 
    kill(kill>N0)=0;
    transfer(transfer>N0) = 0;
    transfer = transfer(kill>0);
    kill = kill(kill>0);
    for ii = 1:length(kill)
        if Floe(kill(ii)).alive>0 
            if transfer(ii)>0 && (Floe(kill(ii)).area>2e4 || Floe(kill(ii)).area>2e4)
                Floe1 = Floe(kill(ii));
                Floe2 = Floe(transfer(ii));
                Floe1.poly = polyshape(Floe1.c_alpha'+[Floe1.Xi Floe1.Yi]);
                Floe2.poly = polyshape(Floe2.c_alpha'+[Floe2.Xi Floe2.Yi]);
                floes  = Fuse_Floes(Floe1,Floe2);
                if isfield(floes,'poly')
                    floes=rmfield(floes,{'poly'});
                end
                if length(floes)>1
                    Floe(transfer(ii)) = floes(1);
                    for jj = 2:length(floes)
                        floenew = [floenew floes(2:end)];
                    end
                elseif isempty(floes)
                    save('FloeTransfer.mat','Floe','kill','transfer','Floe1','Floe2')
                else
                    Floe(transfer(ii)) = floes(1);
                end
            end
            dissolvedNEW = dissolvedNEW+calc_dissolved_mass(Floe(kill(ii)),Nx,Ny,c2_boundary_poly);
            Floe(kill(ii)).alive = 0;
        end
    end
end
Floe = [Floe floenew];
live = cat(1,Floe.alive);
Floe(live==0)=[]; %remove any floes that got dissolved so they do not take up space


Floe=rmfield(Floe,{'potentialInteractions'});


warning('on',id)
warning('on',id3)
end
