function [Floe,Vd] = create_new_ice(Floe,c2_boundary,dhdt,Vd,target, height, min_floe_size, PERIODIC,Nx,Ny,Nb)
%% This function takes in the existing floe state and creates new thin floes in the open space to match the input target concentration
id = 'MATLAB:polyshape:tinyBoundaryDropped';
warning('off',id);
id2 ='MATLAB:polyshape:repairedBySimplify';
warning('off',id2)

alive = cat(1,Floe.alive);
Floe(~logical(alive)) = [];
Fold = Floe;
for kk = 1:length(Floe)
    Fold(kk).poly = polyshape(Fold(kk).c_alpha'+[Fold(kk).Xi Fold(kk).Yi]);
end

N0=length(Floe);
Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));

%If the floe is periodic populate the ghost floes
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
            parent=[parent  i];
            translation = [translation; 0 -2*Ly*sign(y(i))];
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%Caclulate the coarse grid and any potential interactions for the different
%regions of the floe where ice needs to be created
floenew = [];
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
[xx,yy] = meshgrid(0.5*(x(1:end-1)+x(2:end)),0.5*(y(1:end-1)+y(2:end)));
xf = cat(1,Floe.Xi);
yf = cat(1,Floe.Yi);
r_max = sqrt((dx/2)^2+(dy/2)^2);
rmax = cat(1,Floe.rmax);
potentialInteractions = zeros(Ny,Nx,length(Floe));

%Identfy which floes are potentially within the subgrids
for kk = 1:length(Floe)
    Floe(kk).poly = polyshape(Floe(kk).c_alpha'+[Floe(kk).Xi Floe(kk).Yi]);
    pint = sqrt((xx-xf(kk)).^2+(yy-yf(kk)).^2)-(rmax(kk)+r_max);
    pint(pint>0) = 0;
    pint(pint<0) = 1;
    potentialInteractions(:,:,kk) = pint;
end

%function to determine the probablity that ice gets created
SimpMin = @(A) 3*log10(A);

%% Loop through all regions creating sea ice in each if probability criteria is met
for j = 1:Nx*Ny
    FloePack(j).floenew = [];
    kill(j).floes = [];
    floe2(j).floe = [];
end
parfor j = 1:Nx*Ny 

    count = 1;
    
    ii = ceil(j/Ny);
    jj = mod(j,Ny);
    if jj == 0;  jj = Ny; end
    
    %Find coverage of sea ice in region of interest and calculate concentration
    bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
    box = polyshape(bound(1,:), bound(2,:));
    k = find(logical(potentialInteractions(jj,ii,:))==1);
    if isempty(k)
        poly = [];
        A = 0;
        inBox = [];
    else
        poly = intersect(box,[Floe(logical(potentialInteractions(jj,ii,:))).poly]);
        A = area(poly);
        inBox = k(A>0);
        poly = [Floe(inBox).poly];
        Atot = area(poly);
    end
    
    c = sum(A)/area(box);
    
    %if concentration is below target then create new sea ice
    if c<0.999*target
        atarget = target*area(box);
        
        %Break up open water into chunks that will be the new floes
        [Xi,Yi] = centroid(box);
        N = ceil(atarget/(50*min_floe_size));
        if N < 3
            N = 3;
        elseif N > 5
            N = 5;
        end
        
        X = Xi+r_max*(2*rand(N,1)-1);
        Y = Yi+r_max*(2*rand(N,1)-1);
        boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*r_max+ [Xi Yi];
        [~, b,~,~,~] = polybnd_voronoi([X Y],boundingbox);
        
        
        subfloes = [];
        for kk = 1:length(b)
            new = subtract(polyshape(b{kk}),poly);
            new = intersect(new);
            new = intersect(new,box);
            for iii = 1:length(new)
                R = regions(new(iii));
                subfloes = [subfloes; R(area(R)>min_floe_size)];
            end
        end
        
        %Add in new floes until we reach target concentration
        for iii = 1:length(subfloes)
            
            polyFull = rmholes(subfloes(iii));
            for jjj = 1:length(polyFull)
                if subfloes(iii).NumHoles > 0
                    hnew.mean = area(subfloes(iii))*height.mean/area(polyFull(jjj));
                    hnew.delta = height.delta;
                    floe2(j).floe = initialize_floe_values(polyFull(jjj),hnew); %create new floe
                else
                    floe2(j).floe = initialize_floe_values(polyFull(jjj),height); %create new floe
                end
                
                %floes can not have holes in them so find which floes from
                %simulation are in those holes and weld together
                if subfloes(iii).NumHoles > 0
                    poly2 = rmholes(floe2(j).floe.poly);
                    live = cat(1,Floe.alive);
                    if sum(area(intersect(poly2,[Floe(live==1).poly])))>0
                        polyO = intersect(poly2,[Floe(inBox).poly]);
                        A2 = area(polyO);
                        polyin = A2./Atot;
                        in = inBox(polyin>0.99);
                        in = flipud(in);
                        boundaryFloe = [];
                        for kk = 1:length(in)
                            if in(kk) > Nb
                                kill(j).floes = [kill(j).floes; in];
                            else
                                boundaryFloe = [boundaryFloe; in(kk)];
                            end
                        end
                        in(in<Nb+1) = [];
                        %if it is a boundary floe that is causing hole then
                        %split floe through the middle of that boundary
                        if ~isempty(boundaryFloe)
                            subnew = [];
                            for kk = 1:length(boundaryFloe)
                                Xi2 = Floe(boundaryFloe(kk)).Xi; Yi2 = Floe(boundaryFloe(kk)).Yi;
                                L = [ Xi2 Yi2; Xi2+1 Yi2];
                                PcT = cutpolygon(polyFull(jjj).Vertices, L, 1);
                                pnew = polyshape(PcT);
                                pnew = subtract(pnew,union(Floe(boundaryFloe).poly));
                                if ~isempty(subnew); pnew = subtract(pnew,union(subnew)); end
                                if kk==length(boundaryFloe)
                                    PcB = cutpolygon(polyFull.Vertices, L, 2);
                                    pnew2 = polyshape(PcB);
                                    pnew2 = subtract(pnew2,union(Floe(boundaryFloe).poly));
                                    subnew = [subnew; regions(pnew); regions(pnew2)];
                                else
                                    subnew = [subnew regions(pnew)];
                                end
                            end
                        else
                            subnew = polyFull(jjj);
                        end
                        newfloes = [];
                        for kk = 1:length(subnew)
                            tmp = initialize_floe_values(subnew(kk),hnew);
                            newfloes = [newfloes tmp];
                        end
                        floe2(j).floe = newfloes;
                        for kk = 1:length(in)
                            overlap = area(intersect(subnew,Floe(in(kk)).poly))/Floe(in(kk)).area;
                            [~,OverlappingFloe] = max(overlap);
                            tmpFloe = Fuse_Floes(floe2(j).floe(OverlappingFloe),Floe(in(kk)));
                            if length(tmpFloe) > 1
                                [~,I] = max(cat(1,floe2(j).floe.area));
                                tmp = floe2(j).floe(I);
                                floe2(j).floe(I) = [];
                                FloePack(j).floenew = [FloePack(j).floenew floe2(j).floe];
                                floe2(j).floe(OverlappingFloe) = tmp;
                            else
                                floe2(j).floe(OverlappingFloe) = tmpFloe;
                            end
                            count = count+1;
                        end
                        if length(floe2(j).floe) > 1
                            [~,I] = max(cat(1,floe2(j).floe.area));
                            tmp = floe2(j).floe(I);
                            floe2(j).floe(I) = [];
                            FloePack(j).floenew = [FloePack(j).floenew floe2(j).floe];
                            floe2(j).floe = tmp;
                        end
                    end
                end
                if subfloes(iii).NumHoles > 0
                    hnew = height;
                    hnew.mean = floe2(j).floe.h;
                end
                
                FloePack(j).floenew = [FloePack(j).floenew floe2(j).floe];
                floe2(j).floe = [];
            end
        end

    end
    
end

%Remove any floes that were determined to need to be removed from
%simulation
Floe=Floe(1:N0);
killFloes = [];
floenew = [];
for ii = 1:length(kill)
    killFloes = [killFloes; kill(ii).floes];
    floenew = [floenew FloePack(ii).floenew];
end
if ~isempty(killFloes)
    killFloes = unique(killFloes(killFloes>0));
    killFloes = sort(killFloes,'descend');
    for ii = 1:length(killFloes)
        Floe(killFloes(ii)).alive = 0;
    end
end
alive = cat(1,Floe.alive);
Floe(alive==0) = [];
if isfield(floenew,'potentialInteractions')
    floenew=rmfield(floenew,{'potentialInteractions'});
end
if isfield(Floe,'potentialInteractions')
    Floe=rmfield(Floe,{'potentialInteractions'});
end

Floe = [Floe floenew];
warning('on',id)
warning('on',id2)

Floe=rmfield(Floe,{'poly'});
Vd(Vd<0)=0;

end
