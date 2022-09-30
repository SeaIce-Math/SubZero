function [Floe] = weld(Floe,Nb,Fweld,c2_boundary,maxWeld,Nx,Ny)
%%This function takes in two floe field and based upon a specified
%%probability function welds interacting floes together

id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id2 = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id2)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

for ii = 1:length(Floe)
    if abs(Floe(ii).area/area(polyshape(Floe(ii).c_alpha'))-1)>1e-3
        Floe(ii).area = area(polyshape(Floe(ii).c_alpha'));
    end
end

A = cat(1,Floe.area);
Atotal = sum(A);
if Nb > 0
    Fbound = Floe(1:Nb);
else
    Fbound = [];
end
Floe = Floe(1+Nb:length(Floe));

floenew = [];

%Find floes and create bins 
%Create coarse grid and coarse floe variables
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xi=cat(1,Floe.Xi);
Yi=cat(1,Floe.Yi);
Binx = fix((Xi-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Yi-min(y))/(max(y)-min(y))*Ny+1);
Binox = fix(abs(Xi)/max(x));
Binoy = fix(abs(Yi)/max(y));

count = 1;
for i = 1:Nx
    for j = 1:Ny
        liveX = zeros(1,length(Binx)); liveX(Binx == i) = 1;
        liveY = zeros(1,length(Biny)); liveY(Biny == j) = 1;
        FloeBin(count).Floe = Floe(logical( liveX+liveY==2));
        count = count+1;
    end
end
Bino = logical(Binox+Binoy);
FloeOut = [];
if sum(Bino)>0
    FloeOut = Floe(Bino);
end

%Idenfity potential floes to weld with
parfor ii = 1:length(FloeBin)
    FloeBin(ii).floenew = [];
    Xi=cat(1,FloeBin(ii).Floe.Xi);
    Yi=cat(1,FloeBin(ii).Floe.Yi);
    alive = cat(1,FloeBin(ii).Floe.alive);
    rmax = cat(1,FloeBin(ii).Floe.rmax);
    for i =1:length(FloeBin(ii).Floe)
        k=1;
        FloeBin(ii).Floe(i).potentialInteractions = [];
        FloeBin(ii).Floe(i).poly = polyshape(FloeBin(ii).Floe(i).c_alpha'+[FloeBin(ii).Floe(i).Xi FloeBin(ii).Floe(i).Yi]);
        for j=1:length(FloeBin(ii).Floe)
            if alive(j) && sqrt((Xi(i)-Xi(j))^2 + (Yi(i)-Yi(j))^2)>1 && sqrt((Xi(i)-Xi(j))^2 + (Yi(i)-Yi(j))^2)<(rmax(i)+rmax(j)) % if floes are potentially overlapping
                FloeBin(ii).Floe(i).potentialInteractions(k).floeNum=j;
                FloeBin(ii).Floe(i).potentialInteractions(k).c=[FloeBin(ii).Floe(j).c_alpha(1,:)+Xi(j); FloeBin(ii).Floe(j).c_alpha(2,:)+Yi(j)];
                FloeBin(ii).Floe(i).potentialInteractions(k).Ui=FloeBin(ii).Floe(j).Ui;
                FloeBin(ii).Floe(i).potentialInteractions(k).Vi=FloeBin(ii).Floe(j).Vi;
                FloeBin(ii).Floe(i).potentialInteractions(k).h=FloeBin(ii).Floe(j).h;
                FloeBin(ii).Floe(i).potentialInteractions(k).area=FloeBin(ii).Floe(j).area;
                FloeBin(ii).Floe(i).potentialInteractions(k).Xi=Xi(j);
                FloeBin(ii).Floe(i).potentialInteractions(k).Yi=Yi(j);
                FloeBin(ii).Floe(i).potentialInteractions(k).ksi_ice = FloeBin(ii).Floe(j).ksi_ice;
                k=k+1;
            end
        end
    end
end

%Set probability function
hvsd = @(x) [0.5*(x == 0) + (x > 0)];
ramp = @(frac) hvsd(frac)*frac;
SimpMin = @(A) 3*log10(A);

%Loop through all floes to determine which will weld
clear k
clear j
for kk = 1:length(FloeBin)
    Floe = FloeBin(kk).Floe;
    for i = 1:length(Floe)
        if Floe(i).alive && ~isempty(Floe(i).potentialInteractions) && Floe(i).area <maxWeld
            %Identify potential floes to weld with
            c = unique(cat(1,Floe(i).potentialInteractions.floeNum));
            c = c(c<length(Floe));
            c = c(c>i);
            A = cat(1,Floe(i).potentialInteractions.area);
            if ~isempty(c)
                %Find overlapping area between floes
                OverlapA = area(intersect(Floe(i).poly,[Floe(c).poly]));
                OverlapA(A>maxWeld) = 0; %Weld floes can only take up a certain part of computational domain so exclude any that exceed this threshold
                if max(OverlapA)>0
                    live = 0;
                    %Find if any floes meet criteria to weld together
                    p = rand(1,length(OverlapA));
                    weldp = Fweld*OverlapA/Floe(i).area;
                    weld = max(weldp(weldp>p)); 
                    if ~isempty(weld)
                        [~,k]=min(abs(weldp-weld));
                        j = c(k);
                        live = Floe(j).alive;
                    else
                        j = [];
                    end
                    if live && ~isempty(j)
                        if area(union(Floe(i).poly,Floe(j).poly)) < Atotal/5 && area(union(Floe(i).poly,Floe(j).poly)) > 2e4                            
                            floe = Fuse_Floes(Floe(i),Floe(j)); %If probability is met then weld them together
                            if length(floe)>1
                                FloeBin(kk).floenew =[FloeBin(kk).floenew floe(2:length(floe))];
                                floe = floe(1);
                            end
                            
                            %Set one of the two floes to be removed from
                            %the simulation. Also see if this welding
                            %together needs to subsume any other floes.
                            FloeNums1 = cat(1,Floe(i).potentialInteractions.floeNum);
                            FloeNums1(FloeNums1 == j) = [];
                            FloeNums2 = cat(1,Floe(j).potentialInteractions.floeNum);
                            FloeNums2(FloeNums2 == i) = [];
                            PI = unique([FloeNums1; FloeNums2]);
                            PI = PI(PI<length(Floe));
                            for ii = 1:length(PI)
                                if Floe(PI(ii)).alive
                                    Anew = 0;
                                    c = [Floe(PI(ii)).c_alpha(1,:)+Floe(PI(ii)).Xi; Floe(PI(ii)).c_alpha(2,:)+Floe(PI(ii)).Yi]';
                                    [Xnew,Ynew] = polyclip([floe.poly.Vertices],c,'int');
                                    if ~isempty(Xnew)
                                        Xt = Xnew{1}; Yt = Ynew{1};
                                        Anew = polyarea(Xt,Yt);
                                    end
                                    if Anew/Floe(PI(ii)).area>0.4
                                        floe = Fuse_Floes(floe,Floe(PI(ii)));
                                        if length(floe)>1
                                            FloeBin(kk).floenew =[FloeBin(kk).floenew floe(2:length(floe))];
                                            floe = floe(1);
                                        end
                                        Floe(PI(ii)).alive = 0;
                                    end
                                end
                            end
                            
                            Floe(j).alive = 0;
                            for jj = 1:length(floe)
                                if jj == 1
                                    Floe(i) = floe(jj);
                                else
                                    FloeBin(kk).floenew = [FloeBin(kk).floenew floe(jj)];
                                end
                            end
                            
                        end
                    end
                    
                end
            end
        end
    end
    FloeBin(kk).Floe = Floe;
end
%Remove floes designated to be booted from the simulation
Floe = Fbound;
floenew = [];
for ii = 1:length(FloeBin)
    if isfield(FloeBin(ii).Floe,'potentialInteractions')
        FloeBin(ii).Floe=rmfield(FloeBin(ii).Floe,{'potentialInteractions'});
    end
    if isfield(FloeBin(ii).Floe,'poly')
        FloeBin(ii).Floe=rmfield(FloeBin(ii).Floe,{'poly'});
    end
    Floe = [Floe FloeBin(ii).Floe];
    floenew = [floenew FloeBin(ii).floenew];
end

for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end
for ii =1:length(FloeOut)
    FloeOut(ii).poly = polyshape(FloeOut(ii).c_alpha'+[FloeOut(ii).Xi FloeOut(ii).Yi]);
end
if isfield(FloeOut,'potentialInteractions')
    FloeOut=rmfield(FloeOut,{'potentialInteractions'});
end
if isfield(Floe,'potentialInteractions')
    Floe=rmfield(Floe,{'potentialInteractions'});
end
if isfield(floenew,'potentialInteractions')
    floenew=rmfield(floenew,{'potentialInteractions'});
end

Floe = [Floe floenew FloeOut];
live = cat(1,Floe.alive);
Floe(live == 0) = [];

Floe=rmfield(Floe,{'poly'});

warning('on',id)
warning('on',id2)
warning('on',id3)

end

