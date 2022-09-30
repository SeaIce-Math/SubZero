function [floes,Floe0] = FloeSimplify(floe,PACK,Floe0,polyboundary)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
%% Remap the main polygon to a shape with fewer vertices
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id2 = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id2)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
floes = [];
rho_ice = 920;

%Identify potential interactions
x=cat(1,Floe0.Xi);
y=cat(1,Floe0.Yi);
rmax=cat(1,Floe0.rmax);
alive=cat(1,Floe0.alive);
k=1;
floe.potentialInteractions = [];
for j=1:length(Floe0)
    if alive(j) && sqrt((floe.Xi-x(j))^2 + (floe.Yi-y(j))^2) > 1 && sqrt((floe.Xi-x(j))^2 + (floe.Yi-y(j))^2)<(floe.rmax+rmax(j)) % if floes are potentially overlapping
        floe.potentialInteractions(k).floeNum=j;
        floe.potentialInteractions(k).c=[Floe0(j).c_alpha(1,:)+x(j); Floe0(j).c_alpha(2,:)+y(j)];
        floe.potentialInteractions(k).Ui=Floe0(j).Ui;
        floe.potentialInteractions(k).Vi=Floe0(j).Vi;
        floe.potentialInteractions(k).h=Floe0(j).h;
        floe.potentialInteractions(k).area=Floe0(j).area;
        floe.potentialInteractions(k).Xi=x(j);
        floe.potentialInteractions(k).Yi=y(j);
        floe.potentialInteractions(k).ksi_ice = Floe0(j).ksi_ice;
        k=k+1;
    end
end


%Create new simplified polyshape
floenew = floe;
floenew.poly = polyshape(floenew.c_alpha'+[floenew.Xi floenew.Yi]);
[verts] = reducepoly(floe.c0');
pnew = polyshape(verts);
if ~isempty(polyboundary)
    pnew = translate(pnew,[floenew.Xi floenew.Yi]);
    pnew = subtract(pnew,polyboundary); 
    pnew = translate(pnew,-[floenew.Xi floenew.Yi]);
end
Atot = sum(area(pnew));
if isinf(Atot)
  save('floefail.mat','floe','pnew');
elseif isnan(Atot)
  save('floefail.mat','floe','pnew');
elseif Atot ==0
  save('floefail.mat','floe','pnew');
  R = [];
else
  polynew = scale(pnew,sqrt(floe.area/Atot));

  %Align center of old polygon with the enw one
  [x1,y1] = centroid(polynew);
  dx = floe.Xi-x1;
  dy = floe.Yi-y1;

  %Check if simplifcation led to polygon having multiple regions
  polyout = sortregions(polynew,'area','descend');
  R = regions(polyout);
  R = R(area(R)>1e4);
  Atot = sum(area(R));
end

%Check if simplification removed any holes and now needs to weld with smaller
%floes
for jj = 1:length(R)
    for ii = 1:length(floe.potentialInteractions)
        if ~isinf(floe.potentialInteractions(ii).floeNum) & ~isempty(R)
            polyq = rmholes(R(jj));
            [Xi, Yi] = centroid(polyq);
            Xt = Xi+floe.Xi; Yt = Yi+floe.Yi;
            polyq = translate(polyq,[Xt,Yt]);
            poly = polyshape(floe.potentialInteractions(ii).c');
            polyI = intersect(poly,polyq);
            if area(polyI)/area(poly)>0.4
                floenew = Fuse_Floes(floenew,Floe0(floe.potentialInteractions(ii).floeNum));
                Floe0(floe.potentialInteractions(ii).floeNum).alive = 0;
                if length(floenew)>1
                    m = cat(1,floenew.mass); u = cat(1,floenew.Ui); v = cat(1,floenew.Vi);
                    dksi = cat(1,floenew.dksi_ice_p); inertia_moment = cat(1,floenew.inertia_moment);
                    [~,I] = max(m);
                    height.mean = sum(m)/(area(floenew(I).poly)*rho_ice);
                    height.delta = 0;
                    tmp = initialize_floe_values(floenew(I).poly, height);
                    tmp.mass = sum(m); tmp.Ui = sum(m.*u)/tmp.mass; tmp.Vi = sum(m.*v)/tmp.mass;
                    tmp.dksi_ice_p = sum(dksi.*inertia_moment)/tmp.inertia_moment;
                    tmp.FxOA = [];
                    tmp.FyOA = [];
                    floenew = tmp;

                end
            end
        end
    end
end


%% Calculate the new properties associated with this floe since it has a new shape
if length(R) == 1
    floes = floenew;
    R = rmholes(R);
    floes.area = area(R);
    [Xi, Yi] = centroid(R);
    floes.Xi = Xi+floenew.Xi; floes.Yi = Yi+floenew.Yi;
    floes.h = floes.mass/(rho_ice*floes.area);
    floes.c0 = [R.Vertices(:,1)-Xi,R.Vertices(:,2)-Yi; R.Vertices(1,1)-Xi,R.Vertices(1,2)-Yi]';
    floes.angles = polyangles(R.Vertices(:,1),R.Vertices(:,2));
    A_rot=[cos(floes.alpha_i) -sin(floes.alpha_i); sin(floes.alpha_i) cos(floes.alpha_i)]; %rotation matrix
    floes.c_alpha=A_rot*floes.c0;
    floes.inertia_moment = PolygonMoments(floes.c_alpha',floes.h);
    floes.rmax = max(sqrt(floes.c_alpha(1,:).^2+floes.c_alpha(2,:).^2));
    floes.X = floes.rmax*(2*rand(1000,1) - 1);
    floes.Y = floes.rmax*(2*rand(1000,1) - 1);
    floes.A = inpolygon(floes.X,floes.Y,floes.c_alpha(1,:),floes.c_alpha(2,:));
elseif length(R) > 1
    for ii = 1:length(R)
        floenew = floe;
        poly1new = R(ii);
        polya = rmholes(poly1new);
        if ~PACK
            poly1new = polya;
        end
        [Xi,Yi] = centroid(poly1new);
        floenew.area = area(poly1new);
        poly1new = translate(poly1new,[floe.Xi,floe.Yi]);
        floenew.poly = poly1new;
        floenew.mass = floe.mass*area(R(ii))/Atot;
        floenew.h = floenew.mass/(rho_ice*floenew.area);
        floenew.c_alpha = [(polya.Vertices-[Xi Yi])' [polya.Vertices(1,1)-Xi; polya.Vertices(1,2)-Yi]];
        floenew.angles = polyangles(polya.Vertices(:,1),polya.Vertices(:,2));
        floenew.c0 = floenew.c_alpha;
        floenew.inertia_moment = PolygonMoments(floenew.c0',floenew.h);
        floenew.rmax = max(sqrt(floenew.c_alpha(1,:).^2+floenew.c_alpha(2,:).^2));
        floenew.X = floenew.rmax*(2*rand(1000,1) - 1);
        floenew.Y = floenew.rmax*(2*rand(1000,1) - 1);
        floenew.A = inpolygon(floenew.X,floenew.Y,floenew.c_alpha(1,:),floenew.c_alpha(2,:));
        
        floenew.Xi = floe.Xi+Xi; floenew.Yi = floe.Yi+Yi; floenew.alive = 1;
        floenew.alpha_i = 0; floenew.Ui = floe.Ui; floenew.Vi = floe.Vi;
        floenew.dXi_p = floe.dXi_p; floenew.dYi_p = floe.dYi_p;
        floenew.dUi_p = floe.dUi_p; floenew.dVi_p = floe.dVi_p;
        floenew.dalpha_i_p = 0; floenew.ksi_ice = floenew.area/floe.area*floe.ksi_ice;
        floenew.dksi_ice_p = floe.dksi_ice_p;
        floenew.interactions = [];
        floenew.potentialInteractions = [];
        floenew.collision_force = 0;
        floenew.collision_torque = 0;
        floenew.OverlapArea = 0;
        
        floes = [floes floenew];
    end
else
    floes = floe;
    floes.alive = 0;
end
if isfield(floes,'potentialInteractions')
    floes=rmfield(floes,{'potentialInteractions'});
end

warning('on',id)
warning('on',id2)
warning('on',id3)
end

