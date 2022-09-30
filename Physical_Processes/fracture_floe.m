function Floes=fracture_floe(Floe,N,Floe0)
%%This function takes an input of floes and fractures each floe into a
%%specified number of smaller ones using Voronoi Tesselation
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Floes=[]; rho_ice = 920; 
[~,~,Nz] = size(Floe0(1).StressH);

for kk=1:length(Floe)
    
    %Slight permanent deformation to floe upon fracture
    floe = Floe(kk);
    a = floe.interactions;
    XInew = [];
    if ~isempty(a)
        k1 = isinf(a(:,1));
        a(k1,:) = [];
    end
    if ~isempty(a)
        [~,k] = max(a(:,7));
        X1new = floe.c_alpha(1,:)+floe.Xi;
        Y1new = floe.c_alpha(2,:)+floe.Yi;
        clip = a(k,1);
        if clip < length(Floe0)+1
            X2new = Floe0(clip).c_alpha(1,:)+Floe0(clip).Xi;
            Y2new = Floe0(clip).c_alpha(2,:)+Floe0(clip).Yi;
            [XInew,YInew] = polyclip([X1new' Y1new'],[X2new' Y2new'],'int');
        end
    end
    if ~isempty(XInew)
        Xt = XInew{1}; Yt = YInew{1};
        [xm,ym] = centroid(polyshape(Xt,Yt));
        [d_min1] = p_poly_dist(xm, ym,Xt', Yt');
        F = vecnorm(a(k,2:3));
        xs = a(k,2)*abs(d_min1)/2/F;
        ys = a(k,3)*abs(d_min1)/2/F;
        X2new = X2new + xs; Y2new = Y2new+ys;
        [XInew,YInew] = polyclip([X1new' Y1new'],[X2new' Y2new'],'dif');
        if ~isempty(XInew)
            Xt = XInew{1}; Yt = YInew{1};
            Anew = polyarea(Xt,Yt);
            if Anew/floe.area > 0.9
                [xm,ym] = centroid(polyshape(Xt,Yt));
                floe.Xi = xm; floe.Yi = ym;
                floe.area = polyarea(Xt,Yt);
                floe.c_alpha =[ Xt'-xm;Yt'-ym];
            end
        end
    end
    
    in = 0;
    count = 0;
    while sum(in)<0.5
        X = floe.rmax*(2*rand(N,1)-1);
        Y = floe.rmax*(2*rand(N,1)-1);
        in = inpolygon(X,Y,floe.c_alpha(1,:)',floe.c_alpha(2,:));
        count = count+1;
    end
    
    %Create a box to be used that bounds the polyshape
    boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*floe.rmax;
    worked = 1;
    while worked > 0.5
        [~, b,~,~,worked] = polybnd_voronoi([X Y],[floe.c_alpha]');
        if worked == 1
            X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
            Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);
        end
    end
    for i =1:length(b)
        a=regions(intersect(polyshape(floe.c_alpha'),polyshape(b{i})));
        
        %%Loop through all the new shapes to calculate the new properties of
        %%each
        for p=1:length(a)
            FloeNEW.poly = rmholes(a(p));
            [Xi,Yi] = centroid(FloeNEW.poly);
            FloeNEW.area = area(FloeNEW.poly);
            FloeNEW.mass = floe.mass*area(a(p))/floe.area;
            FloeNEW.h = floe.mass*area(a(p))/(rho_ice*FloeNEW.area*floe.area);
            FloeNEW.c_alpha = [(FloeNEW.poly.Vertices-[Xi Yi])' [FloeNEW.poly.Vertices(1,1)-Xi; FloeNEW.poly.Vertices(1,2)-Yi]];
            FloeNEW.c0 = FloeNEW.c_alpha;
            FloeNEW.inertia_moment = PolygonMoments(FloeNEW.c0',FloeNEW.h);
            FloeNEW.angles = polyangles(FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
            FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
            FloeNEW.strain = floe.strain;
       	    FloeNEW.Stress = zeros(2);
      	    FloeNEW.StressH = zeros(2,2,Nz);
            FloeNEW.StressCount = 1;
            FloeNEW.FxOA = 0; FloeNEW.FyOA = 0; FloeNEW.torqueOA = 0;
            
            err = 1;
            while err > 0.1
                FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1);
                FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1);
                FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
                err = (sum(FloeNEW.A)/1000*4*FloeNEW.rmax^2-FloeNEW.area)/FloeNEW.area;
            end
            
            FloeNEW.Xi = floe.Xi+Xi; FloeNEW.Yi = floe.Yi+Yi; FloeNEW.alive = 1;
            FloeNEW.alpha_i = 0; FloeNEW.Ui = floe.Ui; FloeNEW.Vi = floe.Vi;
            FloeNEW.dXi_p = floe.dXi_p; FloeNEW.dYi_p = floe.dYi_p;
            FloeNEW.dUi_p = floe.dUi_p; FloeNEW.dVi_p = floe.dVi_p;
            FloeNEW.dalpha_i_p = 0; FloeNEW.ksi_ice = floe.ksi_ice;
            FloeNEW.dksi_ice_p = floe.dksi_ice_p;
            FloeNEW.interactions = [];
            FloeNEW.potentialInteractions = [];
            FloeNEW.collision_force = 0;
            FloeNEW.collision_torque = 0;
            FloeNEW.OverlapArea = 0;          


            Floes = [Floes FloeNEW];
            clear FloeNEW
        end
        
    end
end

if ~isempty(Floes)
    Floes=rmfield(Floes,{'poly'});
end

warning('on',id)
warning('on',id3)
end
