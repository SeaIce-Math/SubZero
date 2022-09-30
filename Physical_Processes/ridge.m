function [Floe1,Floe2]= ridge(Floe1,Floe2,c2_boundary_poly,PERIODIC,min_floe_size)
%% This function takes in two floes and based upon the thickness of the two floes will perform a ridging operation
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id);
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3);
%Create polyshapes for the boundary in a nonperiodic case as well as polyshapes for the union of all subfloes 
floe1 = Floe1;
floe2 = Floe2;

poly1 = simplify(polyshape(Floe1.c_alpha'+[Floe1.Xi Floe1.Yi]));
poly2 = simplify(polyshape(Floe2.c_alpha'+[Floe2.Xi Floe2.Yi]));
Floe1.poly = poly1;
Floe2.poly = poly2;


%Find area of overlap
aPoly = area(intersect(poly1,poly2));

%Find critical thickness
rho_ice=920;
rho_l = 997;
E = 1e9;
sigma_m = 400000;
nu = 0.29;
g = 9.81;
hc = 0.2;%14.2*(1-nu^2)/(rho_l*g)*sigma_m^2/E;
disolved = 0;
boundary = 0;

%check to make sure one floe is not inside the other and add mass to
%dissolved if it is
if area(poly2)/area(c2_boundary_poly)>0.99
    boundary = 1;
elseif aPoly/area(Floe1.poly)>0.75 || Floe1.area<min_floe_size
    disolved = 1;
    Floe1.alive = 0;
elseif aPoly/area(Floe2.poly)>0.75 || Floe2.area < min_floe_size
    disolved = 1;
    Floe2.alive = 0;
elseif Floe1.alive+Floe2.alive < 2
    disolved = 0;
end
%% If there is enough overlap then allow ridging to happen


if disolved == 0 && aPoly > 500 && ~boundary
    
    %Determine overlap in the different subfloes
    V1 = aPoly*Floe1.h;
    V2 = aPoly*Floe2.h;
    
    %Use the thicknesses to determine how mass will be transfered
    if Floe1.h>= hc && Floe2.h >= hc
        p=1/(1+Floe1.h/Floe2.h);
        if rand(1)>= p
            [Floe1, Floe2] = ridge_values_update(Floe1,Floe2, V2);
        else
            [Floe2, Floe1] = ridge_values_update(Floe2,Floe1, V1);
        end
    elseif Floe1.h>= hc && Floe2.h< hc
        [Floe1, Floe2] = ridge_values_update(Floe1,Floe2, V2);
    elseif Floe1.h < hc && Floe2.h >= hc
        [Floe2, Floe1] = ridge_values_update(Floe2,Floe1, V1);
    end
end
    

%Perform ridging with boundary
if ~PERIODIC && boundary
    Lx= max(c2_boundary_poly.Vertices(:,1));
    Ly= max(c2_boundary_poly.Vertices(:,2));%c2 must be symmetric around x=0 for channel boundary conditions.
    x=[-1 -1 1 1 -1]*Lx*2;
    y=[-1 1 1 -1 -1]*Ly*2;
    polybound = polyshape(x,y);
    c2_poly = subtract(polybound,c2_boundary_poly);
    Abound = area(intersect(Floe1.poly,c2_poly));
    poly1new = subtract(Floe1.poly,c2_poly);
    V1 = Abound*Floe1.h; 

    if area(poly1new) < 1e4
        FloeNEW = Floe1;
        FloeNEW.alive = 0;
    elseif Abound>0 && area(poly1new) >= 1e4
        %Calculate new floe properties after the mass transfer and shape is
        %updated
        polyout = sortregions(poly1new,'area','descend');
        R = regions(polyout);
        R = R(area(R)>1e4);
        R = rmholes(R);
        Atot = 0;
        for ii = 1:length(R)
            Atot = area(R(ii),1)+Atot;
        end
        if ~isempty(area(R))
            FloeNEW = [];
            for kk = 1:length(R)
                FloeLess = Floe1;
                poly2new = rmholes(R(kk));
                FloeLess.poly = poly2new;
                if area(poly2new) < 10
                    FloeLess.alive = 0;
                else
                    FloeLess.area = area(poly2new);
                    [Xi,Yi] = centroid(poly2new);
                    FloeLess.Xi = Xi;
                    FloeLess.Yi = Yi;
                    FloeLess.dalpha_i_p = 0;
                    FloeLess.alpha_i = 0;
                    FloeLess.mass = FloeLess.area/Atot*Floe1.mass;
                    FloeLess.c_alpha = [(poly2new.Vertices-[FloeLess.Xi FloeLess.Yi])' [poly2new.Vertices(1,1)-FloeLess.Xi; poly2new.Vertices(1,2)-FloeLess.Yi]];
                    FloeLess.angles = polyangles(poly2new.Vertices(:,1),poly2new.Vertices(:,2));
                    FloeLess.c0 = FloeLess.c_alpha;
                    FloeLess.h = FloeLess.mass/(rho_ice*FloeLess.area);
                    if poly2new.NumRegions > 1
                        if area(poly2new) > 1000
                            FloeLess.c_alpha = poly2new;
                            FloeLess.h = FloeLess.mass/(rho_ice*Atot);
                        end
                    end
                end

                FloeLess.inertia_moment = PolygonMoments(FloeLess.c_alpha',FloeLess.h);
                FloeLess.rmax = max(sqrt(FloeLess.c_alpha(1,:).^2+FloeLess.c_alpha(2,:).^2));
                FloeLess.X = FloeLess.rmax*(2*rand(1000,1) - 1);
                FloeLess.Y = FloeLess.rmax*(2*rand(1000,1) - 1);
                FloeLess.A = inpolygon(FloeLess.X,FloeLess.Y,FloeLess.c_alpha(1,:),FloeLess.c_alpha(2,:));
                FloeNEW = [FloeNEW FloeLess];
            end
        else
            FloeNEW = Floe1;
            FloeNEW.alive = 0;
        end
    else
        FloeNEW = Floe1;
    end
    Floe1 = FloeNEW;
end


if isfield(Floe1,'poly')
    Floe1=rmfield(Floe1,{'poly'});
end
if isfield(Floe2,'poly')
    Floe2=rmfield(Floe2,{'poly'});
end


end
