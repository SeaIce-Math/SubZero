function [FloeWhole,FloeNEW] = ridge_values_update(FloeWhole,FloeLess, V)
%%This function calculates the values for the floes after their shapes and
%%mass get updated after ridging


%Store a copy of this structure for later reference
Floe1 = FloeWhole;
Floe2 = FloeLess;

%   Calculate values for Floe gaining the new ice
rho_ice = 920;
hold = FloeWhole.h;
FloeWhole.h = FloeWhole.h + V/FloeWhole.area;
if FloeWhole.h > 30
    FloeWhole.h = 30;
end
FloeWhole.mass = FloeWhole.mass+V*rho_ice;
FloeWhole.inertia_moment = FloeWhole.h/hold.*FloeWhole.inertia_moment;

%Calculate values for floe losing ice
[poly2new] = subtract(FloeLess.poly,FloeWhole.poly);
polyout = sortregions(poly2new,'area','descend');
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
        FloeLess = Floe2;
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
            FloeLess.mass = FloeLess.area/Atot*(Floe2.mass-V*rho_ice);
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
    FloeNEW = Floe2;
    FloeNEW.alive = 0;
end

end

