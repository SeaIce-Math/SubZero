function [Floes] = Fuse_Floes(floe1,floe2)
%This function takes two input floes and fuses them together while
%conserving mass and momentum
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
save('FuseFloesArctoc.mat','floe1','floe2');
rho_ice = 920;
mass = cat(1,floe2.mass);

%Create and initialize the new shape
polynew = union([floe1.poly floe2.poly]);
polyout = sortregions(polynew,'area','descend');
[Xt, Yt] = centroid(polynew);
R = regions(polyout);
R = R(area(R)>1e4);
Atot = 0;
for ii = 1:length(R)
    Atot = area(R(ii),1)+Atot;
end
Floes = [];
mtot = floe1.mass + sum(mass);
for ii = 1:length(R)
    poly1new = R(ii);
    poly1new = rmholes(poly1new);
    poly = simplify(poly1new);
    massNew = area(R(ii))/Atot*mtot;
    height.mean = massNew/(area(poly)*rho_ice);
    height.delta = 0;
    floenew = initialize_floe_values(poly, height);
    Floes = [Floes floenew];
end

%% Calculate new properties to conserve the momemtum of the existing floes
Ui = (floe1.Ui*floe1.mass + sum(cat(1,floe2.Ui).*mass))/(mtot);
Vi = (floe1.Vi*floe1.mass + sum(cat(1,floe2.Vi).*mass))/(mtot);
for ii = 1:length(R)
    h2(ii) = (Xt-Floes(ii).Xi).^2+(Yt-Floes(ii).Yi).^2;
end
inertia_ice = sum(cat(1,Floes.inertia_moment))+sum(cat(1,Floes.mass).*h2');
dUi_p = (floe1.dUi_p*floe1.mass + sum(cat(1,floe2.dUi_p).*mass))/(mtot);
dVi_p = (floe1.dVi_p*floe1.mass + sum(cat(1,floe2.dVi_p).*mass))/(mtot);
ksi_ice = (floe1.ksi_ice*floe1.inertia_moment + sum(cat(1,floe2.ksi_ice).*cat(1,floe2.inertia_moment)))/(inertia_ice);%use inertia moment instead of mass
dXi_p = (floe1.dXi_p*floe1.mass + sum(cat(1,floe2.dXi_p).*mass))/(mtot);
dYi_p = (floe1.dYi_p*floe1.mass + sum(cat(1,floe2.dYi_p).*mass))/(mtot);
dksi_ice_p = (floe1.dksi_ice_p*floe1.inertia_moment + sum(cat(1,floe2.dksi_ice_p).*cat(1,floe2.inertia_moment)))/(inertia_ice);


for ii = 1:length(R)
    Floes(ii).Ui = Ui;
    Floes(ii).Vi = Vi;
    Floes(ii).dUi_p = dUi_p;
    Floes(ii).dVi_p = dVi_p;
    Floes(ii).ksi_ice = ksi_ice;
    Floes(ii).dXi_p = dXi_p;
    Floes(ii).dYi_p = dYi_p;
    Floes(ii).dksi_ice_p = dksi_ice_p;
    Floes(ii).FxOA = [];
    Floes(ii).FyOA = [];
    Floes(ii).Stress = floe1.Stress*floe1.mass;
    Floes(ii).StressH = floe1.StressH*floe1.mass;
    for jj = 1:length(floe2)
        Floes(ii).Stress = Floes(ii).Stress+floe2(jj).Stress*floe2(jj).mass;
        Floes(ii).StressH = Floes(ii).StressH+floe2(jj).StressH*floe2(jj).mass;
    end
    Floes(ii).Stress = Floes(ii).Stress/mtot;
    Floes(ii).StressH = Floes(ii).StressH/mtot;
    Floes(ii).StressCount = 1;

end

end

