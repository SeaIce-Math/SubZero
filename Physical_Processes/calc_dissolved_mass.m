function [Vd] = calc_dissolved_mass(Floe,Nx,Ny,c2_boundary_poly)
%%This function takes the floes that are to small to be kept track of and
%%puts it in the correct coarse grained bin

%Create coarse grids
x = min(c2_boundary_poly.Vertices(:,1)):(max(c2_boundary_poly.Vertices(:,1))-min(c2_boundary_poly.Vertices(:,1)))/Nx:max(c2_boundary_poly.Vertices(:,1));
y = min(c2_boundary_poly.Vertices(:,2)):(max(c2_boundary_poly.Vertices(:,2))-min(c2_boundary_poly.Vertices(:,2)))/Ny:max(c2_boundary_poly.Vertices(:,2));

%Find floes and create bins 
Xi=cat(1,Floe.Xi);
Yi=cat(1,Floe.Yi);
Binx = fix((Xi-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Yi-min(y))/(max(y)-min(y))*Ny+1);

% Idenfity floes that are alive
live = cat(1,Floe.alive);

%Place all floes in correct bins
Vd = zeros(Ny,Nx);
for ii = 1:Nx
    for jj = 1:Ny
        Vd(jj,ii) = sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).mass));
    end
end

end

