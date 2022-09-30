function FloeNEW = initialize_floe_values(poly1, height)
%%This function populates all the fields of the floes upon their creation

rho_ice=920;

polyout = sortregions(poly1,'area','descend');
R = regions(polyout);
poly1 = R(1);
polya = rmholes(poly1);
h=height.mean+(-1)^randi([0 1])*rand*height.delta;

FloeNEW.poly = poly1;
[Xi,Yi] = centroid(FloeNEW.poly); 
FloeNEW.area = area(FloeNEW.poly); %Floe area (m^2)
FloeNEW.h = h;  %Floe Thickness (m)
FloeNEW.mass = FloeNEW.area*h*rho_ice; %kg
FloeNEW.c_alpha = [(polya.Vertices-[Xi Yi])' [polya.Vertices(1,1)-Xi; polya.Vertices(1,2)-Yi]]; %Floe boundary rotated about a reference position
FloeNEW.c0 = FloeNEW.c_alpha; %Unroated floe boundary
FloeNEW.inertia_moment = PolygonMoments(FloeNEW.c0',h); %Moment of inertia
FloeNEW.angles = polyangles(polya.Vertices(:,1),polya.Vertices(:,2)); %Interior angles of floe vertices
FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1))); %Distance of vertix farthest from the centroid (m)
FloeNEW.Stress = [0 0; 0 0]; %Homogenized stress over floe (Pa)
FloeNEW.strain = [0 0; 0 0]; %Homogenized strain 
FloeNEW.StressH = zeros(2,2,1000); %Stress history
FloeNEW.StressCount = 1; %place counter for most recent stress 
FloeNEW.FxOA = 0; FloeNEW.FyOA = 0; FloeNEW.torqueOA = 0; %Forces/torques from ocean and atmosphere

err = 1;
count = 1;
while err > 0.1
    FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1); %X coordinate for monte-carlo integration 
    FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1); %Y coordinate for monte-carlo integration
    FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:)); %Identify which are within floe boundary
    err = abs((sum(FloeNEW.A)/1000*4*FloeNEW.rmax^2-area(polya)))/area(polya);
    count = count+1; if count>10; err = 0; FloeNEW.alive = 0; end
end

FloeNEW.Xi = Xi; FloeNEW.Yi = Yi; FloeNEW.alive = 1; %centroid of the floe
FloeNEW.alpha_i = 0;  %Angle rotated from reference c0
FloeNEW.Ui = 0; FloeNEW.Vi = 0;FloeNEW.ksi_ice = 0; %Linear and angular velocities

% Previous values
FloeNEW.dXi_p = 0; FloeNEW.dYi_p = 0;
FloeNEW.dUi_p = 0; FloeNEW.dVi_p = 0;
FloeNEW.dalpha_i_p = 0; 
FloeNEW.dksi_ice_p = 0;

FloeNEW.interactions = []; %Information of interactions with other floes
FloeNEW.potentialInteractions = []; %Informatio of floes that could possibly be interacted with
FloeNEW.collision_force = 0; %Total force from collisions with floes
FloeNEW.collision_torque = 0; %Total torque from collisions with floes
FloeNEW.OverlapArea = 0; %Total overlap area with other floes


end
