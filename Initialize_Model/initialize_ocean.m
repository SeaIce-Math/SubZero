function [ocean, heat_flux, h0]=initialize_ocean(dt,nDTOut)

% defining ocean currents
ocean.fCoriolis=1.4e-4; % Coriolis parameter.

ocean.U = 0.5;

ocean.turn_angle=15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% ocean grid;
dXo=10000; % in meters

Lx = 4e5; kx = pi/Lx; Ly = 4e5; ky = pi/Ly;
Xo=-Lx:dXo:Lx; Yo=-Ly:dXo:Ly; 
[Xocn, Yocn]=meshgrid(Xo,Yo);

%defining ocean streamfunction with some eddies
transport=0.5e4; % horizontal transport, in m^2/s (controls ocean currents) 
psi_ocean=transport/1*(sin(4*kx*Xocn).*sin(4*ky*Yocn));

%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=zeros(size(Xocn));
Uocn(2:end,:)=-(psi_ocean(2:end,:)-psi_ocean(1:end-1,:))/dXo; 
Vocn(:,2:end)=(psi_ocean(:,2:end)-psi_ocean(:,1:end-1))/dXo;

ocean.Xo=Xo;
ocean.Yo=Yo;
ocean.Xocn = Xocn;
ocean.Yocn = Yocn;
ocean.kx = kx;
ocean.ky = ky;
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;

%Calculate heat flux and how much ice would grow between creation of new
%floes
k = 2.14; %Watts/(meters Kelvin)
dt = 10; %seconds
Ta = -20; %Kelvin
To = 0; %Kelvin
rho_ice = 920; %kg/m^3
L = 2.93e5; % Joules/kg

h0 = real(sqrt(2*k*dt*nDTOut*(To-Ta)/(rho_ice*L)));
heat_flux = k*(Ta-To)/(rho_ice*L); 
h0 = mean(h0(:)); %Thickness of newly created sea ice;

end
