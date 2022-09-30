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
Tice = -20; Tocean = 3*ones(size(Xocn));
heat_flux = 7.4*10^(-4)*(Tice-Tocean)/(72); %cm^2/s
heat_flux = heat_flux/100^2; %m^2/s
h0 = real(sqrt(-2*dt*heat_flux*nDTOut));
h0 = mean(h0(:)); %Thickness of newly created sea ic

end
