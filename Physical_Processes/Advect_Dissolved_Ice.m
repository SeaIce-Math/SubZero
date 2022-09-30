function Vd_new = Advect_Dissolved_Ice(Vd,coarseMean,im_num,dissolved_new,c2_boundary,dt)
%This function takes in the the variable containing the mass of dissolved
%sea ice on an eularian grid as well as the coarse mean eularian velocity.
%The dissolved mass of sea ice is then advected around the floe field per
%this information.

%% find current terms
U = squeeze(coarseMean(2,:,:,im_num));
V = squeeze(coarseMean(3,:,:,im_num));
vdcurrent = Vd(:,:,1);

%Add ghost cells for current time step
[Ny,Nx] = size(U);
Nx = Nx+4;
Ny = Ny+4;
u(3:Ny-2,3:Nx-2) = U;
u(1:2,:) = [u(3,:); u(3,:)];
u(Ny-1:Ny,:) = [u(Ny-2,:);u(Ny-2,:)];
u(:,1:2) = [u(:,3) u(:,3)];
u(:,Nx-1:Nx) = [u(:,Nx-2) u(:,Nx-2)];
v(3:Ny-2,3:Nx-2) = V;
v(1:2,:) = [v(3,:); v(3,:)];
v(Ny-1:Ny,:) = [v(Ny-2,:);v(Ny-2,:)];
v(:,1:2) = [v(:,3) v(:,3)];
v(:,Nx-1:Nx) = [v(:,Nx-2) v(:,Nx-2)];
Vdcurrent(3:Ny-2,3:Nx-2) = vdcurrent;
Vdcurrent(1:2,:) = [Vdcurrent(3,:); Vdcurrent(3,:)];
Vdcurrent(Ny-1:Ny,:) = [Vdcurrent(Ny-2,:);Vdcurrent(Ny-2,:)];
Vdcurrent(:,1:2) = [Vdcurrent(:,3) Vdcurrent(:,3)];
Vdcurrent(:,Nx-1:Nx) = [Vdcurrent(:,Nx-2) Vdcurrent(:,Nx-2)];


%define strength of diffusion
diffusion = 1e4;

DissolvedNew = zeros(Ny,Nx);
DissolvedNew(3:Ny-2,3:Nx-2) = dissolved_new;

%% Define grid
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/(Nx-2):max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/(Ny-2):max(c2_boundary(2,:));
delx = abs(x(2)-x(1));
dely = -abs(y(2)-y(1));

%% Find differentiation matrix for x
d2x = zeros(Nx,Nx);
b = 1/(delx^2);
d2x(2:Nx+1:Nx*Nx) = b;
d2x(2+Nx:Nx+1:Nx*Nx) = -2*b;
d2x(2+2*Nx:Nx+1:Nx*Nx) = b;
d2x(1,1) = b/(3/2);
d2x(1,2) = -2*b/(3/2);
d2x(1,3) = b/(3/2);
d2x(Nx,Nx) = b/(3/2);
d2x(Nx,Nx-1) = -2*b/(3/2);
d2x(Nx,Nx-2) = b/(3/2);

%Find differentiation matrix for y
d2y = zeros(Ny,Ny);
b = 1/(dely^2);
d2y(2:Ny+1:Ny*Ny) = b;
d2y(2+Ny:Ny+1:Ny*Ny) = -2*b;
d2y(2+2*Ny:Ny+1:Ny*Ny) = b;
d2y(1,1) = b/(3/2);
d2y(1,2) = -2*b/(3/2);
d2y(1,3) = b/(3/2);
d2y(Ny,Ny) = b/(3/2);
d2y(Ny,Ny-1) = -2*b/(3/2);
d2y(Ny,Ny-2) = b/(3/2);

Ix = eye(Nx);
Iy = eye(Ny);

%% Shift points for calculating Advective Terms
Vxshift = zeros(Ny-4,Nx-3);
Vxshift(:,2:end-1) = 0.5*(Vdcurrent(3:Ny-2,3:Nx-3)+Vdcurrent(3:Ny-2,4:Nx-2));
Ushift = zeros(Ny-4,Nx-3);
Ushift(:,2:end-1) = 0.5*(u(3:Ny-2,3:Nx-3)+u(3:Ny-2,4:Nx-2));
Vyshift = zeros(Ny-3,Nx-4);
Vyshift(2:end-1,:) = 0.5*(Vdcurrent(3:Ny-3,3:Nx-2)+Vdcurrent(4:Ny-2,3:Nx-2));
Vshift = zeros(Ny-3,Nx-4);
Vshift(2:end-1,:) = 0.5*(v(3:Ny-3,3:Nx-2)+v(4:Ny-2,3:Nx-2));
%% Time step new equation for dissolved ice. Adams Basheforth for NL terms and Crank-Nicolsen for linear terms
%Calculate RHS first without advection terms
RHS = Vdcurrent(:)+DissolvedNew(:)+dt*diffusion*(kron(Ix,d2y)+kron(d2x,Iy))*Vdcurrent(:);
Vd_new = reshape(RHS,Ny,Nx);

%calculate advection terms
Advec = diff(Ushift.*Vxshift,1,2)/delx+ diff(Vshift.*Vyshift,1,1)/dely;
Vd_new(3:Ny-2,3:Nx-2) = Vd_new(3:Ny-2,3:Nx-2)-dt*Advec;

%end
Vd_new = Vd_new(3:Ny-2,3:Nx-2);
Vd_new(Vd_new<0) = 0;
end