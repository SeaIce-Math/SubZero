function c2_boundary =initialize_boundaries()

%Adding walls around the domain
Lx=1e5; Ly=1e5;
x=[-1 -1 1 1 -1]*Lx; 
y=[-1 1 1 -1 -1]*Ly;

%%Inertial
% R = 65e4;
% t = 0:pi/50:2*pi;
% x = R*cos(t); y = R*sin(t);

% %%Nares
% Lx=2e5; Ly=1e6;
% x=[-1 -1 1 1 -1]*Lx/2; 
% y=[-1 1 1 -1 -1]*Ly/2;

c2_boundary = [x; y];


end