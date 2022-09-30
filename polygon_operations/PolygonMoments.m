function Izz = PolygonMoments (xy,h)
if nargin==0
error('*** PolygonMoments: missing all arguments')
elseif nargin==1
    mn=[];
    PlotFlag=0;
elseif nargin<=2
    PlotFlag=0;
elseif nargin>3
    error('*** PolygonMoments: Number of arguments > 3')
end
if PlotFlag<0, NoIds=1;PlotFlag=abs(PlotFlag);else NoIds=0;end
%
rho_ice = 920;
N=size(xy,1);
x=xy(:,1);
y=xy(:,2);
% make sure the last point is the first too
x=[x;x(1)];
y=[y;y(1)];
% Compute common values and Area Centroid (AC) always
wi=x(1:N).*y(2:N+1)-x(2:N+1).*y(1:N);
% area of the polygon
Polygon.Area=0.5*sum(wi);
% area moments
Polygon.MAx=1/6 *sum(wi.*(y(1:N) + y(2:N+1)));
Polygon.MAy=1/6 *sum(wi.*(x(1:N) + x(2:N+1)));
% Area moments of second degree (Inertia)
Polygon.Ixx=1/12*sum(wi.*(  (y(1:N)+y(2:N+1)).^2 - y(1:N).*y(2:N+1)));
Polygon.Iyy=1/12*sum(wi.*(  (x(1:N)+x(2:N+1)).^2 - x(1:N).*x(2:N+1)));
Polygon.Ixy=1/24*sum(wi.*( (x(1:N)+x(2:N+1)).*(y(1:N)+y(2:N+1))  + x(1:N).*y(1:N) + x(2:N+1).*y(2:N+1) ));
Izz=abs(Polygon.Ixx+Polygon.Iyy)*h*rho_ice;
% Izz=abs(Polygon.Ixx+Polygon.Iyy)*rho_ice;
% coordinates of the area centroid:
% Polygon.ACx=Polygon.MAy/Polygon.Area;
% Polygon.ACy=Polygon.MAx/Polygon.Area;
% % area moments in the area centroid:
% Polygon.IxxAC=Polygon.Ixx-Polygon.ACy^2*Polygon.Area;
% Polygon.IyyAC=Polygon.Iyy-Polygon.ACx^2*Polygon.Area;
% Polygon.IxyAC=Polygon.Ixy-Polygon.ACx*Polygon.ACy*Polygon.Area;
% IzzAC=abs(Polygon.IxxAC+Polygon.IyyAC)*h*rho_ice;

end