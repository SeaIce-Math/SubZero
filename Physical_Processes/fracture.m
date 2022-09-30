function [Floe,Princ] = fracture(Floe,Nb,min_floe_size,concentration)
%basic isotropic fracture mechanism is implemented based on the stress experienced by 
%floes and fractures a floe into a defined number of smaller pieces when the 
%principal stress values satisfy the specified fracture criteria

A = cat(1,Floe.area);

%   elliptical yield curve that was used in Hibler viscous-plastic rheology
Pstar = 2.25e5; C = 20;
h = mean(cat(1,Floe.h));
P = Pstar*h*exp(-C*(1-concentration));
t = linspace(0,2*pi) ;
a = P*sqrt(2)/2 ; b = a/2 ;
x = a*cos(t) ;
y = b*sin(t) ;
Mohr = polyshape(x,y);
Mohr = rotate(Mohr,45);
Mohr = translate(Mohr,[-P/2, -P/2]);

%Use Mohr's Cone
q = 5.2; SigC = 250e3;
Sig1 = (1/q+1)*SigC/(1/q-q);
Sig2 = q*Sig1+SigC;
Sig11 = -3.375e4;
Sig22 = q*Sig11+SigC;
MohrX = [Sig1; Sig11; Sig22];
MohrY = [Sig2; Sig22; Sig11];
Mohr = polyshape(-MohrX,-MohrY);

%Calculate Principal Stresses
for ii = 1:length(Floe)
    Stress = eig(Floe(ii).Stress);
    Princ(ii,1) = max(Stress);
    Princ(ii,2) = min(Stress);
    Princ(ii,3) = Floe(ii).area;
end
Princ1 = Princ(:,1); Princ2 = Princ(:,2);

%Determine if stresses are inside or outside allowable regions
[in,~] = inpolygon(Princ1,Princ2,Mohr.Vertices(:,1), Mohr.Vertices(:,2));
keep = zeros(length(Floe),1);
keep(in) = 1;
keep(A<min_floe_size)=1;
keep(1:Nb) = ones(Nb,1);
keep = logical(keep);

%Fracture those floes
for ii = 1:length(Floe)
    FracFloes(ii).floenew = [];
end
parfor ii = 1:length(keep)
    if ~keep(ii)
        FracFloes(ii).floenew=fracture_floe(Floe(ii),3,Floe);
    end
end
fracturedFloes =[];
for ii = 1:length(FracFloes)
    fracturedFloes = [fracturedFloes FracFloes(ii).floenew];
end
if isfield(fracturedFloes,'potentialInteractions')
    fracturedFloes=rmfield(fracturedFloes,{'potentialInteractions'});
end
if ~isempty(fracturedFloes)
    Floe=[Floe(keep) fracturedFloes];
end

end

