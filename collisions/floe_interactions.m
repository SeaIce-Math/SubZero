function [ force_1, pcontact, overlap,gam] = floe_interactions(floe1, floe2, c2_boundary,PERIODIC,Modulus,dt)
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id2 = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id2)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

% Calculate equivalent spring constants and shear modulus
h1 = floe1.h; h2 = floe2.h;
r1 = sqrt(floe1.area); r2 = sqrt(floe2.area);
Force_factor=Modulus*(h1*h2)/(h1*r2+h2*r1); overlap = 0;
if isfield(floe2,'alive')
    Force_factor=Modulus*h1/r1;
elseif r1>1e5 || r2>1e5
    r1 = min([r1 r2]);
    h1 = min([h1 h2]);
    Force_factor=Modulus*h1/r1;
end
nu = 0.3;
G = Modulus/(2*(1+nu)); mu = 0.2;%0.75;
gam = 0;

%Check if floes are actually overlapping
c1=[floe1.c_alpha(1,:)+floe1.Xi; floe1.c_alpha(2,:)+floe1.Yi];
if isfield(floe2,'c')
    c2=floe2.c;
    boundary = 0;
    [Xi,Yi] = polyclip(c1',c2','int');
else
    polyb = holes(floe2.poly);
    c2 = [polyb.Vertices]';
    boundary = 1;
    [Xi,Yi] = polyclip(c1',c2','dif');
    if ~isempty(Xi)
        X2 = Xi{1}; Y2 = Yi{1};
        if polyarea(X2,Y2)/floe1.area > 0.75
            overlap = Inf;
        end
    end
end

if isempty(Xi)
    Ar = 0;
else
    for k = 1:length(Xi)
        X = Xi{k}; Y = Yi{k};
        poly = polyshape(X,Y);
        Ar(k) = area(poly);
    end
end

%If floes are overlapping to much set to merge into one floe
if  (max(c1(1,:))<max(c2_boundary(1,:)) && min(c1(1,:))>min(c2_boundary(1,:)) && max(c1(2,:))<max(c2_boundary(2,:)) && min(c1(2,:))>min(c2_boundary(2,:))|| floe2.area<0.95*area(polyshape(c2_boundary')) || PERIODIC) 
    if sum(Ar)/floe1.area > 0.55
        overlap = Inf;
    elseif sum(Ar)/floe2.area > 0.55
        overlap = -Inf;
    end
end

if norm(c1(:,1)-c1(:,end))> 1
    c1(:,length(c1)+1) = c1(:,1);
end
if norm(c2(:,1)-c2(:,end))> 1
    c2(:,length(c2)+1) = c2(:,1);
end

% Calculate interaction forces
P=InterX(c1,c2);
if isempty(P) || size(P,2)<2 || isinf(overlap) || isempty(Xi)
    force_1=[0 0];
    pcenter=[0 0];
    pcontact=[0 0];    
else
    
    %Find number of overlapping regions
    N1 = length(c1)-1; N2 = length(c2)-1;
    Amin =  min([N1,N2])*100/1.75;
    if abs(length(Ar)-length(Xi)) > 0
        save('fail.mat','floe1','floe2','c2_boundary','PERIODIC','Modulus','dt');
    end    
    Xi(Ar<Amin) = []; Yi(Ar<Amin) = []; Ar(Ar<Amin) = [];
    N_contact=length(Xi);
    
    force_1=zeros(N_contact,2);
    
    pcenter=zeros(N_contact,2);
    pcontact=zeros(N_contact,2);
    
        
    for k=1:N_contact
        
        %Identify contact points
        X = Xi{k}; Y = Yi{k};
        poly = polyshape(X,Y);
        [cx, cy] = centroid(poly);
        [verts,dist] = dsearchn([X Y],P');
        p = [X(verts(dist<1),:) Y(verts(dist<1),:)];
        [m,~] = size(p);
        
        %Caclulate normal forces
        if Ar(k) == 0
            force_dir = [0; 0];
            pcontact(k,:) = [0, 0];
            dl = 0;
        elseif m == 2
            pcontact(k,:)= [cx, cy];
            xv = [p(1,1); p(2,1)]; yv = [p(1,2); p(2,2)];
            xgh = xv(2)-xv(1); ygh = yv(2:end)-yv(1:end-1);
            b = sqrt(xgh.^2+ygh.^2); force_dir = [-ygh./b; xgh./b];
            dl = mean(b);
        elseif m ==0
            force_dir = [0; 0];
            pcontact(k,:) = [cx, cy];
            dl = 0;
        else
            xv = [X; X(1)]; yv = [Y; Y(1)];
            xgh = xv(2:end)-xv(1:end-1); xm = (xv(2:end)+xv(1:end-1))/2;
            ygh = yv(2:end)-yv(1:end-1); ym = (yv(2:end)+yv(1:end-1))/2;
            b = sqrt(xgh.^2+ygh.^2); n = [-ygh./b xgh./b];
            xt = xm+n(:,1)/100; yt = ym+n(:,2)/100;
            in = inpolygon(xt,yt,xv,yv);
            n(~in,:) = -n(~in,:);
            Fn = -Force_factor*(b*ones(1,2)).*n;
            [d_min1] = p_poly_dist(xm, ym,c1(1,:)', c1(2,:)');
            on = logical(abs(d_min1)<1e-8);
            if sum(on)<length(d_min1)&& sum(on)>0
                f_dir = sum(Fn(on,:),1);
                force_dir=f_dir'/sqrt(f_dir*f_dir');
                dl = mean(b(on));
            else
                force_dir = [0; 0];
                dl = 0;
            end
            pcontact(k,:) = [cx, cy];
        end

        %Find direction of force
        [mf,nf] = size(force_dir);
        if dl < 0.1
            force_dir = [0; 0];
            X1new = c1(1,:);
            Y1new = c1(2,:);
        elseif mf == 2 && nf == 1
            X1new = c1(1,:)+force_dir(1);
            Y1new = c1(2,:)+force_dir(2);
        else
            f_dir = sum(Fn(on,:),1);
        end
        if boundary
            [XInew,YInew] = polyclip([X1new' Y1new'],c2','dif');
        else
            [XInew,YInew] = polyclip([X1new' Y1new'],c2','int');
        end
        for ii = 1:length(XInew)
            Xt = XInew{ii}; Yt = YInew{ii};
            [Xn,~] = polyclip([Xt Yt],[X Y],'int');
            if ~isempty(Xn)
                Anew = polyarea(Xt,Yt);
                if Anew/Ar(k)-1 > 0
                    force_dir = -force_dir;
                end
            end
        end
        
        force=force_dir*Ar(k)*Force_factor; %proportional to the overlap area
        
        % Calculate tangential forces
        v1 = ([floe1.Ui floe1.Vi]+ floe1.ksi_ice*(pcontact(k,:)-[floe1.Xi floe1.Yi]));
        v2 = ([floe2.Ui floe2.Vi]+ floe2.ksi_ice*(pcontact(k,:)-[floe2.Xi floe2.Yi]));
        v_t = (v1-v2);
        if max(abs(v_t)) == 0
            dir_t =[0 0];
        else
            dir_t = v_t/vecnorm(v_t);
        end
        force_t = -dot(dir_t,v_t)*dl*G*vecnorm(v_t)*dir_t*dt;

        if vecnorm(force_t)>mu*vecnorm(force)
            force_t = -mu*vecnorm(force)*dir_t;
%            force_t = -dot(dir_t,v_t)*mu*vecnorm(force)*dir_t;
        end
        gam(k) = vecnorm(force)/vecnorm(force_t);
        
        force_1(k,:)= force'+force_t;
        overlap(k) = Ar(k);
        
        
    end

end


warning('on',id)
warning('on',id2)
warning('on',id3)

end

