function [floe,FxOA,FyOA] =calc_trajectory(dt,ocean,winds,floe,HFo, doInt)

ext_force=floe.collision_force;
ext_torque=floe.collision_torque;
HFo = mean(HFo(:));

%Calcualte the current stress from floe interactions at this time step and
%average with existing
if ~isempty(floe.interactions)
    a=floe.interactions;
    r=[floe.Xi floe.Yi];
    Stress =1/(2*floe.area*floe.h)*([sum((a(:,4)-r(1)).*a(:,2)) sum((a(:,5)-r(2)).*a(:,2)); sum((a(:,4)-r(1)).*a(:,3)) sum((a(:,5)-r(2)).*a(:,3))]...
        +[sum(a(:,2).*(a(:,4)-r(1))) sum(a(:,3).*(a(:,4)-r(1))); sum(a(:,2).*(a(:,5)-r(2))) sum(a(:,3).*(a(:,5)-r(2)))]);
    [~,~,Nz] = size(floe.StressH);
    if floe.StressCount > Nz
        floe.StressCount = 1;
    end
    floe.StressH(:,:,floe.StressCount) = Stress;
    floe.StressCount = floe.StressCount + 1;
    floe.Stress = mean(floe.StressH,3);
end

if length(ext_force) == 1
    ext_force = [0 0];
end

%Bound values from being to large or small
if floe.h > 10
    floe.h = 10;
elseif floe.mass<100
    floe.mass = 1e3;
    floe.alive=0;
end
while max((abs(ext_force))) > floe.mass/(5*dt)
    ext_force = ext_force/10;
    ext_torque = ext_torque/10;
    if ~isempty(floe.interactions); a = a/10; end
end

Xo=ocean.Xo;
Yo=ocean.Yo;
Uocn=ocean.Uocn;
Vocn=ocean.Vocn;
dXo=Xo(2)-Xo(1);

Xi=floe.Xi;
Yi=floe.Yi;

% ice-ocean parameters
rho_ice=920; % kg/m3
rho0=1027;   % ocean density
Cd=3e-3;

% ice-water drag coefficient
rho_air=1.2;
Cd_atm=1e-3;

fc=ocean.fCoriolis; %coriolis parameter

%% ice floe params

floe_area=floe.area;
floe_mass=floe.mass; % total mass
h = floe.h;
floe_inertia_moment=floe.inertia_moment; % moment of inertia

%% update values based upon thermodynamic growth
dh = HFo*dt./h;
floe_mass = (h-dh)./h.*floe_mass; floe.mass = floe_mass;
floe_inertia_moment = (h-dh)./h.*floe_inertia_moment;
floe.inertia_moment = floe_inertia_moment;
floe.h = h-dh;
R_floe=sqrt(2)*floe.rmax;

%% atmospheric winds
Uwinds=winds.u; 
Vwinds=winds.v; 


%% Update trajectory
if isnan(floe.Xi), disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];
else
    
    % Only need to update interactions with ocean on a shorter time scale
    % so check to see if these needs done here
    if doInt.flag || isempty(floe.FxOA) || floe.h < 0.1
        x = floe.X;
        y = floe.Y;
        A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)]; %rotation matrix
        xr = A_rot*[x';y'];
        
        if sum(floe.A(:))==0
            err = 1;
            count = 1;
            while err > 0.1
                floe.X = floe.rmax*(2*rand(1000,1) - 1);
                floe.Y = floe.rmax*(2*rand(1000,1) - 1);
                floe.A = inpolygon(floe.X,floe.Y,floe.c_alpha(1,:),floe.c_alpha(2,:));
                err = abs((sum(floe.A)/1000*4*floe.rmax^2-floe.area))/floe.area;
                count = count+1; if count>10; err = 0; floe.alive = 0; end
            end
        end
    elseif isempty(floe.FxOA)
        floe.FxOA = 0; floe.FyOA = 0;
        
    end
    
    if  ( max(max(floe.c_alpha(1,:)))+floe.Xi>max(Xo) || min(min(floe.c_alpha(1,:)))+floe.Xi<min(Xo) || max(max(floe.c_alpha(2,:)))+floe.Yi>max(Yo) || min(min(floe.c_alpha(2,:)))+floe.Yi<min(Yo)   )       
        disp('Ice floe sacked: out of ocean grid bounds!'); floe=[];        
    elseif floe.alive == 1
        
        %Calculate forces from ocean/atmospheric stresses
        if doInt.flag  || isempty(floe.FxOA) || floe.h < 0.1
            Xg = xr(1,:)+Xi; Yg = xr(2,:)+Yi;
            
            [theta,rho] = cart2pol(xr(1,:),xr(2,:));
            
            
            Uice=floe.Ui-rho*floe.ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
            Vice=floe.Vi+rho*floe.ksi_ice.*cos(theta); % Y-dir velocity
                    
            % interpolating ocean currents onto ice floe grid.
            x_ind=logical((Xo <= Xi+R_floe+2*dXo).*(Xo >= Xi-R_floe-2*dXo));
            y_ind=logical((Yo <= Yi+R_floe+2*dXo).*(Yo >= Yi-R_floe-2*dXo));
            
            Uocn_interp=interp2(Xo(x_ind),Yo(y_ind), Uocn(y_ind,x_ind),Xg,Yg);
            Vocn_interp=interp2(Xo(x_ind),Yo(y_ind), Vocn(y_ind,x_ind),Xg,Yg);
            Uwinds_interp=interp2(Xo(x_ind),Yo(y_ind), Uwinds(y_ind,x_ind),Xg,Yg);
            Vwinds_interp=interp2(Xo(x_ind),Yo(y_ind), Vwinds(y_ind,x_ind),Xg,Yg);

            U10 = mean(Uwinds_interp(floe.A)); V10 = mean(Vwinds_interp(floe.A));
            Fx_atm=rho_air*Cd_atm*sqrt(U10^2+V10^2)*U10;
            Fy_atm=rho_air*Cd_atm*sqrt(U10^2+V10^2)*V10;
            
            Fx_pressureGrad=-(floe_mass/floe_area)*fc*Vocn_interp; % SSH tilt term
            Fy_pressureGrad=+(floe_mass/floe_area)*fc*Uocn_interp;        
        
            du=Uocn_interp-Uice; dv=Vocn_interp-Vice;            
        
            tau_ocnX=rho0*Cd*sqrt(du.^2+dv.^2).*( cos(ocean.turn_angle)*du-sin(ocean.turn_angle)*dv); % ocean stress with the turning angle
            tau_ocnY=rho0*Cd*sqrt(du.^2+dv.^2).*(sin(ocean.turn_angle)*du+cos(ocean.turn_angle)*dv);
            
            Fx=tau_ocnX+Fx_atm+Fx_pressureGrad; 
            Fy=tau_ocnY+Fy_atm+Fy_pressureGrad;
            
            
            % updating the ice floe vorticity with averaged torques over the ice floe area
            torque=(-Fx.*sin(theta)+Fy.*cos(theta)).*rho;  % torque
            
            %adding the remaining Coriolis force; it has no torque.
            Fx=Fx+(floe_mass/floe_area)*fc*floe.Vi;
            Fy=Fy-(floe_mass/floe_area)*fc*floe.Ui;
            
                                    
            floe.FxOA = mean(Fx(floe.A));
            floe.FyOA = mean(Fy(floe.A));
            floe.torqueOA = mean(torque(floe.A));
        end
            
        
        
        %Using 2nd order time-stepping here, utilizing tendencies calculated at
        %the previos time steps d = 1.5*dt*(d/dt)-0.5*dt*(d/dt)_previous
        
        % updating the ice floe coordinates with velocities
        dx =1.5*dt*floe.Ui -0.5*dt*floe.dXi_p; dy = 1.5*dt*floe.Vi -0.5*dt*floe.dYi_p;
        floe.Xi=floe.Xi+dx;  floe.dXi_p=floe.Ui;
        floe.Yi=floe.Yi+dy;  floe.dYi_p=floe.Vi;
        floe.alpha_i=floe.alpha_i+1.5*dt*floe.ksi_ice-0.5*dt*floe.dalpha_i_p; floe.dalpha_i_p=floe.ksi_ice;
        
            
        % updating the ice floe velocities with mean forces and torques
        dUi_dt=(floe.FxOA*floe_area+ext_force(1))/floe_mass;
        dVi_dt=(floe.FyOA*floe_area+ext_force(2))/floe_mass;
        frac = []; frac1 = [];frac2 = [];
        if abs(dt*dUi_dt) > 0.5*floe.h && abs(dt*dVi_dt) > 0.5*floe.h
            dUi_dt = sign(dUi_dt)*0.5*floe.h/dt;    dVi_dt = sign(dVi_dt)*0.5*floe.h/dt;
            frac1 = dUi_dt/(floe.FxOA*floe_area+ext_force(1))*floe_mass;
            frac2 = dVi_dt/(floe.FyOA*floe_area+ext_force(2))*floe_mass;
            frac = min([frac1, frac2]);
            dUi_dt=(floe.FxOA*floe_area+ext_force(1))/floe_mass;
            dVi_dt=(floe.FyOA*floe_area+ext_force(2))/floe_mass;
            dUi_dt = frac*dUi_dt; dVi_dt = frac*dVi_dt; 
        elseif abs(dt*dUi_dt) > 0.5*floe.h && abs(dt*dVi_dt) < 0.5*floe.h
            dUi_dt = sign(dUi_dt)*0.5*floe.h/dt;   
            frac = dUi_dt/(floe.FxOA*floe_area+ext_force(1))*floe_mass;
            dUi_dt=(floe.FxOA*floe_area+ext_force(1))/floe_mass;
            dVi_dt=(floe.FyOA*floe_area+ext_force(2))/floe_mass;
            dUi_dt = frac*dUi_dt; dVi_dt = frac*dVi_dt; 
        elseif abs(dt*dUi_dt) < 0.5*floe.h && abs(dt*dVi_dt) > 0.5*floe.h
            dVi_dt = sign(dVi_dt)*0.5*floe.h/dt;
            frac = dVi_dt/(floe.FyOA*floe_area+ext_force(2))*floe_mass;
            dUi_dt=(floe.FxOA*floe_area+ext_force(1))/floe_mass;
            dVi_dt=(floe.FyOA*floe_area+ext_force(2))/floe_mass;
            dUi_dt = frac*dUi_dt; dVi_dt = frac*dVi_dt; 
        end
        floe.Ui=floe.Ui+1.5*dt*dUi_dt-0.5*dt*floe.dUi_p;  
        floe.Vi=floe.Vi+1.5*dt*dVi_dt - 0.5*dt*floe.dVi_p;
        floe.dUi_p=dUi_dt;
        floe.dVi_p=dVi_dt;
        
        dksi_ice_dt=(floe.torqueOA*floe_area+ext_torque)/floe_inertia_moment;
        if ~isempty(frac)
            dksi_ice_dt = frac*dksi_ice_dt;
        end
        ksi_ice=floe.ksi_ice+1.5*dt*dksi_ice_dt - 0.5*dt*floe.dksi_ice_p;
        if abs(ksi_ice) > 1e-5
            ksi_ice = sign(ksi_ice)*1e-5;
        end
        floe.ksi_ice = ksi_ice;
        floe.dksi_ice_p=dksi_ice_dt;
        
        A_rot=[cos(floe.alpha_i) -sin(floe.alpha_i); sin(floe.alpha_i) cos(floe.alpha_i)]; %rotation matrix
        floe.c_alpha=A_rot*floe.c0; %rotate floe contour
            
        if doInt.flag
            [theta,rho] = cart2pol(floe.c_alpha(1,:),floe.c_alpha(2,:));
            
            Uice=floe.Ui-rho*floe.ksi_ice.*sin(theta); % X-dir floe velocity (variable within the ice floe)
            Vice=floe.Vi+rho*floe.ksi_ice.*cos(theta); % Y-dir velocity
            du_dx = 0.5*sum(diff([Uice Uice(1)]).*diff([floe.c_alpha(2,:) floe.c_alpha(2,1)]))/floe_area;
            du_dy = 0.5*sum(diff([Uice Uice(1)]).*diff([floe.c_alpha(1,:) floe.c_alpha(1,1)]))/floe_area;
            dv_dx = 0.5*sum(diff([Vice Vice(1)]).*diff([floe.c_alpha(2,:) floe.c_alpha(2,1)]))/floe_area;
            dv_dy = 0.5*sum(diff([Vice Vice(1)]).*diff([floe.c_alpha(1,:) floe.c_alpha(1,1)]))/floe_area;
            floe.strain = 1/2*([du_dx du_dy;dv_dx dv_dy ] + [du_dx dv_dx; du_dy dv_dy]);
        end
    end
end


FxOA = ext_force(1);
FyOA = ext_force(2);
end


