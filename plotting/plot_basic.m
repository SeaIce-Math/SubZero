function [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nb,PERIODIC)
%This function creates plots of the floe state showing the stress and and
%thickness of the floes
Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 
live = cat(1,Floe.alive);
Floe(live == 0) = [];

N0=length(Floe);
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    translation = [];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            parent=[parent  i];
            translation = [translation; -2*Lx*sign(x(i)) 0];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            parent=[parent  i];
            translation = [translation; 0 -2*Ly*sign(y(i))];

        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

%Find length of new Floe variable including the ghost floes
N=length(Floe);

%% Set Up The Plots

ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end
figure(fig)
clf(fig);

dn=1; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;

axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

%title(['Time = ' num2str(Time/3600) ' h'],'fontsize',36);
for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end
xi = cat(1,Floe.Xi); yi = cat(1,Floe.Yi); Ui = cat(1,Floe.Ui); Vi = cat(1,Floe.Vi);
%quiver(xi,yi,Ui, Vi,'linewidth',2)

grid
ax = gca;
ax.GridAlpha = 0.5;

%% Plot the Floes and Ghost Floes
plot([Floe(1:N0).poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
if PERIODIC && N > N0
    plot([Floe(N0+1:N).poly],'FaceColor','k','FaceAlpha',0.5,'EdgeColor',[1 1 1]*0.2);
end
if Nb > 0
    plot([Floe(1:Nb).poly],'FaceColor','k','FaceAlpha',0.3,'EdgeColor',[1 1 1]*0.2);
end
% if ~isempty(Floe(1).interactions)
%    cx = Floe(1).interactions(4); cy = Floe(1).interactions(5);
%    load F_dir
%    quiver(cx,cy,force_dir(1), force_dir(2),'autoscalefactor',500,'linewidth',2)
% end
% xx = 1; xx(1) =[1 2];
set(0,'CurrentFigure',fig);
xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
plot(xb,yb, 'k-','linewidth',2);

colormap('gray'); caxis([0 1]);
axis([-Lx-Lx/10 Lx+Lx/10 -Ly-Ly/10 Ly+Ly/10])
axis([-100000 100000 -100000 100000])
%xlabel('m');ylabel('m');
set(gca,'Ydir','normal');
% set(gca,'xtick',[])
% set(gca,'ytick',[])

% drawnow
end
