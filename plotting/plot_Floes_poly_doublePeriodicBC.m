function fig=plot_Floes_poly_doublePeriodicBC(fig, Time,Floe,ocean,eularian_data,c2_boundary_poly, PERIODIC)

showContactPoints=0;
showCenterOfMass=0;
stressColor=1;
c2_boundary = c2_boundary_poly.Vertices';

for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

N0=length(Floe);

Lx= max(c2_boundary_poly.Vertices(:,1)); %c2 must be symmetric around x=0 for channel boundary conditions.
Ly= max(c2_boundary_poly.Vertices(:,2)); 

ghostFloeX=[];
ghostFloeY=[];
parent=[];

x=cat(1,Floe.Xi);
alive=cat(1,Floe.alive);



if PERIODIC
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,1)))>Lx/2)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).poly=translate(Floe(i).poly,[-2*Lx*sign(x(i)) 0]);
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            parent=[parent  i];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,2)))>Ly/2)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).poly=translate(Floe(i).poly,[0 -2*Ly*sign(y(i))]);
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            parent=[parent  i];
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end


%ratio=max(c2_boundary(:,1:end-1),[],2)-mean(c2_boundary(:,1:end-1),2); ratio=ratio(2)/ratio(1);

ratio=max(ocean.Yo)/max(ocean.Xo);
if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 1000 1000*ratio],'visible','on');  
    set(fig,'PaperSize',12*[1 ratio],'PaperPosition',12*[0 0 1 ratio]);
end

clf(fig);

dn=1; % plot every dn'th velocity vector
hold on;
set(gca,'Color',[0.5843    0.8157    0.9882])

h = cat(1,Floe.h);
h(h>1) = 1;
%h = h/max(h);
%axis([-1 1 -1 1]*7e4);
%axis([nanmin(c2_boundary(1,:)) nanmax(c2_boundary(1,:)) nanmin(c2_boundary(2,:)) nanmax(c2_boundary(2,:))]);
axis([ocean.Xo(1) ocean.Xo(end) ocean.Yo(1) ocean.Yo(end)]);

colormap('gray'); caxis([0 1]);

title(['Time = ' num2str(Time/3600) ' hours'],'fontsize',24);

n=50;
for j=1:length(Floe)
    if Floe(j).alive
        poly=intersect(Floe(j).poly,c2_boundary_poly);
        if PERIODIC, poly_ghost=subtract(Floe(j).poly,c2_boundary_poly); end

        if stressColor==1
            cFact=min(1,Floe(j).OverlapArea/Floe(j).area); 
%             xx = 1; xx(1) =[1 2];
            plot(poly,'FaceColor',[1 1 1]*h(j),'FaceAlpha',1);%(1+h(j)*(n-1))/n);
            if PERIODIC, plot(poly_ghost,'FaceColor','k','FaceAlpha',0.4,'EdgeColor',[1 1 1]*0.4); end

        else
            plot(poly,'FaceColor','k'); 
            if PERIODIC,  plot(poly_ghost,'FaceColor','k','FaceAlpha',0.2,'EdgeColor',[1 1 1]*0.4); end

        end
    end
end

%plot([Floe(logical(cat(1,Floe.alive))).poly]);

if ~PERIODIC 
    xb=c2_boundary_poly.Vertices(:,1); xb(end+1)=xb(1);
    yb=c2_boundary_poly.Vertices(:,2); yb(end+1)=yb(1);
    plot(xb,yb, 'k-','linewidth',2);
end
[Nx,Ny] = size(squeeze(eularian_data.u));
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = flipud(y);
uscale = max(max(sqrt(ocean.Uocn.^2+ocean.Vocn.^2)));
xc = (x(2:end)+x(1:end-1))/2; yc = (y(2:end)+y(1:end-1))/2;
[xx,yy] = meshgrid(xc,yc);
Uocn_interp=interp2(ocean.Xocn,ocean.Yocn, ocean.Uocn,xx,yy);
Vocn_interp=interp2(ocean.Xocn,ocean.Yocn, ocean.Vocn,xx,yy);
quiver(xc,yc,25/uscale*Uocn_interp,25/uscale*Vocn_interp,'r','linewidth',2,'autoscale','off');
quiver(xc,yc,25/uscale*flipud(eularian_data.u),25/uscale*flipud(eularian_data.v),'g','linewidth',1,'autoscale','off');
% xx = 1; xx(1) =[1 2];


%plot(xb,yb,'-','linewidth',2,'color', 'k');
%plot([min(ocean.Xo)  max(ocean.Xo)], [min(yb) min(yb)],'-','linewidth',2,'color', 'k')
%plot([min(ocean.Xo)  max(ocean.Xo)], [max(yb) max(yb)],'-','linewidth',2,'color', 'k')

if showCenterOfMass
    plot(cat(1,Floe.Xi),cat(1,Floe.Yi),'k.'); 
end

if showContactPoints
    a=cat(1,Floe(logical(cat(1,Floe.alive))).interactions);
    if ~isempty(a)
        plot(a(:,4),a(:,5),'ro','markersize',2,'MarkerFaceColor','r','lineWidth',1);
    end
end

colorbar
colormap('gray'); caxis([0 1]);
axis([-Lx-Lx/10 Lx+Lx/10 -Lx-Lx/10 Lx+Lx/10])
%axis([min(ocean.Xo) max(ocean.Xo) min(ocean.Yo) max(ocean.Yo)])
xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow;

sacked_floes=sum(~cat(1,Floe.alive));
if sacked_floes>0, display(['Sacked floes: ' num2str(sacked_floes)]); end


end