function fig=plot_Floes(fig, Time,Floe,ocean,c2_boundary)

showContactPoints=0;
showCenterOfMass=0;

if (fig==0 || ~isvalid(fig))
    fig=figure('Position',[100 100 500 500],'visible','off');  
    set(fig,'PaperSize',[8 8],'PaperPosition',[0 0 8 8]);
end

clf(fig);

dn=2; % plot every dn'th velocity vector
quiver(ocean.Xo(1:dn:end),ocean.Yo(1:dn:end),ocean.Uocn(1:dn:end,1:dn:end),ocean.Vocn(1:dn:end,1:dn:end));
hold on;
axis([-1 1 -1 1]*7e4);
%axis([nanmin(c2_boundary(1,:)) nanmax(c2_boundary(1,:)) nanmin(c2_boundary(2,:)) nanmax(c2_boundary(2,:))]);
colormap('gray'); caxis([0 1]);
title(['Time = ' num2str(Time) ' s']);

for j=1:length(Floe)
    if Floe(j).alive
         c=[Floe(j).c_alpha(1,:)+Floe(j).Xi; Floe(j).c_alpha(2,:)+Floe(j).Yi];
         ind=~isnan(c(1,:));
         fill(c(1,ind),c(2,ind),[1 1 1]*0.6);
         axis ij
    end
end

plot(c2_boundary(1,:),c2_boundary(2,:),'linewidth',2,'color',[0 102 51]/255);

if showCenterOfMass
    plot(cat(1,Floe.Xi),cat(1,Floe.Yi),'k.'); 
end

if showContactPoints
    a=cat(1,Floe(logical(cat(1,Floe.alive))).interactions);
    if ~isempty(a)
        plot(a(:,4),a(:,5),'ro','markersize',2,'MarkerFaceColor','r','lineWidth',1);
    end
end

colormap('gray'); caxis([0 1]);
xlabel('m');ylabel('m');
set(gca,'Ydir','normal');

drawnow;

sacked_floes=sum(~cat(1,Floe.alive));
if sacked_floes>0, display(['Sacked floes: ' num2str(sacked_floes)]); end


end
