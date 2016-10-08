function stan_fig2ij(CORR_DATA)
% Generates Fig 2i-j, correlation time-course

fontsize=7;
ci_inv=.01; % alpha level
ylimit_rounding=1e-1;
markersize=5;

% spike and rms plot

fig=figure();

ax(1)=subplot(2,1,1);

b(1,:)=stan_plot_regress(CORR_DATA.days_since(:),CORR_DATA.spikes(:),CORR_DATA.birdid(:),'markersize',markersize);
yh=ylabel('FR');
set(yh,'position',get(yh,'position')+[5 -.2 0])
set(gca,'XTick',[],'TickLength',[0 0],'YTick',[0 1],'fontsize',fontsize);

pos=get(ax(1),'position')
asp_ratio=pos(3)/pos(4);
width=.3;
new_width=width/asp_ratio;
h_offset=.25;
v_offset=.27;

new_axis(1)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);
stan_plot_regress(CORR_DATA.days_since(:),CORR_DATA.spikes(:),CORR_DATA.birdid(:),'markersize',markersize);

ylim([0 1]);
set(gca,'XTick',[],'YTick',[],'TickDir','out','TickLength',[ 0 0 ],'fontsize',fontsize);

ax(2)=subplot(2,1,2);
b(2,:)=stan_plot_regress(CORR_DATA.days_since(:),CORR_DATA.rms(:),CORR_DATA.birdid(:),'markersize',markersize);
yh=ylabel('RMS');
set(yh,'position',get(yh,'position')+[5 -.2 0])
set(gca,'TickLength',[0 0],'YTick',[0 1],'fontsize',fontsize);
ylim([0 1]);

pos=get(ax(2),'position')

new_axis(2)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);

stan_plot_regress(CORR_DATA.days_since(:),CORR_DATA.rms(:),CORR_DATA.birdid(:),'markersize',markersize);
ylim([0 1]);
set(gca,'XTick',[],'YTick',[],'TickDir','out','TickLength',[ 0 0 ],'fontsize',fontsize);
[r,p]=corrcoef(CORR_DATA.days_since(:),CORR_DATA.rms(:))
linkaxes(ax,'xy');
linkaxes(new_axis,'xy');
linkaxes([ax(:);new_axis(:)],'y');
ylim([0 1]);

linkaxes(ax,'x');
set(ax(1),'xlim',[0 40]);
xh=xlabel(ax(end),'Days');
set(ax(end),'xtick',get(ax(end),'xlim'));
set(xh,'position',get(xh,'position')+[0 .01 0]);

linkaxes(new_axis,'x');
set(new_axis(1),'xlim',[0 120]);
set(new_axis(end),'xtick',get(new_axis(end),'xlim'));

set(fig,'units','centimeters','position',[4 4 4 6],'paperpositionmode','auto');
ax=findall(fig,'type','axes');

for i=1:length(ax)
	set(ax(i),'fontsize',7);
	set(get(ax(i),'ylabel'),'fontsize',7);
	set(get(ax(i),'xlabel'),'fontsize',7);
end
