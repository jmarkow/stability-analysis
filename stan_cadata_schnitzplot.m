function stan_cadata_schnitzplot()
%
%
%
%
%

[options,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.fig_dir,'AwesomeNess.mat'));

newdata=stan_cadata_format(data1,data2,data3,data4,data5);

fig.sortnorm=figure();
stan_cadata_sortmat(newdata,'scaling','s','sort_day',1,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6);
ax=findall(gcf,'type','axes')
linkaxes(ax,'xy')
xlim([1.4 3])
ylim([4.5 max(ylim())]);
ylimits=ylim();
set(ax(end),'YTick',ylim(),'YTickLabel',[ylim()-min(ylim())]+1,'FontSize',8);
yh=ylabel(ax(end),'Cell');
set(yh,'position',get(yh,'position')+[.3 0 0]);
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.4],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(1))
set(h,'clipping','off');
title(ax(end-2),'Sort day scaling');


fig.roinorm=figure();
stan_cadata_sortmat(newdata,'scaling','r','sort_day',1,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6);
ax=findall(gcf,'type','axes')
linkaxes(ax,'xy')
xlim([1.4 3])
ylim([4.5 max(ylim())]);
ylimits=ylim();
set(ax(end),'YTick',ylim(),'YTickLabel',[ylim()-min(ylim())]+1,'FontSize',8);
yh=ylabel(ax(end),'Cell');
set(yh,'position',get(yh,'position')+[.3 0 0]);
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.4],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(1))
set(h,'clipping','off');
title(ax(end-2),'Across day scaling');

names=fieldnames(fig);

for i=1:length(names)
	set(fig.(names{i}),'units','centimeters','position',[10 10 16 5],'paperpositionmode','auto');
	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'schnitzplot_' names{i} ],'eps,png,fig,pdf');
end
