function stan_cadata_schnitzplot()
%
%
%
%
%

[options,dirs]=stan_preflight;
%load(fullfile(dirs.agg_dir,dirs.fig_dir,'AwesomeNess.mat'));
load(fullfile(dirs.agg_dir,dirs.fig_dir,'combined_rois.mat'));

newdata=stan_cadata_format(data1,data2,data3,data4,data5);

fig.localnorm=figure();
stan_cadata_sortmat(newdata,'scaling','l','sort_day',1,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6);
colormap(jet);
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
title(ax(end-2),'Same day/roi scaling');

fig.localnorm_rev=figure();
stan_cadata_sortmat(newdata,'scaling','l','sort_day',5,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6,'chk_day',5);
colormap(jet);
ax=findall(gcf,'type','axes')
linkaxes(ax,'xy')
xlim([1.4 3])
%ylim([4.5 max(ylim())]);
ylimits=ylim();
set(ax(end),'YTick',ylim(),'YTickLabel',[ylim()-min(ylim())]+1,'FontSize',8);
yh=ylabel(ax(end),'Cell');
set(yh,'position',get(yh,'position')+[.3 0 0]);
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.4],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(1))
set(h,'clipping','off');
title(ax(end-2),'Same day/roi scaling (rev)');

fig.sortnorm=figure();
stan_cadata_sortmat(newdata,'scaling','s','sort_day',1,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6);
colormap(jet);
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

fig.sortnorm_rev=figure();
stan_cadata_sortmat(newdata,'scaling','s','sort_day',5,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6,'chk_day',5);
colormap(jet);
ax=findall(gcf,'type','axes')
linkaxes(ax,'xy')
xlim([1.4 3])
%ylim([4.5 max(ylim())]);
ylimits=ylim();
set(ax(end),'YTick',ylim(),'YTickLabel',[ylim()-min(ylim())]+1,'FontSize',8);
yh=ylabel(ax(end),'Cell');
set(yh,'position',get(yh,'position')+[.3 0 0]);
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.4],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(1))
set(h,'clipping','off');
title(ax(end-2),'Sort day scaling (rev)');

fig.roinorm=figure();
stan_cadata_sortmat(newdata,'scaling','r','sort_day',1,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6);
colormap(jet);
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

fig.roinorm_rev=figure();
stan_cadata_sortmat(newdata,'scaling','r','sort_day',5,'smoothing',.15,'smooth_kernel','g',...
	'dff_check',.5,'padding',.6,'chk_day',5);
colormap(jet);
ax=findall(gcf,'type','axes')
linkaxes(ax,'xy')
xlim([1.4 3])
%ylim([4.5 max(ylim())]);
ylimits=ylim();
set(ax(end),'YTick',ylim(),'YTickLabel',[ylim()-min(ylim())]+1,'FontSize',8);
yh=ylabel(ax(end),'Cell');
set(yh,'position',get(yh,'position')+[.3 0 0]);
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.4],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(1))
set(h,'clipping','off');
title(ax(end-2),'Across day scaling (rev)');

names=fieldnames(fig);

for i=1:length(names)
	set(fig.(names{i}),'units','centimeters','position',[10 10 16 5],'paperpositionmode','auto');
	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'schnitzplot_' names{i} ],'eps,png,fig,pdf');
end
