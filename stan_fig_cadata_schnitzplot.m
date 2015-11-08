function stan_cadata_schnitzplot()
%
%
%
%
%

[options,dirs]=stan_preflight;
%load(fullfile(dirs.agg_dir,dirs.fig_dir,'AwesomeNess.mat'));
load(fullfile(dirs.agg_dir,dirs.ca_dir,'lw76.mat'),'roi_data');


fig.localnorm=figure();

stan_cadata_sortmat(roi_data,'scaling','l','sort_day',1,'smoothing',.15,'smooth_kernel','g',...
	'padding',[1 1],'movie_fs',22,'chk_day',1,'fig_row',1,'fig_nrows',2);

%stan_cadata_sortmat(roi_data,'scaling','l','sort_day',3,'smoothing',.15,'smooth_kernel','g',...
%	'padding',[1 1],'movie_fs',22,'chk_day',1,'fig_row',2,'fig_nrows',3)

stan_cadata_sortmat(roi_data,'scaling','l','sort_day',5,'smoothing',.15,'smooth_kernel','g',...
	'padding',[1 1],'movie_fs',22,'chk_day',1,'fig_row',2,'fig_nrows',2)
colormap(jet);

ax=findall(gcf,'type','axes')
linkaxes(ax,'xy')

ylimits=ylim();
set(ax(end),'YTick',ylim(),'YTickLabel',[ylim()-min(ylim())]+1,'FontSize',8);
yh=ylabel(ax(end),'Cell');
set(yh,'position',get(yh,'position')+[.3 0 0]);
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.4],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(1))
set(h,'clipping','off');

names=fieldnames(fig);

for i=1:length(names)
	set(fig.(names{i}),'units','centimeters','position',[10 10 16 5],'paperpositionmode','auto');
	%markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'schnitzplot_' names{i} ],'eps,png,fig,pdf');
end
