function stan_fig_barecarbon_timecourses()
%
%
%
%

[options,dirs]=stan_preflight;

load custom_colormaps;
load(fullfile(dirs.agg_dir,dirs.datastore_dir,'ephys_baseline_stats.mat'),'BASELINE_STATS');
fig=stan_ephys_plot_timecourse(BASELINE_STATS);
set(fig,'units','centimeters','position',[4 4 5 6],'paperpositionmode','auto');
ax=findall(fig,'type','axes');

for i=1:length(ax)
	set(ax(i),'fontsize',7);
	set(get(ax(i),'ylabel'),'fontsize',7);
	set(get(ax(i),'xlabel'),'fontsize',7);
end

markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_2c'],'eps,png,fig,pdf');
