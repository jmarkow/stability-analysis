% REQUIRES plotSpread for beeswarm plot
%
%
%

[options,dirs]=stan_preflight;
load custom_colormaps;
load(fullfile(dirs.agg_dir,dirs.fig_dir,'ephys_baseline_stats.mat'),'BASELINE_STATS');
load(fullfile(dirs.agg_dir,dirs.fig_dir,'ephys_nervecut_stats.mat'),'NERVECUT_STATS');
fig=stan_ephys_plot_correlation_comparison_smooth(BASELINE_STATS,NERVECUT_STATS);
set(fig,'units','centimeters','position',[5 5 4.6 6.4],'paperpositionmode','auto');
markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'figure_4b','eps,png,fig,pdf');

