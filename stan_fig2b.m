function stan_fig2b
%
%
%
%

[options,dirs]=stan_preflight;

load custom_colormaps;
load(fullfile(dirs.agg_dir,dirs.fig_dir,'ephys_baseline_stats.mat'),'BASELINE_STATS');
fig=stan_ephys_plot_correlation_baseline(BASELINE_STATS);
set(fig,'units','centimeters','position',[4 4 4.5 8],'paperpositionmode','auto');
