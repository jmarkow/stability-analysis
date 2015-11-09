function stan_fig3b
% REQUIRES plotSpread for beeswarm plot
%
%
%

[options,dirs]=stan_preflight;

post_z=stan_audio_analyze_nervecut;
[fig,p]=stan_audio_nervecut_plot(post_z);

fig_names=fieldnames(fig);

set(fig.swarm,'units','centimeters','position',[4 4 4.5 6],'paperpositionmode','auto');
set(fig.box_plot,'units','centimeters','position',[4 4 4.5 4],'paperpositionmode','auto');

markolab_multi_fig_save(fig.swarm,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_3b_swarm'],'eps,png,fig,pdf','renderer','painters');
markolab_multi_fig_save(fig.box_plot,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_3b_boxplot'],'eps,png,fig,pdf');
