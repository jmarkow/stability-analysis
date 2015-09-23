function stan_fig1a
%
%
%

[options,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.fig_dir,'btlfp_data.mat'));
fig=figure();
stan_ephys_plot_btlfp(BTLFP);
set(fig,'units','centimeters','position',[4 4 3 7],'paperpositionmode','auto');
markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'btlfp_fig','eps,png,fig,pdf','renderer','painters');

