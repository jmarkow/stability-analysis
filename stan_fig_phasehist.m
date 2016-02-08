function stan_fig_phasehist()
%
%
%

[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.fig_dir,'ephys_lfp_spikes_phasehist.mat'),'phasehist');
fig=figure();
stan_ephys_plot_phasehist(phasehist);
set(fig,'units','centimeters','position',[4 4 3 7],'paperpositionmode','auto');
markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'phasehist_fig','eps,png,fig,pdf','renderer','painters');
