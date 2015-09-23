%
%
%
%

[options,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.fig_dir,'ephys_baseline_lfp_data.mat'));

fig=figure();stan_plot_lfp_angdiff(LFP_DATA);
ylim([0 .3])
xlim([0 pi])
set(gca,'XTick',[0:pi/2:pi],'XTickLabel',{'0','pi/2','pi'},'FontName','symbol','YTick',[0 .3]);
xlabel('Phase diff.','FontName','Helvetica')
ylabel('P','FontName','Helvetica')

set(fig,'units','centimeters','paperpositionmode','auto')
set(fig,'position',[5 5 5 4])

markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'lfp_phase_diff','eps,png,fig,pdf','renderer','painters');

% log-normal mle estimates, mode most intuitive here 


