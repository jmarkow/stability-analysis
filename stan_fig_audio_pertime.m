function stan_fig_cadata_driftanalysis_all()
%
%
%
%

[options,dirs]=stan_preflight;
[fig,fig_stats]=stan_audio_analyze_sdi_pertime;

%%
names=fieldnames(fig);
width=6;

for i=1:length(names)

	set(fig.(names{i}),'units','centimeters','position',[10 10 4.5 4],'paperpositionmode','auto');
  set(fig.(names{i}),'paperpositionmode','auto');

	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'audio_analysis_revision' names{i} ],'eps,png,fig',...
		'renderer','painters');
end
