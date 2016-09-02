function stan_fig_cadata_driftanalysis_all()
%
%
%
%

[options,dirs]=stan_preflight;
ext='con';
[fig,fig_stats]=stan_cadata_drift_analyze_pertime_ave(ext);

%%
names=fieldnames(fig);
width=6;

for i=1:length(names)
	set(fig.(names{i}),'units','centimeters','position',[10 10 4.5 4],'paperpositionmode','auto');
	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'cadata_analysis_revision-' ext '_' names{i} ],'eps,png,fig',...
		'renderer','painters');
end
