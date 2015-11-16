[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats.mat'),'stats');
fig=stan_plot_ca_nightday(stats);
names=fieldnames(fig);

for i=1:length(names)

	if strcmp(names{i},'daynightcompare')
		width=5.2;
	else
		width=4.5;
	end

	set(fig.(names{i}),'units','centimeters','position',[10 10 width 5],'paperpositionmode','auto');
  set(fig.(names{i}),'paperpositionmode','auto');

	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'cadata_analysis_' names{i} ],'eps,png,fig',...
		'renderer','painters');
end
