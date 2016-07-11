function stan_fig_cadata_nightdayanalysis()
%
[options,dirs]=stan_preflight;

%load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats_new-lib.mat'),'stats');
ext='con';
[fig,overnight_stats]=stan_plot_ca_nightday(ext);
names=fieldnames(fig);

for i=1:length(names)

	if strcmp(names{i},'daynightcompare')
		width=4.5;
	elseif strfind(names{i},'overnight')
		width=2.8;
    else
        ylim([0 50]);
		width=4.5;
	end

	set(fig.(names{i}),'units','centimeters','position',[10 10 width 4],'paperpositionmode','auto');
    set(fig.(names{i}),'paperpositionmode','auto');

	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'cadata_analysis_revision-' ext '_' names{i} ],'eps,png,fig',...
		'renderer','painters');

end
