function stan_fig_cadata_driftanalysis_all()
%
%
%
%

[options,dirs]=stan_preflight;
ext='con';
[fig,fig_stats]=stan_mu_ca_timecourse(ext);
names=fieldnames(fig);

for i=1:length(names)
	if length(strfind(names{i},'frac'))>0
		set(fig.(names{i}),'units','centimeters','position',[10 10 4.5 2],'paperpositionmode','auto');
    elseif length(strfind(names{i},'bee'))>0
		set(fig.(names{i}),'units','centimeters','position',[10 10 4.5 2],'paperpositionmode','auto');
    else
        set(fig.(names{i}),'units','centimeters','position',[10 10 4.5 4],'paperpositionmode','auto');
    end

  set(fig.(names{i}),'paperpositionmode','auto');
  markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'driftanalysis_revision-' ext '_' names{i} ],'eps,png,fig',...
 		'renderer','painters');
  save(fullfile(dirs.agg_dir,dirs.datastore_dir,['mu_ca_timecourse-' ext '.mat']),'fig_stats');
end
