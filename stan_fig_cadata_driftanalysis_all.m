%
%
%
%

[options,dirs]=stan_preflight;
fig=stan_plot_mu_ca_timecourse();
names=fieldnames(fig);

for i=1:length(names)
	if length(strfind(names{i},'frac'))>0
		set(fig.(names{i}),'units','centimeters','position',[10 10 4.5 2],'paperpositionmode','auto');
	else
		set(fig.(names{i}),'units','centimeters','position',[10 10 4.5 4],'paperpositionmode','auto');
	end

  set(fig.(names{i}),'paperpositionmode','auto');
	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'driftanalysis_' names{i} ],'eps,png,fig',...
		'renderer','painters');
end
