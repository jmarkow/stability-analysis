function stan_ephys_getdata_lfp_phasehist_all
%
%
%
%

[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_lfp_spikes_data.mat']),'lfp_spikes');

cell_types=fieldnames(lfp_spikes);
for i=1:length(cell_types)

	if strcmp(cell_types{i},'mu')
		[phasehist.(cell_types{i}).mag phasehist.(cell_types{i}).bins]=stan_ephys_stats_phasehist(lfp_spikes.(cell_types{i}),1,0);
	else
		[phasehist.(cell_types{i}).mag phasehist.(cell_types{i}).bins]=stan_ephys_stats_phasehist(lfp_spikes.(cell_types{i}),0,100);
	end

end

save(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_lfp_spikes_phasehist.mat']),'phasehist');
