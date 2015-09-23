function stan_ephys_getdata_lfp_spikes_all
%
%
%
%

[options,dirs]=stan_preflight;

lfp_spikes.mu=stan_ephys_getdata_lfp_mu;
lfp_spikes.int=stan_ephys_getdata_lfp_spikes_legacy(fullfile(dirs.legacy_dir,dirs.legacy_int));
lfp_spikes.pn=stan_ephys_getdata_lfp_spikes_legacy(fullfile(dirs.legacy_dir,dirs.legacy_pn));

save(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_lfp_spikes_data.mat']),'lfp_spikes')
