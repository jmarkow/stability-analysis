function BASELINE_STATS=stan_ephys_stats(EPHYS_DATA)
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collapse data into plot-able vectors

BASELINE_STATS.rms_corr={};
BASELINE_STATS.rms_mu={};
BASELINE_STATS.rms_var={};
BASELINE_STATS.spikes_corr={};
BASELINE_STATS.spikes_mu={};
BASELINE_STATS.spikes_var={};
BASELINE_STATS.threshold={};
BASELINE_STATS.days_since={};

padding_smps=round([options.padding-.1]*options.spike_fs);

mufun_corr=@(x) median(x,2);
mufun=@(x) mean(x);
varfun=@(x) std(x);

for i=1:length(EPHYS_DATA.dates)

	% over time plot:
	%
	%	(1) rate correlation
	%	(2) rms correlation
	%	(3) average rms
	%	(4) rms modulation (variance?)
	%	(5) threshold

	% get corr

	rms_mu_corr=cellfun(mufun_corr,EPHYS_DATA.rms{i},'uniformoutput',0);
	spikerate_mu_corr=cellfun(mufun_corr,EPHYS_DATA.spike_rate{i},'uniformoutput',0);

	rms_mu_corr=cat(2,rms_mu_corr{:});
	spikerate_mu_corr=cat(2,spikerate_mu_corr{:});

	rms_mu_corr=zscore(rms_mu_corr);
	spikerate_mu_corr=zscore(spikerate_mu_corr);

	BASELINE_STATS.rms_corr{i}=corr(rms_mu_corr);
	BASELINE_STATS.spikes_corr{i}=corr(spikerate_mu_corr);

	% average rms, rate, rms and rate modulation

	BASELINE_STATS.rms_mu{i}=cellfun(mufun,EPHYS_DATA.rms{i},'uniformoutput',0);
	BASELINE_STATS.spikes_mu{i}=cellfun(mufun,EPHYS_DATA.spike_rate{i},'uniformoutput',0);

	BASELINE_STATS.rms_var{i}=cellfun(varfun,EPHYS_DATA.rms{i},'uniformoutput',0);
	BASELINE_STATS.spikes_var{i}=cellfun(varfun,EPHYS_DATA.spike_rate{i},'uniformoutput',0);

	% get threshold

	BASELINE_STATS.threshold{i}=EPHYS_DATA.spike_threshold{i};
	BASELINE_STATS.days_since{i}=EPHYS_DATA.days_since{i};

end

save(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_baseline_stats.mat']),'BASELINE_STATS')
