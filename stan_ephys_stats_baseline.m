function STATS=stan_ephys_stats(EPHYS_DATA)
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collapse data into plot-able vectors

STATS.rms_corr={};
STATS.rms_mu={};
STATS.rms_var={};
STATS.spikes_corr={};
STATS.spikes_mu={};
STATS.spikes_var={};
STATS.threshold={};
STATS.days_since={};

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

	STATS.rms_corr{i}=corr(rms_mu_corr);
	STATS.spikes_corr{i}=corr(spikerate_mu_corr);

	% average rms, rate, rms and rate modulation

	STATS.rms_mu{i}=cellfun(mufun,EPHYS_DATA.rms{i},'uniformoutput',0);
	STATS.spikes_mu{i}=cellfun(mufun,EPHYS_DATA.spike_rate{i},'uniformoutput',0);

	STATS.rms_var{i}=cellfun(varfun,EPHYS_DATA.rms{i},'uniformoutput',0);
	STATS.spikes_var{i}=cellfun(varfun,EPHYS_DATA.spike_rate{i},'uniformoutput',0);

	% get threshold

	STATS.threshold{i}=EPHYS_DATA.spike_threshold{i};
	STATS.days_since{i}=EPHYS_DATA.days_since{i};

end
