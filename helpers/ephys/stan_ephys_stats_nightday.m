function store=stan_ephys_stats_nightday(EPHYS_DATA)
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collapse data into plot-able vectors

padding_smps=round([options.padding-.1]*options.spike_fs);

mufun_corr_day=@(x) mean(zscore(x(:,1:floor(size(x,2)/2))),2);
mufun_corr_night=@(x) mean(zscore(x(:,end-(floor(size(x,2)/2)-1):end)),2);
mufun_corr_all=@(x) mean(zscore(x),2);

%mufun_corr_day=@(x) mean(zscore(x(:,1:100)),2);
%mufun_corr_night=@(x) mean(zscore(x(:,end-100:end)),2);
mufun=@(x) mean(x);
varfun=@(x) std(x);
dist_mat={};

store.within_day=[];
store.between_day=[];
store.lags=[];

for i=1:length(EPHYS_DATA.dates)

	% over time plot:
	%
	%	(1) rate correlation
	%	(2) spike_rate correlation
	%	(3) average spike_rate
	%	(4) spike_rate modulation (variance?)
	%	(5) threshold

	% get corr

	disp([num2str(i)]);

	ntrials=cellfun(@(x) size(x,2),EPHYS_DATA.spike_rate{i});

	spikerate_mu_corr_day=cellfun(mufun_corr_day,EPHYS_DATA.spike_rate{i},'uniformoutput',0);
	spikerate_mu_corr_night=cellfun(mufun_corr_night,EPHYS_DATA.spike_rate{i},'uniformoutput',0);
	spikerate_mu_corr_all=cellfun(mufun_corr_all,EPHYS_DATA.spike_rate{i},'uniformoutput',0);
	spikerate_mu_corr_day=cat(2,spikerate_mu_corr_day{:});
	spikerate_mu_corr_night=cat(2,spikerate_mu_corr_night{:});
	spikerate_mu_corr_all=cat(2,spikerate_mu_corr_all{:});

	ndays=size(spikerate_mu_corr_day,2)
	daydiff=diff(EPHYS_DATA.days_since{i})
	lags=EPHYS_DATA.days_since{i}

	for j=1:ndays
		
		lags(j)

		for k=j:ndays
			curlag=lags(k)-lags(j);

			tmp_within=corr(spikerate_mu_corr_night(:,k),spikerate_mu_corr_night(:,j),'type','pearson');
			tmp_between=corr(spikerate_mu_corr_day(:,k),spikerate_mu_corr_night(:,j),'type','pearson');

			store.between_day=[store.between_day tmp_between];
			store.within_day=[store.within_day tmp_within];
			store.lags=[store.lags curlag];

		end

	end

end

%save(fullfile(dirs.agg_dir,dirs.datastore_dir,['mu_variance.mat']),'dist_mat');
