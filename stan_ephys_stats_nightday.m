function dist_mat=stan_ephys_stats_mu_variance(EPHYS_DATA)
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collapse data into plot-able vectors

padding_smps=round([options.padding-.1]*options.spike_fs);

mufun_corr_day=@(x) mean(zscore(x(:,1:floor(size(x,2)/2))),2);
mufun_corr_night=@(x) mean(zscore(x(:,end-(floor(size(x,2)/2)-1):end)),2);
%mufun_corr_day=@(x) mean(zscore(x(:,1:100)),2);
%mufun_corr_night=@(x) mean(zscore(x(:,end-100:end)),2);
mufun=@(x) mean(x);
varfun=@(x) std(x);
dist_mat={};

within_day=[];
between_day=[];

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

	spikerate_mu_corr_day=cat(2,spikerate_mu_corr_day{:});
	spikerate_mu_corr_night=cat(2,spikerate_mu_corr_night{:});

	ndays=size(spikerate_mu_corr_day,2);

	daydiff=diff(EPHYS_DATA.days_since{i});

	tmp_within=corr(spikerate_mu_corr_day,spikerate_mu_corr_night,'type','pearson');
	tmp_within=tmp_within(find(diag(ones(ndays,1),0)));

	idx=find(daydiff==1);

	if isempty(idx)
		continue;
	end

	tmp_between=corr(spikerate_mu_corr_night(:,idx),spikerate_mu_corr_day(:,idx+1),'type','pearson');
	size(tmp_between)
	length(idx)
	tmp_between=tmp_between(find(diag(ones(length(idx),1),0)));

	%tmp_between=tmp_between-(tmp_within(idx));
	%tmp_within=tmp_within(2:end)-tmp_within(1:end-1);

	within_day=[within_day;tmp_within(:)];
	between_day=[between_day;tmp_between(:)];

end

d{1}=within_day;
d{2}=between_day;
figure();plotSpread(d,'binwidth',.01);
figure();markolab_boxplot(d);
ranksum(within_day,between_day)

%save(fullfile(dirs.agg_dir,dirs.datastore_dir,['mu_variance.mat']),'dist_mat');
