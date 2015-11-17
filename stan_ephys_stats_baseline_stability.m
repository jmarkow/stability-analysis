function [teststats]=stan_ephys_stats_baseline(BASELINE_DATA,varargin)
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

scaling='l';
nperms=1e3; % as expected, permutation and ranksum give roughly the same answer
method='r'; % (r)anksum, (t)test, (p)ermutation (note that permutation is dog slow)
nparams=length(varargin);
tail='right';
data_type='spike_rate';

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'scaling'
			scaling=varargin{i+1};
		case 'tail'
			tail=varargin{i+1};
		case 'method'
			method=varargin{i+1};
	end
end

% collapse data into plot-able vectors
% use all control data, estimate whether RMS or SPIKE RATE corr significantly changes

use_data=BASELINE_DATA.(data_type);
counter=1;

for i=1:length(use_data)

	% take uppper triangle of correlation matrix from first day, compare with
	% pairwise on subsequent days, etc. etc.

	ntrials=cellfun(@(x) size(x,2),use_data{i});
	trials_to_use=find(ntrials>options.ephys_trial_limit);

	if length(trials_to_use)<2
		continue;
	end

	tmp=corr(zscore(use_data{i}{trials_to_use(1)}));
	baseline_data=tmp(find(triu(ones(size(tmp)),1))); % get upper triangle of corr

	i

	for j=trials_to_use

		mu1=mean(zscore(use_data{i}{trials_to_use(j)}),2);

		for k=trials_to_use

			%tmp=corr(zscore(use_data{i}{trials_to_use(1)}),zscore(use_data{i}{j}));
			mu2=mean(zscore(use_data{i}{k}),2);
			tmp_mu=corr(mu1(:),mu2(:),'type','pearson');

			%new_data=tmp(find(triu(ones(size(tmp)),0))); % get upper triangle of corr
			%teststats.val(counter)=(mean(new_data(:))-mean(baseline_data(:)))./(sqrt(std(baseline_data(:))*std(new_data(:))));

			teststats.days_since(counter)=BASELINE_DATA.days_since{i}(k)-BASELINE_DATA.days_since{i}(j);
			teststats.birdid(counter)=i;
			teststats.val_mu(counter)=tmp_mu;

			%[teststats.pval{i}(counter),h]=ranksum(baseline_data(:),new_data(:),'tail','right')
			counter=counter+1;

		end

	end

end

save(fullfile(dirs.agg_dir,dirs.datastore_dir,['mu_baseline_stability']),'teststats');
