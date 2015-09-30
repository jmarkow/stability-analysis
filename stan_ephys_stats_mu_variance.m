function dist_mat=stan_ephys_stats_mu_variance(EPHYS_DATA)
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collapse data into plot-able vectors

padding_smps=round([options.padding-.1]*options.spike_fs);

mufun_corr=@(x) median(x,2);
mufun=@(x) mean(x);
varfun=@(x) std(x);
dist_mat={};

for i=1:length(EPHYS_DATA.dates)

	% over time plot:
	%
	%	(1) rate correlation
	%	(2) rms correlation
	%	(3) average rms
	%	(4) rms modulation (variance?)
	%	(5) threshold

	% get corr

	disp([num2str(i)]);

	%rms_mu_corr=cellfun(mufun_corr,EPHYS_DATA.rms{i},'uniformoutput',0);
	
	spikerate_mu_corr=cellfun(mufun_corr,EPHYS_DATA.spike_rate{i},'uniformoutput',0);

	% loop through spike_rates, estimate variance via peak-to-peak distance across trials 
	% zscore, within trials take all turning points, take closest peak relative to mean peaks

	% typecast, etc. etc.

	cur_rates=double(EPHYS_DATA.spike_rate{i}{1});
	template=(mean(cur_rates,2));

	[~,template_locs]=findpeaks(template,'minpeakheight',100,'minpeakdistance',round(.005*options.spike_fs));

	% loop, find closest peak 
	
	[nsamples,ntrials]=size(cur_rates);

	nlocs=length(template_locs);

	dist_mat{i}=zeros(nlocs,ntrials);

	for j=1:ntrials

		[~,trial_locs]=findpeaks((cur_rates(:,j)),'minpeakheight',100,'minpeakdistance',round(.005*options.spike_fs));

		for k=1:nlocs
			tmp=min(abs(trial_locs-template_locs(k)));
			if isempty(tmp)
				tmp=NaN;
			end
			dist_mat{i}(k,j)=tmp;
		end
	end

end
