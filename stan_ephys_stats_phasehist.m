function [phasehist,phasebins]=stan_ephys_stats_lfp_mu(LFP_DATA,NORMALIZE,THRESH)
%
%
%
%

if nargin<3 | isempty(THRESH)
	THRESH=100;
end

if nargin<2 | isempty(NORMALIZE)
	NORMALIZE=0;
end

[options,dirs]=stan_preflight;
phasebins=[-pi:pi/4:pi];
phasehist=zeros(1,length(phasebins));
% get n spikes

for i=1:length(LFP_DATA)
	
	disp([num2str(i)]);
	
	[b,a]=sfield_filt_coeffs(LFP_DATA(i).lfp.fs,2);
	lfp_data=filtfilt(b,a,LFP_DATA(i).lfp.data);

	% grab the phases.spikes associated with the spikes
	% convert spikes to appropriate sample

	hil_data=hilbert(mean(lfp_data'));
	ang_data=angle(hil_data);
	mag_data=abs(hil_data).^2;
	[nsamples,ntrials]=size(LFP_DATA(i).lfp.data);

	if NORMALIZE
		LFP_DATA(i).spikes.smooth_rate=zscore(LFP_DATA(i).spikes.smooth_rate);
	end

	mu=mean(LFP_DATA(i).spikes.smooth_rate');
	[vals,locs]=findpeaks(mu,'minpeakheight',THRESH,'minpeakdistance',round(.01*LFP_DATA(i).spikes.smooth_fs));
	
	locs(locs<200|locs>nsamples-200)=[];
	locs=round((locs/LFP_DATA(i).spikes.smooth_fs)*LFP_DATA(i).lfp.fs);
	
	%locs(locs-win_smps(1)<1|locs>nsamples-win_smps(2))=[];

	for j=1:length(locs)

		% get the lfp data at the spike time

		phases=ang_data(locs(j));
		mags=mag_data(locs(j));
		[n,bins]=histc(phases,phasebins);
		phasehist(bins)=phasehist(bins)+mags;

	end
end

