function phases=stan_ephys_stats_lfp_spikes(LFP_DATA)
%
%
%
%

[options,dirs]=stan_preflight;
% get n spikes

mintrials=inf;
smooth_fs=1e3;
smooth_sig=.005;
lfp_win=[.1 .1];

sigma=.005;

kernedges=[-3*sigma:1/smooth_fs:3*sigma];
kernel=(1/(sigma*sqrt(2*pi)))*exp((-(kernedges-0).^2)./(2*sigma^2));

for i=1:length(LFP_DATA)
	ntrials=size(LFP_DATA(i).lfp.data,2);

	if ntrials<mintrials
		mintrials=ntrials;
	end
end

phases.lfp_win=[];
phases.hil_win=[];
phases.cell_id=[];
phases.peak_id=[];
counter=1;

for i=1:length(LFP_DATA)
	
	disp([num2str(i)]);
	
	[b,a]=sfield_filt_coeffs(LFP_DATA(i).lfp.fs,2);
	lfp_data=filtfilt(b,a,LFP_DATA(i).lfp.data);

	% grab the phases.spikes associated with the spikes
	% convert spikes to appropriate sample

	hil_data=hilbert(lfp_data);
	[nsamples,ntrials]=size(LFP_DATA(i).lfp.data);

	mu=mean(zscore(LFP_DATA(i).spikes.smooth_rate)');
	%mu=mean(LFP_DATA(i).spikes.smooth_rate');
	[vals,locs]=findpeaks(mu,'minpeakheight',0,'minpeakdistance',round(.01*LFP_DATA(i).spikes.smooth_fs));
	
	%locs(locs<200|locs>nsamples-200)=[];
	locs=round((locs/smooth_fs)*LFP_DATA(i).lfp.fs);
	win_smps=round(lfp_win*LFP_DATA(i).lfp.fs);
	lfp_data=zscore(lfp_data);

	%lfp_data=zscore(lfp_data(200:end-200,:));
	%hil_data=hil_data(200:end-200,:);
	%nsamples=size(lfp_data,1);
	%locs=locs-200;
	locs(locs-win_smps(1)<1|locs>nsamples-win_smps(2))=[];

	for j=1:length(locs)
		phases.lfp_win=[phases.lfp_win lfp_data(locs(j)-win_smps(1):locs(j)+win_smps(2),:)];
		phases.hil_win=[phases.hil_win hil_data(locs(j)-win_smps(1):locs(j)+win_smps(2),:)];
		phases.cell_id=[phases.cell_id ones(1,size(lfp_data,2))*i];
		phases.peak_id=[phases.peak_id ones(1,size(lfp_data,2))*counter];
		counter=counter+1;
	end
end

%binvec=[-pi:pi/8:pi];
%[~,bins]=histc(phases.angles,binvec);

%phases.hist=zeros(1,length(unique(bins)));

%for i=1:length(bins)
%	phases.hist(bins(i))=phases.hist(bins(i))+phases.mags(i).^2;
%end

