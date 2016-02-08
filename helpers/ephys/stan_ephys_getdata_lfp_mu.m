function LFP_DATA=stan_ephys_getdata_lfp_mu
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;
options.lfp_fs=1e3;

kernedges=[-3*options.smooth_sig:1/options.smooth_fs:3*options.smooth_sig];
kernel=(1/(options.smooth_sig*sqrt(2*pi)))*exp((-(kernedges-0).^2)./(2*options.smooth_sig^2));

tmp=dir(fullfile(dirs.agg_dir,dirs.lfp_mu_dir,'*.mat'));
lfp_files={tmp(:).name};

LFP_DATA=struct();
rem_bird=[];

% filter setting

for i=1:length(lfp_files)

	disp([lfp_files{i}]);

	load(fullfile(dirs.agg_dir,dirs.lfp_mu_dir,lfp_files{i}),'adc');

	% get fields, process for spikes, etc. etc.

	ch_idx=cellfun(@length,strfind(lower(adc.names),'unfiltered'))>0;

	[b,a]=ellip(3,.2,40,[600 4e3]/(adc.fs/2),'bandpass');
	spike_data=filtfilt(b,a,double(adc.data(:,:,ch_idx)));
	
	spikethreshold=options.sigma_t*median(abs(spike_data)/.6745);
	LFP_DATA(i).spikes=spikoclust_spike_detect_mu(spike_data,spikethreshold,adc.fs,'visualize','n','method','b');

	[nsamples,ntrials]=size(spike_data);

	spikes=LFP_DATA(i).spikes.times;
	spikes=spikes/LFP_DATA(i).spikes.fs;
	trials=LFP_DATA(i).spikes.trial;

	uniq_trials=unique(trials);

	smooth_samples=ceil((nsamples/adc.fs)*options.smooth_fs);
	smooth_rate=zeros(smooth_samples,ntrials);

	for j=1:length(uniq_trials)
		spikes_smps=round(spikes(trials==uniq_trials(j))*options.smooth_fs);
		spikes_smps(spikes_smps==0)=1;
		smooth_rate(spikes_smps,uniq_trials(j))=1;
		smooth_rate(:,uniq_trials(j))=conv(smooth_rate(:,j),kernel,'same');
	end

	LFP_DATA(i).spikes.smooth_rate=smooth_rate;
	LFP_DATA(i).spikes.smooth_fs=options.smooth_fs;
	LFP_DATA(i).spikes.smooth_sig=options.smooth_sig;

	% downsample for lfp processing
	% median filter spikes
	
	downfact=round(adc.fs/options.lfp_fs);
	[b,a]=ellip(4,.2,40,(options.lfp_fs/2)/(adc.fs/2),'low');

	lfp_data=double(adc.data(:,:,ch_idx));
	lfp_data=medfilt1(lfp_data,round(.001*adc.fs));
	lfp_data=filtfilt(b,a,lfp_data);
	lfp_data=downsample(lfp_data,downfact);

	LFP_DATA(i).lfp.data=lfp_data;
	LFP_DATA(i).lfp.fs=options.lfp_fs;
	LFP_DATA(i).lfp.filt_b=b;
	LFP_DATA(i).lfp.filt_a=a;

end
