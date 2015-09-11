function LFP_DATA=stan_ephys_getdata_baseline()
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;
options.lfp_fs=1e3;

tmp=dir(fullfile(dirs.agg_dir,dirs.lfp_mu_dir,'*.mat'));
lfp_files={tmp(:).name};

LFP_DATA=struct();
rem_bird=[];

% filter setting

%[n,Wn,beta,ftype]=kaiserord([10 25 35 50],[0 1 0],[.01 .05 .01],options.lfp_fs);
%b=fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
%a=1;
%[b,a]=ellip(4,.2,40,[25 35]/(options.lfp_fs/2),'bandpass');

for i=1:length(lfp_files)

	disp([lfp_files{i}]);

	load(fullfile(dirs.agg_dir,dirs.lfp_mu_dir,lfp_files{i}),'adc');

	% get fields, process for spikes, etc. etc.

	ch_idx=cellfun(@length,strfind(lower(adc.names),'unfiltered'))>0;

	[b,a]=ellip(3,.2,40,[600 4e3]/(adc.fs/2),'bandpass');
	spike_data=filtfilt(b,a,double(adc.data(:,:,ch_idx)));
	
	spikethreshold=options.sigma_t*median(abs(spike_data)/.6745);
	LFP_DATA(i).spikes=spikoclust_spike_detect_mu(spike_data,spikethreshold,adc.fs,'visualize','n','method','b');

	% downsample for lfp processing

	% median filter spikes
	
	downfact=round(adc.fs/options.lfp_fs);
	[b,a]=ellip(4,.2,40,(options.lfp_fs/2)/(adc.fs/2),'low');

	lfp_data=double(adc.data(:,:,ch_idx));
	lfp_data=medfilt1(lfp_data,round(.0025*adc.fs));
	lfp_data=filtfilt(b,a,lfp_data);
	lfp_data=downsample(lfp_data,downfact);

	LFP_DATA(i).lfp.data=lfp_data;
	LFP_DATA(i).lfp.fs=options.lfp_fs;
	LFP_DATA(i).lfp.filt_b=b;
	LFP_DATA(i).lfp.filt_a=a;

end

save(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_lfp_mu_data.mat']),'LFP_DATA')
