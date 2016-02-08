function LFP_DATA=stan_ephys_getdata_lfp_spikes_legacy(DIR)
%
%
%
%
% short script to estimate the coherence between spikes and LFPs collapsed across time
%
% methodology;
%
% 1) Compute multi-taper coherence
% 2) Compute either asymptotic or jackknife significance (jackknife will be much more time-consuming,
%    worth doing to validate distributional assumptions)

if nargin<1 | isempty(DIR), DIR=pwd; end

[options,dirs]=stan_preflight;
options.lfp_fs=1e3;

kernedges=[-3*options.smooth_sig:1/options.smooth_fs:3*options.smooth_sig];
kernel=(1/(options.smooth_sig*sqrt(2*pi)))*exp((-(kernedges-0).^2)./(2*options.smooth_sig^2));

[status,result]=unix(['find ' DIR ' -type f -name "spikedata_ch*.mat"']);
pkfile=regexp(result,'\n','split');

for i=1:length(pkfile)-1

	[path,file,ext]=fileparts(pkfile{i});

	fid=fopen(fullfile(path,'cellinfo.txt'),'r');

	readdata=textscan(fid,'%s%[^\n]','commentstyle','#',...
		'delimiter','\t','MultipleDelimsAsOne',1);	

	% close the file

	fclose(fid);

	% read in cluster number from cellinfo.txt

	clusternum=str2num(readdata{2}{find(strcmpi(readdata{1},'cluster:'))});
	channel=str2num(readdata{2}{find(strcmpi(readdata{1},'channel:'))});

	disp([pkfile{i}]);

	load(pkfile{i},'cluster');
	load(fullfile(path,'lfpdata.mat'),'EPHYS_DATA','CHANNELS','fs');

	ch_idx=CHANNELS==channel;

	LFP_DATA(i).spikes.times=cluster.times{clusternum};
	LFP_DATA(i).spikes.trial=cluster.trials{clusternum};
	LFP_DATA(i).spikes.fs=cluster.parameters.fs;
	LFP_DATA(i).spikes.thresholds=cluster.parameters.threshold;

	[nsamples,ntrials,~]=size(EPHYS_DATA);

	spikes=LFP_DATA(i).spikes.times;
	spikes=spikes/LFP_DATA(i).spikes.fs;
	trials=LFP_DATA(i).spikes.trial;

	uniq_trials=unique(trials);

	smooth_samples=ceil((nsamples/fs)*options.smooth_fs);
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
	
	downfact=round(fs/options.lfp_fs);
	[b,a]=ellip(4,.2,40,(options.lfp_fs/2)/(fs/2),'low');

	lfp_data=double(EPHYS_DATA(:,:,ch_idx));
	lfp_data=medfilt1(lfp_data,round(.001*fs));
	lfp_data=filtfilt(b,a,lfp_data);
	lfp_data=downsample(lfp_data,downfact);

	LFP_DATA(i).lfp.data=lfp_data;
	LFP_DATA(i).lfp.fs=options.lfp_fs;
	LFP_DATA(i).lfp.filt_b=b;
	LFP_DATA(i).lfp.filt_a=a;
	
end
