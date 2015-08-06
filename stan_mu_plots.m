function MUFIGS=fluolab_fb_plots(RMS,SPIKES,AUDIO,varargin)

blanking=[.2 .2];
facecolor=[0 0 1];
edgecolor=[0 0 1];
spect_colors='jet';
rms_colors='parula';
cbar_width = .025; 
cbar_dist=.01;
visible='on';
datenums=[];
padding=[.3 .3];
ch_names='';
rms_tau=.005;
song_band=[3e3 7e3];

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'facecolor'
			facecolor=varargin{i+1};
		case 'edgecolor'
			edgecolor=varargin{i+1};
		case 'visible'
			visible=varargin{i+1};
		case 'datenums'
			datenums=varargin{i+1};
		case 'padding'
			padding=varargin{i+1};
		case 'ch_names'
			ch_names=varargin{i+1};
		case 'rms_tau'
			rms_tau=varargin{i+1};
	end
end

MUFIGS=struct();

% take first catch or all trial for AUDIO sample


sample=AUDIO.data(:,1);

padding_smps=round(padding*AUDIO.fs);

sample(1:padding_smps(1))=0;
sample(end-padding_smps(2):end)=0;

[spect.s,spect.f,spect.t]=zftftb_pretty_sonogram(sample,...
	AUDIO.fs,'filtering',300,'clipping',[-2 2],...
	'len',80,'overlap',79,'zeropad',0,'norm_amp',1);
spect.f=spect.f/1e3;

name='';

if ~iscell(SPIKES)
	tmp{1}=SPIKES;
	clear SPIKES;
	SPIKES=tmp;
	clear tmp;

	tmp{1}=RMS;
	clear RMS;
	RMS=tmp;
	clear tmp;
end

for i=1:length(SPIKES)
	
	if ~isempty(ch_names) & iscell(ch_names)
		name=ch_names{i};
	elseif ~isempty(ch_names) & ischar(ch_names)
		name=ch_names;
	end

	if max(SPIKES{i}.trial)<2
		continue;
	end

	fig_name=['songalign_raster_ch' sprintf('%i',i) ];
	MUFIGS.(fig_name)=figure('paperpositionmode','auto','visible',visible);
	stan_songalign_raster(spect,SPIKES{i},'name',name,'colors',spect_colors);

	fig_name=['songalign_spikerate_ch' sprintf('%i',i) ];
	MUFIGS.(fig_name)=figure('paperpositionmode','auto','visible',visible);
	stan_songalign_rms(spect,SPIKES{i}.smooth_rate,SPIKES{i}.fs,'name',name,'spect_colors',spect_colors,'rms_colors',rms_colors,'clim_order',10,...
		'cbar_label','Firing Rate (Hz)','offset',length(SPIKES{i}.smooth_params.filt)/SPIKES{i}.fs);

end


for i=1:length(RMS)

	if ~isempty(ch_names) & iscell(ch_names)
		name=ch_names{i};
	elseif ~isempty(ch_names) & ischar(ch_names)
		name=ch_names;
	end

	if size(RMS{i}.data,2)<2
		continue;
	end

	fig_name=['songalign_rms_ch' sprintf('%i',i) ];
	MUFIGS.(fig_name)=figure('paperpositionmode','auto','visible',visible);	
	stan_songalign_rms(spect,RMS{i}.data,RMS{i}.fs,'name',name,'spect_colors',spect_colors,'rms_colors',rms_colors,...
		'clim_order',.01,'cbar_label','RMS (uV)','offset',RMS{i}.params.tau_smps/RMS{i}.fs);

end

if size(AUDIO.data,2)<2
	return;
end

% get song rms
fig_name='song_rms';
MUFIGS.(fig_name)=figure('paperpositionmode','auto','visible',visible);
[b,a]=ellip(3,.2,40,[song_band]/(AUDIO.fs/2),'bandpass');
AUDIO.data=filtfilt(b,a,double(AUDIO.data));

tau_smps=round(rms_tau*AUDIO.fs);
rms_filter=ones(1,tau_smps)/tau_smps;

AUDIO.data=sqrt(filter(rms_filter,1,AUDIO.data.^2));
stan_songalign_rms(spect,AUDIO.data,AUDIO.fs,'name','Song RMS','spect_colors',spect_colors,'rms_colors',rms_colors,'clim_order',.01,...
	'cbar_label','RMS (V)','offset',rms_tau);
