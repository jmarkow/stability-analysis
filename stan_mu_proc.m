function [RMS,SPIKES]=stan_mu_proc(DATA,FS,varargin)

if nargin<2, PARAMETER_FILE=[]; end

if ~isa(DATA,'double')
	DATA=double(DATA);
end

SPIKES=[];
RMS=[];

if isempty(DATA)
	warning('No data found...');
	return;
end

tau=.005;
nmads=10;
sigma_t=2.5;
gauss_sd=.005;

param_names=who('-regexp','^[a-z]');

% parameter are all values that begin with lowercase letters

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

% scan for intan_frontend files, prefix songdet1

for i=1:2:nparams
	switch lower(varargin{i})
		case 'colors'
			colors=varargin{i+1};
		case 'blanking'
			blanking=varargin{i+1};
		case 'tau'
			tau=varargin{i+1};
		case 'sigma_t'
			sigma_t=varargin{i+1};
		case 'nmads'
			nmads=varargin{i+1};
		case 'gauss_sd'
			gauss_sd=varargin{i+1};
	end
end

% re-filter data?

[b,a]=ellip(3,.2,40,[600 4e3]/(FS/2),'bandpass');
DATA=filtfilt(b,a,DATA);

% rms

tau_smps=round(tau*FS);
rms_filter=ones(tau_smps,1)/tau_smps;

RMS.data=sqrt(filter(rms_filter,1,DATA.^2));
RMS.fs=FS;
RMS.params.tau_smps=tau_smps;

spikethreshold=sigma_t*median(abs(DATA)/.6745);
SPIKES=spikoclust_spike_detect_mu(DATA,spikethreshold,FS,'visualize','n','method','b');

% estimate smooth rate at lower sampling rate, convolve with standard kernel (normal, 5 ms variance)

SPIKES.smooth_rate=zeros(size(DATA));
idx=sub2ind(size(DATA),SPIKES.times,SPIKES.trial);
SPIKES.smooth_rate(idx)=1;

% probability per sample, convert to Hz

smooth_filter=normpdf([-gauss_sd*6:1/FS:gauss_sd*6],0,gauss_sd);
smooth_filter=smooth_filter./sum(smooth_filter);

SPIKES.smooth_rate=(filter(smooth_filter,1,SPIKES.smooth_rate))*FS;
SPIKES.smooth_params.filt=smooth_filter;
SPIKES.smooth_params.sd=gauss_sd;
