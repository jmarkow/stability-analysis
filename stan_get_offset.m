function [SIG2_SHIFT]=stan_get_offset(SIG1,SIG2,varargin)
%
%
%
%

audio_proc=1;
song_band=[2e3 7e3];
rms_tau=.002;
fs=20e3;

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'audio_proc'
			audio_proc=varargin{i+1};
		case 'song_band'
			song_band=varargin{i+1};
		case 'rms_tau'
			rms_tau=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
	end
end

rms_smps=round(rms_tau*fs);
rms_filt=ones(rms_smps,1)/rms_smps;

[b,a]=ellip(3,.2,40,[song_band]/(fs/2),'bandpass');

if audio_proc
	pow1=filtfilt(b,a,SIG1).^2;
	pow2=filtfilt(b,a,SIG2).^2;
	rms1=sqrt(filter(rms_filt,1,pow1));
	rms2=sqrt(filter(rms_filt,1,pow2));
else
	rms1=SIG1;
	rms2=SIG2;
end

[c,lags]=xcorr(zscore(rms1),zscore(rms2));
[~,loc]=max(c);
lag2=lags(loc);

SIG2_SHIFT=lag2;
