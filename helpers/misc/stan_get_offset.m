function [SHIFT,ID]=stan_get_offset(SIG1,SIG2,varargin)
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
	pow{1}=filtfilt(b,a,SIG1).^2;
	pow{2}=filtfilt(b,a,SIG2).^2;
	rms{1}=sqrt(filter(rms_filt,1,pow{1}));
	rms{2}=sqrt(filter(rms_filt,1,pow{2}));
else
	rms{1}=SIG{1};
	rms{2}=SIG{2};
end

[c,lags]=xcorr(zscore(rms{1}),zscore(rms{2}));
[~,loc]=max(c);
lag2=lags(loc);

% which signal is smaller

len(1)=length(SIG1);
len(2)=length(SIG2);

[~,template]=min(len)
[~,compare]=max(len);

% slide template over other signal

step_size=1;

len1=len(template);
len2=len(compare);

template_rms=rms{template};
compare_rms=rms{compare};

steps=1:len2-len1;
dist=zeros(1,length(steps));

for i=1:length(steps)
	% euclidean distance match

	dist(i)=sum(abs(template_rms-compare_rms(i:i+len1-1)));
	
end

ID=compare;
[~,SHIFT]=min(dist)



