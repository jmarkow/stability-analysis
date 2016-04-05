function [peak_check,peak_ispeak,inc_rois]=stan_cadata_drift_analyze(DATA,varargin)
% takes data from stan_format_cadata and generates a series of panels for each time point
%
%
%
%
%
%
pmat=[];

movie_fs=22; % sampling rate of camera
upsample=1; % upsample factor (set to 1 for no upsampling)
upsample_method='spline'; % upsample method (spline and linear work fine)
peak_check_pad=0; % check for peak consistency
peak_thresh=.05; % if closest peak is >peak_thresh, exclude roi
dff_check=1;
chk_day=1; % check for dff peak day
scaling='r'; % scaling ('r' for within roi across days, 's' for within roi sort day, 'l' for within roi and day)
smoothing=0; % smooth ca trace (not working yet)
smooth_kernel='g'; % gauss smoothing kernel (b for boxcar)
padding=1; % padding before and after song
compare_day=1; % day to use as basis for comparison
nperms=1e3; % as expected, permutation and ranksum give roughly the same answer
method='r'; % (r)anksum, (t)test, (p)ermutation (note that permutation is dog slow)
nparams=length(varargin);
tail='right';
maxlag=.02;
lag_corr=0;
realign=1;
nboots=1e4;
dist_thresh=.1;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'movie_fs'
			movie_fs=varargin{i+1};
		case 'upsample'
			upsample=varargin{i+1};
		case 'peak_check_pad'
			peak_check_pad=varargin{i+1};
		case 'peak_thresh'
			peak_thresh=varargin{i+1};
		case 'dff_check'
			dff_check=varargin{i+1};
		case 'chk_day'
			chk_day=varargin{i+1};
		case 'scaling'
			scaling=varargin{i+1};
		case 'smoothing'
			smoothing=varargin{i+1};
		case 'smooth_kernel'
			smooth_kernel=varargin{i+1};
		case 'compare_day'
			compare_day=varargin{i+1};
    case 'padding'
      padding=varargin{i+1};
		case 'tail'
			tail=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'lag_corr'
			lag_corr=varargin{i+1};
		case 'realign'
			realign=varargin{i+1};
		case 'maxlag'
			maxlag=varargin{i+1};
		case 'dist_thresh'
			dist_thresh=varargin{i+1};
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

if ~iscell(DATA)
	error('Wrong data format...');
end

ndays=length(DATA);
[DATA,phase_shift,inc_rois]=stan_cadata_preprocess(DATA,'peak_check_pad',peak_check_pad,'peak_thresh',peak_thresh,'movie_fs',movie_fs,...
	'smoothing',smoothing,'smooth_kernel',smooth_kernel,'padding',padding,'realign',realign,'maxlag',maxlag);

% get the sort indices

[nsamples,nrois,ntrials]=size(DATA{1});
pad_smps=round(padding*movie_fs);

if pad_smps(1)==0;
	pad_smps(1)=1;
end

% not necessary here, but left in just in case
% upsample ave_mat, get peak times, yadda yadda
if upsample>1

	for i=1:ndays
		size(DATA{i})
		interp_x=[0:1/upsample:nsamples-1]/movie_fs;
		DATA{i}=interp1(1:nsamples,DATA{i},[1:1/upsample:nsamples],'spline');
	end

	pad_smps=padding*(movie_fs*upsample);
	nsamples=size(DATA{1},1);
	phase_shift=round((phase_shift/movie_fs)*(movie_fs*upsample));

end
pad_smps
for i=1:ndays
	ave_mat{i}=mean(DATA{i}(pad_smps(1):end-pad_smps(2),:,:),3);
end


[template.peaks,template.vals]=fb_compute_peak_simple(ave_mat{compare_day},...
	'thresh_t',.1,'debug',0,'onset_only',0,'thresh_hi',.5,'thresh_int',5,'thresh_dist',.2,...
	'fs',movie_fs*upsample); % thresh_int previously


peak_check=zeros(ndays,nrois);
peak_ispeak=zeros(ndays,nrois);

for i=1:ndays
	[tmp.peaks,tmp.vals]=fb_compute_peak_simple(ave_mat{i},...
		'thresh_t',.1,'debug',0,'onset_only',0,'thresh_hi',.5,'thresh_int',5,'thresh_dist',.2,...
		'fs',movie_fs*upsample); % thresh_int previously
		for j=1:nrois
			% cycle through all peaks
			mindist=inf;
			for k=1:length(tmp.peaks{j})
					tmp2=min(abs(tmp.peaks{j}(k)-template.peaks{j}))/(movie_fs*upsample);
					if tmp2<mindist
						mindist=tmp2;
					end
			end

			% leave a flag to check if peak exists

			if mindist<dist_thresh
				peak_check(i,j)=1;
			end
			if length(tmp.peaks{j})>0
				peak_ispeak(i,j)=1;
			end
		end
end

% rifle through remaining days, is there an ROI within n msec of peak time
% take the shortest distance for the most conservative measure
% return whether we have a hit or miss, ndays x nrois, easy to compute stats on
