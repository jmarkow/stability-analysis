function [corrvals,comparemat,pmat]=stan_format_cadata(DATA,varargin)
% takes data from stan_format_cadata and generates a series of panels for each time point
%
%
%
%
%
%

movie_fs=22; % sampling rate of camera
upsample=10; % upsample factor (set to 1 for no upsampling)
upsample_method='spline'; % upsample method (spline and linear work fine)
peak_check=0; % check for peak consistency 
peak_thresh=.05; % if closest peak is >peak_thresh, exclude roi
chk_day=1; % check for dff peak day
scaling='r'; % scaling ('r' for within roi across days, 's' for within roi sort day, 'l' for within roi and day)
smoothing=0; % smooth ca trace (not working yet)
smooth_kernel='g'; % gauss smoothing kernel (b for boxcar)
padding=1; % padding before and after song
compare_day=1;
dff_check=.5;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'movie_fs'
			movie_fs=varargin{i+1};
		case 'upsample'
			upsample=varargin{i+1};
		case 'peak_check'
			peak_check=varargin{i+1};
		case 'peak_thresh'
			peak_thresh=varargin{i+1};
		case 'dff_check'
			dff_check=varargin{i+1};
		case 'scaling'
			scaling=varargin{i+1};
		case 'smoothing'
			smoothing=varargin{i+1};
        	case 'padding'
            		padding=varargin{i+1};
		case 'compare_day'
			compare_day=varargin{i+1};
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

if ~iscell(DATA)
	error('Wrong data format...');
end

ndays=length(DATA);
[nsamples,nrois,ntrials]=size(DATA{1});

[DATA,phase_shift]=stan_cadata_preprocess(DATA,'smoothing',smoothing,'smooth_kernel',smooth_kernel,'peak_check',peak_check,...
	'peak_thresh',peak_thresh,'movie_fs',movie_fs);

% get the sort indices

movie_x=[0:nsamples-1]/movie_fs;
pad_smps=padding*movie_fs;

if upsample>1

	for i=1:ndays
		interp_x=[0:1/upsample:nsamples-1]/movie_fs;
		DATA{i}=interp1(1:nsamples,DATA{i},[1:1/upsample:nsamples],'spline');
	end

	pad_smps=padding*(movie_fs*upsample);
	nsamples=size(DATA{1},1);
end


for i=1:ndays
	ave_mat{i}=mean(DATA{i},3);
end

% check for high enough dff

inc_rois=find(any(ave_mat{compare_day}>dff_check));
nrois=length(inc_rois);

for i=1:ndays
	DATA{i}=zscore(DATA{i});
end

% get peak locations, anything in the pads is removed


% null distribution for each day

ntrials=size(DATA{compare_day},3)
nrois=size(DATA{compare_day},2);

corrvals=zeros(nchoosek(ntrials,2),nrois);
for i=1:nrois
	corrmat=corr(squeeze(DATA{compare_day}(:,i,:)));
	corrvals(:,i)=corrmat(find(triu(ones(size(corrmat)),1)));
end

pmat=zeros(ndays-1,nrois);

for i=2:ndays
	for j=1:nrois
		corrmat=corr(squeeze(DATA{compare_day}(:,j,:)),squeeze(DATA{i}(:,j,:)));
		comparemat{i-1,j}=corrmat(find(triu(ones(size(corrmat)),1)));
		pmat(i-1,j)=ranksum(corrvals(:,j),comparemat{i-1,j},'tail','right','method','exact');
	end
end

% now run test between our comparison day and all comparison days

