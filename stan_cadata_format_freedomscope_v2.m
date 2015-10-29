function [FORM_DATA,FORM_T]=stan_cadata_format_freedomscope_v2(CADATA,TIME,THRESH,THRESH2,NEWFS,MINT,MAXT,PADDING)
% Using Bill's newformat, use the following command:
%
% [form_data,t]=stan_cadata_format_freedomscope_v2(roi_ave.RAWdat,roi_ave.RawTime);
%
% takes data where cell arrays correspond to separate songs, rows to samples, and columns to rois
% and reformats for stan_cadata_sortmat
%
% returns a 3d matrix samples x rois x trials 

if nargin<8
	PADDING=[];
end

if nargin<7
	MAXT=[];
end

if nargin<6
	MINT=[];
end

if nargin<5 | isempty(NEWFS)
	%NEWFS=29.97; % sampling rate of interpolated data
	NEWFS=100;
end

if nargin<4 | isempty(THRESH2)
	THRESH2=30; % threshold on raw pixel values
end

if nargin<3 | isempty(THRESH)
	THRESH=50; % threshold on deriv of raw pixel values
end

if nargin<2
	error('Need at least calcium data and frame times to continue.');
end

% first concatenate to sensible format

[nrois,ntrials]=size(CADATA);

FORM_DATA=cell(1,ntrials);
old_t=cell(1,ntrials);

for i=1:ntrials
	FORM_DATA{i}=cat(1,CADATA{:,i})';
	old_t{i}=cat(1,TIME{:,i})';
end

% high diff or low raw values are used to exclude trials w/ camera off

if ~isempty(PADDING)
	for i=1:ntrials
		FORM_DATA{i}=FORM_DATA{i}(PADDING(1):end-PADDING(2),:);
	end
end

cam_change=cellfun(@(x) any(any(abs(diff(x))>THRESH)),FORM_DATA);
cam_off=cellfun(@(x) any(any(x<THRESH2)),FORM_DATA);

FORM_DATA(cam_change|cam_off)=[];
old_t(cam_change|cam_off)=[];

% skip this section until we've debugged everything

% interpolate and cat

if isempty(MAXT)
	maxtime=inf;

	for i=1:length(old_t)
		tmp=min(old_t{i}(end,:));

		if tmp<maxtime
			maxtime=tmp;
		end
	end
else
	maxtime=MAXT;
end

if isempty(MINT)
	mintime=-inf;

	for i=1:length(old_t)
		tmp=max(old_t{i}(1,:));

		if tmp>mintime
			mintime=tmp;
		end
	end
else
	mintime=MINT;
end


%maxtime=1.6; % debug with fixed time to start
%mintime=.033;
mintime=0;
newtime=[mintime:1/NEWFS:maxtime]';

to_del=[];

% spline interpolate with extrapolation

for i=1:length(FORM_DATA)
	FORM_DATA{i}=interp1(old_t{i}(:,1),FORM_DATA{i},newtime,'spline','extrap');
end

FORM_T=newtime;

% convert to dff

for i=1:length(FORM_DATA)
	FORM_DATA{i}=fluolab_detrend(FORM_DATA{i},'fs',NEWFS,'method','prctile','win',0,'per',12);
end

% cat to 3D matrix, samples x rois x trials

FORM_DATA=cat(3,FORM_DATA{:});
