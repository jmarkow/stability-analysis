function [FORM_DATA,FORM_T,FORM_DATE,FORM_MOTIF,FORM_TRIALS,FILENAMES]=stan_cadata_format_freedomscope_v2(...
	CADATA,TIME,THRESH,THRESH2,NEWFS,PADDING,SONG_LEN,OFFSET,FILENAMES,METHOD)
% Using Bill's new format, use the following command:
%
% [form_data,t]=stan_cadata_format_freedomscope_v2(roi_ave.RAWdat,roi_ave.RawTime);
%
% takes data where cell arrays correspond to separate songs, rows to samples, and columns to rois
% and reformats for stan_cadata_sortmat
%
% returns a 3d matrix samples x rois x trials

if nargin<10
	METHOD=1;
end

if nargin<9
  FILENAMES=[];
end

if nargin<8
	OFFSET=[];
end

if nargin<7
	SONG_LEN=.7;
end

if nargin<6
	PADDING=[];
end

if nargin<5 | isempty(NEWFS)
	NEWFS=29.97; % sampling rate of interpolated data
	%NEWFS=100;
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
sig_t=6;
[nrois,ntrials]=size(CADATA);

FORM_DATA=cell(1,ntrials);
FORM_DATE=zeros(1,ntrials);
FORM_MOTIF=zeros(1,ntrials);
FORM_TRIALS=1:ntrials;
old_t=cell(1,ntrials);

for i=1:ntrials
	FORM_DATA{i}=cat(1,CADATA{:,i})';
	old_t{i}=cat(1,TIME{:,i})';
	if ~isempty(OFFSET)
		offset=cat(1,OFFSET{:})';
	end
	if ~isempty(FILENAMES)
		if iscell(FILENAMES{i})
			FILENAMES{i}=FILENAMES{i}{1};
		end
		[~,filename,ext]=fileparts(FILENAMES{i});
		tmp=regexp(filename,'\d+-\d+-\d+ \d+ \d+ \d+','match');
		FORM_DATE(i)=datenum(tmp{1},'yyyy-mm-dd HH MM SS');
		tmp=regexp(filename,'\_(\d+)','tokens');
		FORM_MOTIF(i)=str2num(tmp{1}{1});
	end
end

% time points closest to pad

to_del=[];

for i=1:length(FORM_DATA)

	% apply offset correction if it exists

	if ~isempty(OFFSET)
		old_t{i}=old_t{i}-offset(i)/1e3;
	end

	% find the timepoint just past the pad, and the timepoint just to the left of the pad+songlen

	left_edge=min(find(old_t{i}(:,1)>PADDING(1)));
	right_edge=max(find(old_t{i}(:,1)<(SONG_LEN-PADDING(2))));

	% check for gain shifts or LED off

	%  left_edge
	%  right_edge
	%  size(FORM_DATA{i})

	flag1=any(any(abs(diff(FORM_DATA{i}(left_edge:right_edge,:)))>THRESH));
	flag2=any(any(FORM_DATA{i}(left_edge:right_edge,:)<THRESH2));


	% add to trials for deletion if either flag trips

	if flag1|flag2
		to_del=[to_del i];
	end

end

FORM_DATA(to_del)=[];
FORM_DATE(to_del)=[];
FORM_MOTIF(to_del)=[];
FILENAMES(to_del)=[];
FORM_TRIALS(to_del)=[];
old_t(to_del)=[];

mintime=0;
maxtime=SONG_LEN;
newtime=[mintime:1/NEWFS:maxtime]';

% spline interpolate with extrapolation, remove all data outside of pads

% two methods for stitching the data:
% 1) interpolate within trials, average
% 2) bin data, interpolate over empty bins

if METHOD==1
	for i=1:length(FORM_DATA)
		%old_t{i}=old_t{i}-PADDING(1);
		FORM_DATA{i}=interp1(old_t{i}(:,1),FORM_DATA{i},newtime,'spline');
	end
else

	bins_x=newtime;
	bins_y=cell(1,length(bins));

	for i=1:length(bins_y)
		bins_y{i}=[];
	end

	% collect time points in bin

	for i=1:length(FORM_DATA)

		% map to closest time idx

		ca_t=old_t{i}(:,1);

		for j=1:length(ca_t)
			[~,loc]=min(bins_x-ca_t(j));
			bins_y{loc}(end+1)=FORM_DATA{i}(j,:);
		end

	end

	% average non-empty bins, rest get NaNs, interpolate through the NaNs

	bins_y
	pause();

end

FORM_T=newtime;

for i=1:length(FORM_DATA)
	FORM_DATA{i}=fluolab_detrend(FORM_DATA{i},'fs',NEWFS,'method','prctile','win',0,'per',12);
end

% cat to 3D matrix, samples x rois x trials

FORM_DATA=cat(3,FORM_DATA{:});

% finally hampel filter to remove any remaining outliers...

[nsamples,nrois,ntrials]=size(FORM_DATA);

to_del=[];

for i=1:nsamples
	tmp=squeeze(FORM_DATA(i,:,:));
	tmp=mean(tmp);
	mu=median(tmp);
	sig=mad(tmp,1);
	thresh=mu+sig_t*sig; % yes it's hand-tuned
	to_del=[to_del find(tmp>thresh)];
end

to_del=unique(to_del);
FORM_DATA(:,:,to_del)=[];
FORM_DATE(to_del)=[];
FORM_MOTIF(to_del)=[];
FILENAMES(to_del)=[];
FORM_TRIALS(to_del)=[];
