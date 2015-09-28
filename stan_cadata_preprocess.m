function [DATA,PHASE_SHIFT]=stan_format_cadata(DATA,varargin)
% takes data from stan_format_cadata and generates a series of panels for each time point
%
%
%
%
%
%

movie_fs=22; % sampling rate of camera
sort_day=1; % day to use for sorting
peak_check=0; % check for peak consistency 
peak_thresh=.05; % if closest peak is >peak_thresh, exclude roi
dff_check=.5; % check for dff peak 
scaling='r'; % scaling ('r' for within roi across days, 's' for within roi sort day, 'l' for within roi and day)
smoothing=0; % smooth ca trace (not working yet)
smooth_kernel='g'; % gauss smoothing kernel (b for boxcar)
chk_day=1;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'movie_fs'
			movie_fs=varargin{i+1};
		case 'sort_day'
			sort_day=varargin{i+1};
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
        case 'chk_day'
            chk_day=varargin{i+1};
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

if ~iscell(DATA)
	error('Wrong data format...');
end

ndays=length(DATA);
[nsamples,nrois,ntrials]=size(DATA{1});

if peak_check
	exclude=[];
	sortca.peaks=cell(1,ntrials);
	sortca.mergevals=cell(1,ntrials);

	% requires that the fb toolbox is in your PATH

	for i=1:ntrials
		[sortca.peaks{i},sortca.vals{i}]=fb_compute_peak_simple(DATA{chk_day}(:,:,i),...
			'thresh_t',.2,'debug',0,'onset_only',0,'thresh_hi',1,'thresh_int',8,'thresh_dist',.2); % thresh_int previously 5
	end

	% now flip through criteria to include
	% get timing distance relative to first trial, exclude inconsistent rois

	for i=1:nrois

		% concatenate times, exclude double peaks
	
		template=sortca.peaks{1}{i};

		if	isempty(template)
			exclude(end+1)=i;
			continue;
		end

		% otherwise, check other times
	
		for j=1:ntrials
			mindist=inf;
			for k=1:length(template)
				tmp=min(abs(template(k)-sortca.peaks{j}{i}))/movie_fs;
				if tmp<mindist
					mindist=tmp;
				end			
			end
		end

		if mindist>peak_thresh
			exclude(end+1)=i;
		end

	end

	exclude
	inc_rois=setdiff(1:nrois,exclude);

else

	inc_rois=1:nrois;

end

ave_mat=cell(1,ndays);

if smoothing>0

	ts=round(smoothing*movie_fs);
	
	if strcmp(smooth_kernel,'b')	
		kernel=ones(ts,1)/ts;
	elseif strcmp(smooth_kernel,'g')
		
		kernx=[-3*smoothing:1/movie_fs:3*smoothing];
		kernel=normpdf(kernx,0,smoothing);

		%kernel=kernel./sum(kernel);
	end

	for i=1:ndays

		% boxcar smoothing, swap out kernel here (gauss more appropriate?)
		% zero pad or repeat extreme values (or wrap)

		disp(['Note causal filter length ' num2str(length(kernel))]);
		
		tmp=DATA{i};

		% wrap around for smoothing
	
		zeropad_len=round(length(kernel)/2);
	
		%zeropad=repmat(DATA{i}(1,:,:),[zeropad_len 1 1]);	
		zeropad=DATA{i}(end-(zeropad_len-1):end,:,:);

		tmp=[zeropad;tmp];
		
		tmp=filter(kernel,1,tmp);
		DATA{i}=tmp(zeropad_len+1:end,:,:);

	end

end

for i=1:ndays
	DATA{i}=DATA{i}(:,inc_rois,:);
end

% get peak locations, anything in the pads is removed

if smoothing
	PHASE_SHIFT=round(length(kernel)/2);
else
	PHASE_SHIFT=0;
end


