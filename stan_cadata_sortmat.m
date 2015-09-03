function FORM_DATA=stan_format_cadata(DATA,varargin)
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
sort_day=3; % day to use for sorting
peak_check=0; % check for peak consistency 
peak_thresh=.05; % if closest peak is >peak_thresh, exclude roi
dff_check=.5; % check for dff peak 
chk_day=1; % check for dff peak day
scaling='r'; % scaling ('r' for within roi across days, 's' for within roi sort day, 'l' for within roi and day)
smoothing=0; % smooth ca trace (not working yet)
smooth_kernel='g'; % gauss smoothing kernel (b for boxcar)
fig_row=1; % subplot row 
fig_nrows=1; % total number of subplot rows
padding=1; % padding before and after song
bin_fluo=0; % discretize fluorescence 
nbins=10; % number of bins for discretization

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
		case 'fig_row'
			fig_row=varargin{i+1};
		case 'fig_nrows'
			fig_nrows=varargin{i+1};
		case 'chk_day'
			chk_day=varargin{i+1};
		case 'bin_fluo'
			bin_fluo=varargin{i+1};
		case 'nbins'
			nbins=varargin{i+1};
        case 'padding'
            padding=varargin{i+1};
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
	ave_mat{i}=mean(DATA{i},3);
end

% check for high enough dff

inc_rois=find(any(ave_mat{chk_day}>dff_check));
nrois=length(inc_rois);

for i=1:ndays
	ave_mat{i}=ave_mat{i}(:,inc_rois);
end

% renormalize within each day, using min/max from sort day or within each day, or
% within each roi across days
%

if strcmp(lower(scaling(1)),'l')

	disp('ROI within day');

	for i=1:ndays
		minca=repmat(min(ave_mat{i}),[nsamples 1]);
		maxca=repmat(max(ave_mat{i}),[nsamples 1]);
		ave_mat{i}=(ave_mat{i}-minca)./(maxca-minca);
	end

elseif strcmp(lower(scaling(1)),'r')

	disp('ROI across days');

	% cycle through each roi, agg ca data across days, normalize using 
	% across day min/max

	tmp=cat(1,ave_mat{:});
	minca=repmat(min(tmp),[nsamples 1]);
	maxca=repmat(max(tmp),[nsamples 1]);

	for i=1:ndays
		ave_mat{i}=(ave_mat{i}-minca)./(maxca-minca);
	end

elseif strcmp(lower(scaling(1)),'s')

	disp('ROI sort day');

	% use sort day min/max for normalization

	tmp=ave_mat{sort_day};

	minca=min(tmp);
	maxca=max(tmp);

	mincamat=repmat(minca,[nsamples 1]);
	maxcamat=repmat(maxca,[nsamples 1]);

	% requires clipping in case other days exceed limits

	for i=1:ndays
		for j=1:nrois

			% set <min to 0
			% set >max to max

			minidx=(ave_mat{i}(:,j)-minca(j))<0;
			maxidx=(ave_mat{i}(:,j))>maxca(j);

			ave_mat{i}(minidx,j)=0;
			ave_mat{i}(maxidx,j)=maxca(j);

		end

		ave_mat{i}=(ave_mat{i}-mincamat)./(maxcamat-mincamat);

	end

else

end

% get the sort indices

movie_x=[0:nsamples-1]/movie_fs;
pad_smps=padding*movie_fs;

if upsample>1

	for i=1:ndays
		interp_x=[0:1/upsample:nsamples-1]/movie_fs;
		ave_mat{i}=interp1(1:nsamples,ave_mat{i},[1:1/upsample:nsamples],'spline');
	end

	pad_smps=padding*(movie_fs*upsample);
	nsamples=size(ave_mat{1},1);
end

% get peak locations, anything in the pads is removed

if smoothing
	phase_shift=round(length(kernel)/2);
end

[~,peakloc]=max(ave_mat{sort_day});
del=(peakloc<pad_smps|peakloc>nsamples-(pad_smps-phase_shift));

for i=1:ndays
	ave_mat{i}(:,del)=[];
end

% get peak locations again

[~,peakloc]=max(ave_mat{sort_day});
[~,idx]=sort(peakloc)

% plot

if bin_fluo
	bin_edges=0:1/nbins:1
	for i=1:ndays
		[~,ave_mat{i}]=histc(ave_mat{i},bin_edges);
		ave_mat{i}
	end
end

for i=1:ndays
	ax(i)=subplot(fig_nrows,ndays,(fig_row-1)*ndays+i);
	imagesc(interp_x,[],ave_mat{i}(:,idx)');

	if bin_fluo
		caxis([0 nbins+1]);
	else
		caxis([0 1]);
	end
	set(gca,'xtick',[],'ytick',[]);
end

linkaxes(ax,'xy');
