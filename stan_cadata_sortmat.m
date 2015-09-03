function FORM_DATA=stan_format_cadata(DATA,varargin)
%
%
%
%
%
%
%

movie_fs=22;
upsample=10;
thresh=2;
sort_day=3;
peak_check=0;
dff_check=.5;
scaling='r';
smoothing=0;
fig_row=1;
fig_nrows=1;
padding=.8;
chk_day=1;

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
		case 'thresh'
			thresh=varargin{i+1};
		case 'sort_day'
			sort_day=varargin{i+1};
		case 'peak_check'
			peak_check=varargin{i+1};
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
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

if ~iscell(DATA)
	error('Wrong data format...');
end

ndays=length(DATA);
[nsamples,nrois,ntrials]=size(DATA{1});

if peak_check
	sortca.peaks=cell(1,ntrials);
	sortca.mergevals=cell(1,ntrials);

	for i=1:ntrials
		[sortca.peaks{i},sortca.vals{i}]=fb_compute_peak_simple(DATA{1}(:,:,i),...
			'thresh_t',.2,'debug',0,'onset_only',0,'thresh_hi',1,'thresh_int',8,'thresh_dist',.2); % thresh_int previously 5
	end

	% now flip through criteria to include
	%
	% inc_rois= etc. etc.
	%
else

	inc_rois=1:nrois;

end

ave_mat=cell(1,ndays);

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

	tmp=ave_mat{sort_day};

	minca=min(tmp);
	maxca=max(tmp);

	mincamat=repmat(minca,[nsamples 1]);
	maxcamat=repmat(maxca,[nsamples 1]);

	for i=1:ndays
		for j=1:nrois

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
		ave_mat{i}=interp1(1:nsamples,ave_mat{i},[1:1/upsample:nsamples]);
	end
	pad_smps=padding*(movie_fs*upsample);
	nsamples=size(ave_mat{1},1);
end

[~,peakloc]=max(ave_mat{sort_day});
del=(peakloc<pad_smps|peakloc>nsamples-pad_smps);

for i=1:ndays
	ave_mat{i}(:,del)=[];
end

[~,peakloc]=max(ave_mat{sort_day});
[~,idx]=sort(peakloc)

for i=1:ndays
	ax(i)=subplot(fig_nrows,ndays,(fig_row-1)*ndays+i);
	imagesc(interp_x,[],ave_mat{i}(:,idx)');
	caxis([0 1]);
	set(gca,'xtick',[],'ytick',[]);
end

linkaxes(ax,'xy');
