function [corrvals,comparevals,pmat,zmat]=stan_format_cadata(DATA,varargin)
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
sort_day=1; % day to use for sorting
peak_check=0; % check for peak consistency 
peak_thresh=.05; % if closest peak is >peak_thresh, exclude roi
chk_day=1; % check for dff peak day
scaling='r'; % scaling ('r' for within roi across days, 's' for within roi sort day, 'l' for within roi and day)
smoothing=0; % smooth ca trace (not working yet)
smoothing_kernel='g'; % gauss smoothing kernel (b for boxcar)
padding=1; % padding before and after song
compare_day=1;
nperms=1e5;

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
		case 'smoothing_kernel'
			smoothing_kernel=varargin{i+1};
		case 'compare_day'
			compare_day=varargin{i+1};
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

[DATA,phase_shift]=stan_cadata_preprocess(DATA,'peak_check',peak_check,'peak_thresh',peak_thresh,'movie_fs',movie_fs,...
	'smoothing',smoothing,'smoothing_kernel',smoothing_kernel);
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
	phase_shift=round((phase_shift/movie_fs)*(movie_fs*upsample))
end

% get peak locations, anything in the pads is removed

pad_smps

for i=1:ndays
	DATA{i}=zscore(DATA{i});
	DATA{i}=DATA{i}(pad_smps:end-pad_smps,:,:);
end

% for each day build corr vals

corrvals=zeros(nchoosek(ntrials,2),nrois);

for i=1:nrois
	tmp=corr(squeeze(DATA{1}(:,i,:)));
	corrvals(:,i)=tmp(find(triu(ones(size(tmp)),1)));
end

comparevals=cell(ndays-1,nrois);
pmat=zeros(ndays,nrois);
zmat=zeros(ndays,nrois);

for i=1:ndays
	for j=1:nrois
	
		j

		tmp=corr(squeeze(DATA{1}(:,j,:)),squeeze(DATA{i}(:,j,:)));
		inc_vals=tmp(find(triu(ones(size(tmp)),1)));
		comparevals{i,j}=inc_vals;
		%[~,pmat(i,j)]=ttest2(corrvals(:,j),inc_vals,'tail','right');
		pmat(i,j)=ranksum(corrvals(:,j),inc_vals,'tail','right');
		nulldist=zeros(nperms,1);
		allvals=[corrvals(:,j);inc_vals(:)];

		idx1=1:size(corrvals,1);
		idx2=size(corrvals,1)+1:size(allvals);

		[~,rndidx]=sort(rand(nperms,length(allvals)),2);
		rndvals=allvals(rndidx);
		pop1=rndvals(:,idx1)';
		pop2=rndvals(:,idx2)';

		nulldist=(mean(pop1)-mean(pop2))./(sqrt(std(pop1).*std(pop2)));

		%for k=1:nperms
		%	pop1=rndvals(k,idx1);
		%	pop2=rndvals(k,idx2);
		%	nulldist(k)=(mean(pop1)-mean(pop2));
		%end

		obs=(mean(corrvals(:,j))-mean(inc_vals))/(sqrt(std(corrvals(:,j))*std(inc_vals)));
		tmp=sum(obs<nulldist); % permutation p-val
		pmat(i,j)=(tmp+1)/(nperms+1);

		% permutation test
		
		zmat(i,j)=mean((inc_vals-mean(corrvals(:,j)))/std(corrvals(:,j)));

	end
end

%
%[~,peakloc]=max(ave_mat{sort_day});
%del=(peakloc<pad_smps|peakloc>nsamples-(pad_smps-phase_shift));
%
%for i=1:ndays
%	ave_mat{i}(:,del)=[];
%end
%
%% get peak locations again
%
%[~,peakloc]=max(ave_mat{sort_day});
%[~,idx]=sort(peakloc)
