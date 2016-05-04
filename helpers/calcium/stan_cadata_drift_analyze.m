function [rmat_mu,pmat]=stan_cadata_drift_analyze(DATA,IDX,varargin)
% takes data from stan_format_cadata and generates a series of panels for each time point
%
%
%
%
%
%

if nargin<2 | isempty(IDX)
	IDX=1:length(DATA);
end

pmat=[];

movie_fs=22; % sampling rate of camera
peak_check_pad=0; % check for peak consistency
peak_thresh=.05; % if closest peak is >peak_thresh, exclude roi
dff_check=1;
chk_day=1; % check for dff peak day
smoothing=0; % smooth ca trace (not working yet)
smooth_kernel='g'; % gauss smoothing kernel (b for boxcar)
padding=1; % padding before and after song
compare_day=1; % day to use as basis for comparison
nparams=length(varargin);
tail='right';
maxlag=.02;
lag_corr=0;
realign=1;
nboots=1e4;
frac=2;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'movie_fs'
			movie_fs=varargin{i+1};
		case 'peak_check_pad'
			peak_check_pad=varargin{i+1};
		case 'peak_thresh'
			peak_thresh=varargin{i+1};
		case 'dff_check'
			dff_check=varargin{i+1};
		case 'chk_day'
			chk_day=varargin{i+1};
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
		case 'lag_corr'
			lag_corr=varargin{i+1};
		case 'realign'
			realign=varargin{i+1};
		case 'maxlag'
			maxlag=varargin{i+1};
  	case 'nboots'
      nboots=varargin{i+1};
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

maxlag_smps=round(movie_fs*maxlag);

if ~iscell(DATA)
	error('Wrong data format...');
end

ndays=length(DATA);

[DATA,phase_shift]=stan_cadata_preprocess(DATA,'peak_check_pad',peak_check_pad,'peak_thresh',peak_thresh,'movie_fs',movie_fs,...
	'smoothing',smoothing,'smooth_kernel',smooth_kernel,'padding',padding,'realign',realign,'maxlag',maxlag);

% get the sort indices

for i=1:ndays
	ave_mat{i}=mean(DATA{i},3);
end

[nsamples,nrois,ntrials]=size(DATA{1});
pad_smps=round(padding*movie_fs);

if pad_smps(1)==0;
	pad_smps(1)=1;
end

[~,peakloc]=max(ave_mat{compare_day});
del=(peakloc<pad_smps(1)|peakloc>nsamples-(pad_smps(2)));

%for i=1:ndays
%	DATA{i}(:,del,:)=[];
%end

%inc_rois=find(any(ave_mat{chk_day}(pad_smps(1):end-pad_smps(2),:)>dff_check))
%pad_smps=round(padding*movie_fs);

% trim pads and zscore

for i=1:ndays
	DATA{i}=zscore(DATA{i}(pad_smps(1):end-pad_smps(2),:,:));
end

% get corr values from the comparison day

ntrials=size(DATA{compare_day},3);
pairs=nchoosek(1:ntrials,2);
npairs=size(pairs,1);
corrvals=zeros(npairs,nrois);

% get the correlation difference using means (probably the most reliable here,
% jitter should average out assuming it's symmetric)

mu1=mean(zscore(DATA{compare_day}),3);

for i=1:ndays

	mu2=mean(zscore(DATA{i}),3);
	corrmat=zeros(nrois,nrois);

	if lag_corr
		for j=1:nrois
			corrmat(j,j)=max(xcorr(mu1(:,j),mu2(:,j),maxlag_smps,'coeff'));
		end
	else
		corrmat=corr(mu1,mu2,'type','pearson');
	end

	rmat_mu.acrossday(i,:)=corrmat(find(diag(ones(nrois,1),0)));

end

%rmat_mu.withinday=nan(ndays,nrois);


ntrials=size(DATA{1},3);

pool1=1:floor(ntrials/frac);
pool2=floor(ntrials/frac)+1:ntrials;
diag_idx=find(diag(ones(nrois,1),0));
maxlag=(max(IDX)-min(IDX))+1;

rmat_mu.lag.day=cell(1,maxlag);
rmat_mu.lag.night=cell(1,maxlag);
rmat_mu.lag.all=cell(1,maxlag);

rmat_mu.bootstrap.lag.day=cell(1,maxlag);
rmat_mu.bootstrap.lag.night=cell(1,maxlag);
rmat_mu.bootstrap.lag.all=cell(1,maxlag);

for i=1:ndays

	ntrials1=size(DATA{i},3);
	pool_compare=(ntrials1-(round(ntrials1/frac)-1)):ntrials1;

	compare=mean(zscore(DATA{i}(:,:,pool_compare)),3);
	compare_all=mean(zscore(DATA{i}),3);

	for j=i:ndays

		ntrials2=size(DATA{j},3);

		pool1=1:round(ntrials2/frac);
		pool2=(ntrials2-(round(ntrials2/frac)-1)):ntrials2;

		mu_day=mean(zscore(DATA{j}(:,:,pool1)),3);
		mu_night=mean(zscore(DATA{j}(:,:,pool2)),3);
		mu_all=mean(zscore(DATA{j}),3);

		corrmat_day=zeros(nrois,nrois);
		corrmat_night=zeros(nrois,nrois);
		corrmat_all=zeros(nrois,nrois);

		if lag_corr
			for k=1:nrois
				corrmat_day(k,k)=max(xcorr(compare(:,k),mu_day(:,k),maxlag_smps,'coeff'));
				corrmat_night(k,k)=max(xcorr(compare(:,k),mu_night(:,k),maxlag_smps,'coeff'));
				corrmat_all(k,k)=max(xcorr(compare_all(:,k),mu_all(:,k),maxlag_smps,'coeff'));
			end
		else
			corrmat_day=corr(compare,mu_day,'type','pearson');
			corrmat_night=corr(compare,mu_night,'type','pearson');
			corrmat_all=corr(compare_all,mu_all,'type','pearson');
		end

		% assign the correct lag here, not simply (j-i)+1

		lag=(IDX(j)-IDX(i))+1;
		idx=size(rmat_mu.lag.day{lag},1);

		rmat_mu.lag.day{lag}(end+1,:)=corrmat_day(diag_idx);
		rmat_mu.lag.night{lag}(end+1,:)=corrmat_night(diag_idx);
		rmat_mu.lag.all{lag}(end+1,:)=corrmat_all(diag_idx);

		% monte carlo permutation test for each comparison

		rmat_mu.bootstrap.lag.day{lag}{end+1}=zeros(nboots,nrois);
		rmat_mu.bootstrap.lag.night{lag}{end+1}=zeros(nboots,nrois);
		rmat_mu.bootstrap.lag.all{lag}{end+1}=zeros(nboots,nrois);

		if nboots>0

			% monte carlo for day,night and all

			ntrials1=size(DATA{i},3);
			ntrials2=size(DATA{j},3);

			compare1=ntrials1-(floor(ntrials1/frac)-1):ntrials1;
			pool1=1:floor(ntrials2/frac);
			pool2=ntrials2-(floor(ntrials2/frac)-1):ntrials2;

			all_stitch=zscore(cat(3,DATA{i},DATA{j}));
			day_stitch=zscore(cat(3,DATA{i}(:,:,compare1),DATA{j}(:,:,pool1)));
			night_stitch=zscore(cat(3,DATA{i}(:,:,compare1),DATA{j}(:,:,pool2)));

			% get random permutations

			all_ntrials=size(all_stitch,3);
			day_ntrials=size(day_stitch,3);
			night_ntrials=size(night_stitch,3);

			[~,rndidx_all]=sort(rand(all_ntrials,nboots));
			[~,rndidx_day]=sort(rand(day_ntrials,nboots));
			[~,rndidx_night]=sort(rand(night_ntrials,nboots));

			all_split=floor(all_ntrials/2);
			day_split=floor(day_ntrials/2);
			night_split=floor(night_ntrials/2);

			tic
			for k=1:nboots

				mu_all1=mean(all_stitch(:,:,rndidx_all(1:all_split,k)),3);
				mu_all2=mean(all_stitch(:,:,rndidx_all(all_split+1:end,k)),3);

				mu_day1=mean(day_stitch(:,:,rndidx_day(1:day_split,k)),3);
				mu_day2=mean(day_stitch(:,:,rndidx_day(day_split+1:end,k)),3);

				mu_night1=mean(night_stitch(:,:,rndidx_night(1:night_split,k)),3);
				mu_night2=mean(night_stitch(:,:,rndidx_night(night_split+1:end,k)),3);

				corrmat_day=zeros(nrois,nrois);
				corrmat_night=zeros(nrois,nrois);
				corrmat_all=zeros(nrois,nrois);

				if lag_corr
					for l=1:nrois
						corrmat_day(l,l)=max(xcorr(mu_day1(:,l),mu_day2(:,l),maxlag_smps,'coeff'));
						corrmat_night(l,l)=max(xcorr(mu_night1(:,l),mu_night2(:,l),maxlag_smps,'coeff'));
						corrmat_all(l,l)=max(xcorr(mu_all1(:,l),mu_all2(:,l),maxlag_smps,'coeff'));
					end
				else
					corrmat_day=corr(mu_day1,mu_day2,'type','pearson');
					corrmat_night=corr(mu_night1,mu_night2,'type','pearson');
					corrmat_all=corr(mu_all1,mu_all2,'type','pearson');
				end

				rmat_mu.bootstrap.lag.day{lag}{end}(k,:)=corrmat_day(diag_idx);
				rmat_mu.bootstrap.lag.night{lag}{end}(k,:)=corrmat_night(diag_idx);
				rmat_mu.bootstrap.lag.all{lag}{end}(k,:)=corrmat_all(diag_idx);

			end
			toc
		end

	end

end
