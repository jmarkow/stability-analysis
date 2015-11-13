function [corrvals,comparevals,pmat,zmat,rmat_mu,rmat_mu_withinday,rmat_mu_nightday]=stan_cadata_drift_analyze(DATA,varargin)
% takes data from stan_format_cadata and generates a series of panels for each time point
%
%
%
%
%
%

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
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

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

% not necessary here, but left in just in case

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

% trim pads and zscore

for i=1:ndays
	DATA{i}=zscore(DATA{i}(pad_smps(1):end-pad_smps(2),:,:));
end

% get corr values from the comparison day

ntrials=size(DATA{compare_day},3);
pairs=nchoosek(1:ntrials,2);
npairs=size(pairs,1);
corrvals=zeros(npairs,nrois);
maxlag_smps=round(maxlag*(movie_fs*upsample));

% grab values from upper triangle of corr matrix

for i=1:nrois

	if lag_corr
		tmp=zeros(ntrials,ntrials);
		[x,y]=find(triu(ones(size(tmp)),1));

		% fill upper triangle only

		for j=1:length(x)
			tmp(x(j),y(j))=max(xcorr(DATA{compare_day}(:,i,x(j)),DATA{compare_day}(:,i,y(j)),maxlag_smps,'coeff'));
		end

	else

		tmp=corr(squeeze(DATA{compare_day}(:,i,:)));

	end

	corrvals(:,i)=tmp(find(triu(ones(size(tmp)),1)));

end

comparevals=cell(ndays,nrois);
pmat=zeros(ndays,nrois);
zmat=zeros(ndays,nrois);
zmat_mu=zeros(ndays,nrois);

for i=1:ndays

	for j=1:nrois

		ntrials2=size(DATA{i},3);

		if lag_corr

			tmp=zeros(ntrials,ntrials2);

			[x,y]=find(triu(ones(size(tmp)),0));
			for k=1:length(x)
				tmp(x(k),y(k))=max(xcorr(DATA{compare_day}(:,j,x(k)),DATA{i}(:,j,y(k)),maxlag_smps,'coeff'));
			end

		else
			tmp=corr(squeeze(DATA{compare_day}(:,j,:)),squeeze(DATA{i}(:,j,:)));
		end

		inc_vals=tmp(find(triu(ones(size(tmp)),0)));
		comparevals{i,j}=inc_vals;

		% use either t-test, ranksum or permutation (all yield similar answers)

		switch lower(method(1))
			case 't'
				[~,pmat(i,j)]=ttest2(corrvals(:,j),inc_vals,'tail',tail);
			case 'r'
				pmat(i,j)=ranksum(corrvals(:,j),inc_vals,'tail',tail);
			case 'p'
				nulldist=zeros(nperms,1);
				allvals=[corrvals(:,j);inc_vals(:)];

				idx1=1:size(corrvals,1);
				idx2=size(corrvals,1)+1:size(allvals);

				% random permutation nperms

				[~,rndidx]=sort(rand(nperms,length(allvals)),2);

				% permuted data

				rndvals=allvals(rndidx);

				% data splits

				pop1=rndvals(:,idx1)';
				pop2=rndvals(:,idx2)';

				% test statistic

				nulldist=(mean(pop1)-mean(pop2))./(sqrt(std(pop1).*std(pop2)));
				obs=(mean(corrvals(:,j))-mean(inc_vals))/(sqrt(std(corrvals(:,j))*std(inc_vals)));

				% check appropriate tail

				if strcmp(lower(tail(1)),'r')
					tmp=sum(obs<nulldist); % right tail
				elseif strcmp(lower(tail(1)),'l')
					tmp=sum(obs>nulldist); % left tail
				else
					tmp=sum(abs(obs)>abs(nulldist)); % two tail
				end

				pmat(i,j)=(tmp+1)/(nperms+1);

			otherwise
				error('Did not understand method selection');
			end

			zmat(i,j)=(mean(inc_vals)-mean(corrvals(:,j)))/(sqrt(std(corrvals(:,j))*std(inc_vals)));

		end
	end

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

	rmat_mu(i,:)=corrmat(find(diag(ones(nrois,1),0)));
end


rmat_mu_withinday=zeros(ndays,nrois);
rmat_mu_nightday=zeros(ndays-3,nrois);

for i=1:ndays

	corrmat=zeros(nrois,nrois);
	corrmat_nightday=zeros(nrois,nrois);

	ntrials=size(DATA{i},3);

	pool1=1:floor(ntrials/2);
	pool2=ntrials-(floor(ntrials/2)-1):ntrials;

	mu1=mean(zscore(DATA{i}(:,:,pool1)),3);
	mu2=mean(zscore(DATA{i}(:,:,pool2)),3);

	if lag_corr
		for j=1:nrois
			corrmat(j,j)=max(xcorr(mu1(:,j),mu2(:,j),maxlag_smps,'coeff'));
		end
	else
		corrmat=corr(mu1,mu2,'type','pearson');
	end

	rmat_mu_withinday(i,:)=corrmat(find(diag(ones(nrois,1),0)));

	% modify to compute with separate lags...

	if i+2<(ndays)

		ntrials2=size(DATA{i+3},3);
		pool3=(ntrials2-floor(ntrials2/2)):ntrials2
		pool3=1:floor(ntrials2/5);

		%pool3=1:ntrials2;

		mu3=mean(zscore(DATA{i+3}(:,:,pool3)),3);

		corrmat=zeros(nrois,nrois);
		compare=mu2;

		if lag_corr
			for j=1:nrois
				corrmat(j,j)=max(xcorr(compare(:,j),mu3(:,j),maxlag_smps,'coeff'));
			end
		else
			corrmat=corr(compare,mu3,'type','pearson');
		end

		rmat_mu_nightday(i,:)=corrmat(find(diag(ones(nrois,1),0)));

	end

end
