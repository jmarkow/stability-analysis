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
peak_check=0; % check for peak consistency 
peak_thresh=.05; % if closest peak is >peak_thresh, exclude roi
dff_check=1;
chk_day=1; % check for dff peak day
scaling='r'; % scaling ('r' for within roi across days, 's' for within roi sort day, 'l' for within roi and day)
smoothing=0; % smooth ca trace (not working yet)
smoothing_kernel='g'; % gauss smoothing kernel (b for boxcar)
padding=1; % padding before and after song
compare_day=1; % day to use as basis for comparison
nperms=1e5; % as expected, permutation and ranksum give roughly the same answer
method='r'; % (r)anksum, (t)test, (p)ermutation (note that permutation is dog slow)
nparams=length(varargin);
tail='right';

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
		case 'chk_day'
			chk_day=varargin{i+1};
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
		case 'tail'
			tail=varargin{i+1};
		case 'method'
			method=varargin{i+1};
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

if ~iscell(DATA)
	error('Wrong data format...');
end

ndays=length(DATA);

[DATA,phase_shift]=stan_cadata_preprocess(DATA,'peak_check',peak_check,'peak_thresh',peak_thresh,'movie_fs',movie_fs,...
	'smoothing',smoothing,'smoothing_kernel',smoothing_kernel);


% get the sort indices

for i=1:ndays
	ave_mat{i}=mean(DATA{i},3);
end

% check for high enough dff

inc_rois=find(any(ave_mat{chk_day}>dff_check))

[nsamples,nrois,ntrials]=size(DATA{1});

for i=1:ndays
	ave_mat{i}=ave_mat{i}(:,inc_rois);
	DATA{i}=DATA{i}(:,inc_rois,:);
end

movie_x=[0:nsamples-1]/movie_fs;
pad_smps=padding*movie_fs;

% any peaks outside of the pads?

[~,peakloc]=max(ave_mat{chk_day});
del=(peakloc<pad_smps|peakloc>nsamples-(pad_smps-phase_shift))

for i=1:ndays
	DATA{i}(:,del,:)=[];
end

[nsamples,nrois,ntrials]=size(DATA{1});

% not necessary here, but left in just in case

if upsample>1

	for i=1:ndays
		interp_x=[0:1/upsample:nsamples-1]/movie_fs;
		DATA{i}=interp1(1:nsamples,DATA{i},[1:1/upsample:nsamples],'spline');
	end

	pad_smps=padding*(movie_fs*upsample);
	nsamples=size(DATA{1},1);
	phase_shift=round((phase_shift/movie_fs)*(movie_fs*upsample));

end

% trim pads and zscore

for i=1:ndays
	DATA{i}=zscore(DATA{i});
	DATA{i}=DATA{i}(pad_smps:end-pad_smps,:,:);
end

% get corr values from the comparison day

ntrials=size(DATA{compare_day},3);
corrvals=zeros(nchoosek(ntrials,2),nrois);

% grab values from upper triangle of corr matrix

for i=1:nrois
	tmp=corr(squeeze(DATA{compare_day}(:,i,:)));
	corrvals(:,i)=tmp(find(triu(ones(size(tmp)),1)));
end

comparevals=cell(ndays,nrois);
pmat=zeros(ndays,nrois);
zmat=zeros(ndays,nrois);

for i=1:ndays
	for j=1:nrois
	
		tmp=corr(squeeze(DATA{compare_day}(:,j,:)),squeeze(DATA{i}(:,j,:)));
		inc_vals=tmp(find(triu(ones(size(tmp)),1)));
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
