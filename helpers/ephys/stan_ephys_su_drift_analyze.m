function [pmat,zmat]=stan_ephys_su_drift_analyze(DATA,varargin)
% takes data from stan_format_cadata and generates a series of panels for each time point
%
%
%
%
%
%

compare_day=1; % day to use as basis for comparison
nperms=1e3; % as expected, permutation and ranksum give roughly the same answer
method='r'; % (r)anksum, (t)test, (p)ermutation (note that permutation is dog slow)
nparams=length(varargin);
tail='right';

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'compare_day'
			compare_day=varargin{i+1};
        case 'padding'
            padding=varargin{i+1};
		case 'tail'
			tail=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'nperms'
			nperms=varargin{i+1};
	end
end

% take sort day, clean up according to criteria (consistent peak? high ave?)

nbirds=length(DATA);

pmat={};
zmat={};

for i=1:nbirds

	ndays=length(DATA(i).smooth_rate);
	
	tmp=corr(DATA(i).smooth_rate{compare_day});
	corrvals=tmp(find(triu(ones(size(tmp)),1)));

	for j=1:ndays
	
		tmp=corr(DATA(i).smooth_rate{compare_day},DATA(i).smooth_rate{j});
		inc_vals=tmp(find(triu(ones(size(tmp)),1)));
	
		% use either t-test, ranksum or permutation (all yield similar answers)

		switch lower(method(1))
			case 't'
				[~,pmat{i}(j)]=ttest2(corrvals,inc_vals,'tail',tail);
			case 'r'
				pmat{i}(j)=ranksum(corrvals,inc_vals,'tail',tail);
			case 'p'
				nulldist=zeros(nperms,1);
				allvals=[corrvals;inc_vals(:)];

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

				obs=(mean(corrvals)-mean(inc_vals))/(sqrt(std(corrvals)*std(inc_vals)));

				% check appropriate tail

				if strcmp(lower(tail(1)),'r')
					tmp=sum(obs<nulldist); % right tail 
				elseif strcmp(lower(tail(1)),'l')
					tmp=sum(obs>nulldist); % left tail
				else
					tmp=sum(abs(obs)>abs(nulldist)); % two tail
				end

				pmat{i}(j)=(tmp+1)/(nperms+1);  
	
			otherwise
				error('Did not understand method selection');		
		end	
				
		zmat{i}(j)=(mean(inc_vals)-mean(corrvals))/(sqrt(std(corrvals)*std(inc_vals))); 

	end
end
