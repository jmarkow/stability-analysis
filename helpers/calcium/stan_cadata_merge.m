function COMBINED_DATA=stan_cadata_merge(varargin)
% merges rois from multiple files
%
%
%
%
%

if length(varargin)<2
	error('Need more than one file to continue');
end;

nfiles=length(varargin);

for i=1:nfiles
	merge.(['set' num2str(i)])=load(varargin{i});
end


fnames=fieldnames(merge.set1);

tmp=regexp(fnames,'data\d+','match');
hits={};

for i=1:length(tmp)
	if ~isempty(tmp{i}{1})
		hits{end+1}=tmp{i}{1};
	end
end

nhits=length(hits);

for i=1:nhits
	
	ntrials=length(merge.set1.(hits{i}));
	COMBINED_DATA.(hits{i})=cell(1,ntrials);

	for j=1:nfiles
		for k=1:ntrials
			COMBINED_DATA.(hits{i}){k}=[ COMBINED_DATA.(hits{i}){k} merge.(['set' num2str(j)]).(hits{i}){k} ];
		end
	end
end
