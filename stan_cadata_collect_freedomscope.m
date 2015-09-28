function COLLECT_DATA=stan_format_cadata(DIR)
% takes data where cell arrays correspond to separate songs, rows to samples, and columns to rois
% and reformats for stan_cadata_sortmat
%
% e.g. newdata=stan_format_cadata(data1,data2,data3);
%
% returns a new cell array where each element is a sample x roi x trial matrix

if nargin<1 | isempty(DIR)
	DIR=pwd;
end

listing=dir(fullfile(DIR,'*.mat'));
COLLECT_DATA=cell(1,length(listing));

for i=1:length(listing)
	load(listing(i).name,'roi_ave');
	COLLECT_DATA{i}=stan_cadata_format_freedomscope(roi_ave.raw);
end

maxsamples=inf;

for i=1:length(COLLECT_DATA)
	nsamples=size(COLLECT_DATA{i},1);

	if nsamples<maxsamples
		maxsamples=nsamples;
	end
end

for i=1:length(COLLECT_DATA)
	COLLECT_DATA{i}=COLLECT_DATA{i}(1:maxsamples,:,:);
end

