function COLLECT_DATA=stan_cadata_collect_freedomscope_v2(DIR)
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

MINT=[];
MAXT=[];
for i=1:length(listing)
	load(listing(i).name,'roi_ave');
	[COLLECT_DATA{i},TIME]=stan_cadata_format_freedomscope_v2(roi_ave.RAWdat,roi_ave.RawTime,50,30,100,MINT,MAXT);
	MINT=min(TIME);
	MAXT=max(TIME);
end
