function FORM_DATA=stan_cadata_format_freedomscope(DATA,THRESH)
% takes data where cell arrays correspond to separate songs, rows to samples, and columns to rois
% and reformats for stan_cadata_sortmat
%
% e.g. newdata=stan_format_cadata(data1,data2,data3);
%
% returns a new cell array where each element is a sample x roi x trial matrix

if nargin<2 | isempty(THRESH)
	THRESH=50;
end

% get sample sizes, take max

trial_num=cellfun(@(x) size(x,2),DATA);
maxtrial=max(trial_num);
idx=trial_num==maxtrial;

DATA=DATA(idx);

cam_off=cellfun(@(x) any(any(abs(diff(x'))>THRESH)),DATA);
DATA(cam_off)=[];

for i=1:length(DATA)
	DATA{i}=fluolab_detrend(DATA{i}(:,1:end,:)','fs',30,'method','prctile','win',.4,'per',12);
end

FORM_DATA=cat(3,DATA{:});
