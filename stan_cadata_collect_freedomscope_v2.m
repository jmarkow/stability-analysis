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
SONG_LEN=.8; % .59 for lny13 .625 for lny18 (.8 for each then pad later?)

% remove pad first here? 

for i=1:length(listing)

	load(fullfile(DIR,listing(i).name),'roi_ave');

	% if not offset data, set to empty

	disp([listing(i).name]);

	if ~isfield(roi_ave,'Offset')
		roi_ave.Offset=[];
	else
		len=cellfun(@length,roi_ave.Offset);
		roi_ave.Offset(len==0)=[];
	end

	%roi_ave.Offset=[];

	[COLLECT_DATA{i},TIME]=stan_cadata_format_freedomscope_v2(roi_ave.RAWdat,...
		roi_ave.RawTime,...
		60,... % threshold on derivative (check for gain shift, in raw px values)
		30,... % threshold for camera on (px values below this considered LED off)
		100,... % new sampling rate
		MINT,... % minimum time point for new time frame
		MAXT,... % maximum time point for new tie frame
		roi_ave.padding,... % padding for extraction 
		SONG_LEN,... % length of song
		roi_ave.Offset); % offset

	MINT=min(TIME);
	MAXT=max(TIME);

end
