function [COLLECT_DATA,COLLECT_DATES,TIME]=stan_cadata_collect_freedomscope_v2(RAW)
% takes data where cell arrays correspond to separate songs, rows to samples, and columns to rois
% and reformats for stan_cadata_sortmat
%
% e.g. newdata=stan_format_cadata(data1,data2,data3);
%
% returns a new cell array where each element is a sample x roi x trial matrix

COLLECT_DATA=cell(1,length(RAW));
COLLECT_DATES=cell(1,length(RAW));
MINT=[];
MAXT=[];
SONG_LEN=1.8; % .59 for lny13 .625 for lny18 (.8 for each then pad later?)

% remove pad first here?

for i=1:length(RAW)

	if ~isfield(RAW{i},'Offset')
		RAW{i}.Offset=[];
	else
		len=cellfun(@length,RAW{i}.Offset);
		RAW{i}.Offset(len==0)=[];
	end

	if isfield(RAW{i},'filename')
		len=cellfun(@length,RAW{i}.filename);
		RAW{i}.filename(len==0)=[];
	else
		RAW{i}.filename=[];
	end

	if ~isfield(RAW{i},'padding')
		RAW{i}.padding=[.25 .75];
	end

	[COLLECT_DATA{i},TIME{i},COLLECT_DATES{i}]=stan_cadata_format_freedomscope_v2(RAW{i}.raw_dat,...
		RAW{i}.raw_time,...
		60,... % threshold on derivative (check for gain shift, in raw px values)
		30,... % threshold for camera on (px values below this considered LED off)
		100,... % new sampling rate
		MINT,... % minimum time point for new time frame
		MAXT,... % maximum time point for new tie frame
		RAW{i}.padding,... % padding for extraction
		SONG_LEN,... % length of song
		RAW{i}.Offset,... % offset
    RAW{i}.filename); % filename (for parsing trial times)

	MINT=min(TIME{i});
	MAXT=max(TIME{i});

end
