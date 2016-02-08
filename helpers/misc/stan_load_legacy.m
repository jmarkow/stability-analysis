function [DATA,CHANNELS,FS]=stan_load_legacy(FILENAME)
%
%
%
%

vars=whos('-file',FILENAME);

var_names={vars(:).name};
islegacy=any(strcmp(var_names,'EPHYS_DATA'));

if islegacy
	load(FILENAME,'EPHYS_DATA','fs','CHANNELS');
	if ~exist('fs','var')
		fs=25e3;
	end
	DATA=EPHYS_DATA;
	FS=fs;
else
	load(FILENAME,'agg_ephys');
	DATA=agg_ephys.data;
	CHANNELS=agg_ephys.labels;
	FS=agg_ephys.fs;
end


