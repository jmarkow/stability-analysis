function [LOG_MAP MAP]=stan_read_config(FILENAME,varargin)
% script for reading config files/sorting for processing/chopping
% 
% takes logfile as input 
%
%
%

fid=fopen(FILENAME,'r');
readdata=textscan(fid,'%s%[^\n]','commentstyle','#');
fclose(fid);

% ID to number map

bird_idx=strcmpi(readdata{1},'bird');

if isempty(bird_idx)
	error('No bird names found, bailing...');
end

bird_str=readdata{2}{bird_idx};

% parse bird_string, oh boy regex 

[bird_tokens]=regexpi(bird_str,'((?:\w*\s*|[0-9\/]*)+)(?:\[|\(|\s)((?:bird|\s*)\s*(?:\d|\#\d))(?:\]|\)|\s*)','tokens');

for i=1:length(bird_tokens)
	bird_tokens{i}{1}=regexprep(bird_tokens{i}{1},'(b|B)ird','');
	bird_tokens{i}=strtrim(bird_tokens{i});
end

for i=1:length(bird_tokens)
	LOG_MAP(i).name=lower(bird_tokens{i}{1});
	tmp=regexp(bird_tokens{i}{2},'.*(\d+)','tokens');

	if isempty(tmp)
		LOG_MAP(i).idx=1;
	else
		tmp=str2num(tmp{1}{1});
		LOG_MAP(i).idx=tmp;
	end

	% check for date number
	
	tmp=regexp(bird_tokens{i}{1},'(\d+\/\d+\/\d+)','tokens');

	if length(tmp)>0 
		LOG_MAP(i).date_num=datenum(tmp{1}{1});
	else
		LOG_MAP(i).date_num=[];
	end

	tmp=regexp(LOG_MAP(i).name,'^(\d+\/\d+) ','tokens');

	[pathname,~,~,]=fileparts(FILENAME);

	path_tokens=regexp(pathname,filesep,'split');

	if length(tmp)>0
		LOG_MAP(i).date_num=datenum([ tmp{1}{1} '/' path_tokens{end-3 }]);
		LOG_MAP(i).name=strtrim(regexprep(LOG_MAP(i).name,'(\d+\/)+\d+','')); % remove date from name
	end

	tmp=regexp(LOG_MAP(i).name,' (\d+\/\d+)$','tokens');

	if length(tmp)>0
		LOG_MAP(i).date_num=datenum([ tmp{1}{1} '/' path_tokens{end-3 }]);
		LOG_MAP(i).name=strtrim(regexprep(LOG_MAP(i).name,'(\d+\/)+\d+','')); % remove date from name
	end

end

% get the channel mapping first word is either column or nidaq

ch_idx=strcmpi(readdata{1},'column');

if ~any(ch_idx)
	ch_idx=strcmpi(readdata{1},'nidaq');
end


ch_tokens=readdata{2}(ch_idx);
cat_idx=cat(1,LOG_MAP(:).idx);

map_idx=ones(1,length(ch_tokens));
ch_to_bird=ones(1,length(ch_tokens));
ismic=zeros(1,length(ch_tokens));
ch_names=cell(1,length(ch_tokens));

if length(LOG_MAP)>1
	ch_flag=1;
else
	ch_flag=0;
end

for i=1:length(ch_tokens)

	tokens=regexp(ch_tokens{i},'>','split');
	map_token=tokens{end};
	tmp=regexp(map_token,'.+(\d+).+$','tokens');
	ismic(i)=~isempty(regexpi(map_token,'mic','match'));

	if ch_flag
		map_idx(i)=str2num(tmp{1}{1});
		ch_to_bird(i)=find(map_idx(i)==cat_idx);
	end

	ch_names{i}=map_token;

end

MAP.names=ch_names;

for i=1:length(LOG_MAP)
	if sum(ch_to_bird==i)>0
		LOG_MAP(i).ch.idx=find(ch_to_bird==i);
		LOG_MAP(i).ch.name=ch_names(LOG_MAP(i).ch.idx);
		LOG_MAP(i).ch.ismic=ismic(LOG_MAP(i).ch.idx)==1; % typecast to log
	end
end

% reconstruct date token from the folder structure, then we're finished



