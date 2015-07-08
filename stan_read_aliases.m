function [TARGET,SOURCE,DATE_TARGET,DATE_SOURCE]=stan_read_aliases(FILENAME)
%
%
%

SOURCE={};
TARGET={};
DATE_TARGET={};
DATE_SOURCE=[];

% read in alias file

fid=fopen(FILENAME,'r');
tmp=textscan(fid,'%s','commentstyle','#','delimiter','=');
fclose(fid);

if length(tmp{1})>=2
	TARGET=tmp{1}(1:2:end);
	SOURCE=tmp{1}(2:2:end);
end

isdate=zeros(1,length(SOURCE));

for i=1:length(SOURCE)
	isdate(i)=~isempty(regexp(SOURCE{i},'\d+\/\d+\/\d+','tokens'));
end

isdate=find(isdate);
tmp=SOURCE(isdate);
DATE_TARGET=TARGET(isdate);

DATE_SOURCE=zeros(1,length(tmp));

for i=1:length(tmp)
	DATE_SOURCE(i)=datenum(tmp{i});
end

TARGET(isdate)=[];
SOURCE(isdate)=[];


