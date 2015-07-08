function OPTIONS=stan_read_options(FILENAME,varargin)
% script for reading config files/sorting for processing/chopping
% 
% takes logfile as input 
%
%
%

fid=fopen(FILENAME,'r');
readdata=textscan(fid,'%s%[^\n]','commentstyle','#','delimiter','=');
fclose(fid);

OPTIONS=struct();

for i=1:length(readdata{1})
	OPTIONS.(readdata{1}{i})=readdata{2}{i};
end
