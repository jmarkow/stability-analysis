function [BIRD,MOTIF]=stan_read_filename(FILENAME)
% read bird and motif from filename
%

[pathname,filename,ext]=fileparts(FILENAME);

tmp=regexp(filename,'^([a-zA-Z0-9]+)_.*','tokens');
BIRD=tmp{1}{1};
tmp=regexp(filename,'^[a-zA-Z0-9]+_(.*)$','tokens');
MOTIF=tmp{1}{1};
