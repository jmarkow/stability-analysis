function [STARTIDX,STOPIDX,SIG_FILTER]=stan_get_offset(SIG1,SIG2,varargin)
%
%
%
%


% use the smaller signal as a filter

len1=length(SIG1);
len2=length(SIG2);

