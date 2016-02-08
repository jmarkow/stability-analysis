function FORM_DATA=stan_cadata_format(varargin)
% takes data where cell arrays correspond to separate songs, rows to samples, and columns to rois
% and reformats for stan_cadata_sortmat
%
% e.g. newdata=stan_format_cadata(data1,data2,data3);
%
% returns a new cell array where each element is a sample x roi x trial matrix

FORM_DATA=cell(1,length(varargin));

for i=1:length(varargin)
	FORM_DATA{i}=cat(3,varargin{i}{:});
end
