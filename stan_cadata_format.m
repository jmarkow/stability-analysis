function FORM_DATA=stan_format_cadata(varargin)
%
%
%


FORM_DATA=cell(1,length(varargin));

for i=1:length(varargin)
	FORM_DATA{i}=cat(3,varargin{i}{:});
end
