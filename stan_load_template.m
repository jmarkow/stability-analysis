function [TEMPLATE,FS]=stan_load_legacy(FILENAME)
%
%
%
%

vars=whos('-file',FILENAME);

var_names={vars(:).name};
islegacy=any(strcmp(var_names,'TEMPLATE'));

if islegacy
	load(FILENAME,'TEMPLATE');
	FS=25e3;
else
	load(FILENAME,'template');
	FS=template.fs;
	TEMPLATE=template.data;
end


