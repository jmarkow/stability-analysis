function stan_ephys_getdata_su_all
%

[~,dirs]=system('find . -type d -name ''to_copy''');
hits=regexp(dirs,'\n','split');
hits(end)=[];

basedir=pwd;

for i=1:length(hits)
	cd(fullfile(basedir,hits{i}(2:end)));
	stan_ephys_getdata_su;
end

