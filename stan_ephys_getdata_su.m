function stan_ephys_get_singleunit()
%
%
%
%

% scan for single unit files

[options,dirs]=stan_preflight;
files=robofinch_dir_recurse(pwd,'sua_channels*.mat');

for i=1:length(files)

	files(i).name

	% scan the metadata
	
	[pathname,filename,ext]=fileparts(files(i).name);

	su_info=stan_read_metadata_su(fullfile(pathname,'..','cellinfo.txt'));

	idx=0;
	
	new_filename=[ su_info.birdid '_' su_info.date ...
		'_idx' num2str(idx) '_ch' num2str(su_info.channel) '_cl' num2str(su_info.cluster) '.mat' ];
	new_fullpath=fullfile(dirs.agg_dir,dirs.su_dir,new_filename);

	while exist(new_fullpath,'file')

		idx=idx+1;
		
		new_filename=[ su_info.birdid '_' su_info.date ...
		'_idx' num2str(idx) '_ch' num2str(su_info.channel) '_cl' num2str(su_info.cluster) '.mat' ];
		new_fullpath=fullfile(dirs.agg_dir,dirs.su_dir,new_filename);

	end

	load(files(i).name,'cluster','proc_data');
	save(new_fullpath,'cluster','proc_data','su_info');

end
