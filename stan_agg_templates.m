function stan_agg_templates()
%
%
%

template_file='template_data.mat';
dirs_name='dirs.txt';

disp('Collecting files...');

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);
dirs=stan_read_options(fullfile(cur_path,dirs_name));

temp_files=robofinch_dir_recurse(dirs.data_dir,template_file,4);
save_dir=fullfile(dirs.agg_dir,dirs.template_dir);

% now split and get the first directory for all files

first_dir=cell(1,length(temp_files));

if ~exist(save_dir,'dir'), mkdir(save_dir); end

for i=1:length(temp_files)

	disp([temp_files(i).name]);

	%tokens=regexp(all_files(i).name,filesep,'split');
	[pathname,filename,ext]=fileparts(temp_files(i).name);

	% grab configuration
	
	parameters=robofinch_read_config(fullfile(pathname,'robofinch_parameters.txt'));

	%ntokens=length(regexp(DIR,filesep,'split')); % first token after DIR

	% take two directories above path
	tokens=regexp(pathname,filesep,'split');
	use_tokens=tokens(2:end-2);

	motif_name=tokens{end};
	bird_name=tokens{end-(temp_files(i).depth-2)}

	load(temp_files(i).name,'template');

	new_filename=fullfile(save_dir,[ bird_name '_' motif_name '.mat' ]);
	save(new_filename,'template','parameters');
	
	%copyfile(temp_files(i).name,new_filename);

end

