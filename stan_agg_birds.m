function stan_agg_birds()
% searches for mat files and associated log files, processes and copies to target directory
%

options_name='options.txt';
dry_run=0;
save_dir='mu_agg';
max_depth=4;

% get options

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);
options=stan_read_options(fullfile(cur_path,options_name));

filenames=robofinch_dir_recurse(pwd,'*stan_withinbird.mat',max_depth);

% all proc data goes into the same directory

save_dir=fullfile(options.agg_dir,save_dir);
if ~exist(save_dir,'dir'), mkdir(save_dir); end

for i=1:length(filenames)

	disp([filenames(i).name]);
	copyfile(filenames(i).name,save_dir);

end
