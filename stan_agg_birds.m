function stan_agg_birds(SAVE_DIR)
% searches for mat files and associated log files, processes and copies to target directory
%

if nargin<1 | isempty(SAVE_DIR)
	SAVE_DIR='mu_agg';
end

max_depth=4;

% get options

[options,dirs]=stan_preflight;
filenames=robofinch_dir_recurse(pwd,'*stan_withinbird.mat',max_depth);

% all proc data goes into the same directory

SAVE_DIR=fullfile(dirs.agg_dir,SAVE_DIR);
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end

for i=1:length(filenames)

	disp([filenames(i).name]);
	copyfile(filenames(i).name,SAVE_DIR);

end
