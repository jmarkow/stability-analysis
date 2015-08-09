function stan_control_raster()
% searches for mat files and associated log files, processes and copies to target directory
%

options_name='options.txt';
dirs_name='dirs.txt';

time_order=1e-2;
% get options

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);
options=stan_read_options(fullfile(cur_path,options_name));
dirs=stan_read_options(fullfile(cur_path,dirs_name));

% all proc data goes into the same directory

option_names=fieldnames(options);

% which options specify control raster

idx=regexp(option_names,'nervecut_raster(\d+)_(\w+)','tokens'); 

% gather all rasters specified in options

ctrl=struct();

for i=1:length(idx)
	if ~isempty(idx{i})
		raster_number=str2num(idx{i}{1}{1});
		ctrl(raster_number).(idx{i}{1}{2})=options.(option_names{i});
	end
end

template_key=stan_read_templates();

for i=1:length(ctrl)
	
	%ctrl(i).path=fullfile(dirs.agg_dir,ctrl(i).path);
	
	%load(ctrl(i).path,'store');

	% get template
	
	%motif_list={store(:).motif_name};
	%idx=find(strcmp(ctrl(i).motif_name,motif_list));

	ctrl(i)

	% set ticks for pub

	% smooth rate, raster




	% rms, raster


end
