function store_data=stan_compute_baseline()
%
% stability analysis--baseline data
% first take all of the data from control

options_name='options.txt';
dirs_name='dirs.txt';

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);
options=stan_read_options(fullfile(cur_path,options_name));
dirs=stan_read_options(fullfile(cur_path,dirs_name));

tmp=dir(fullfile(dirs.agg_dir,dirs.control_dir,'*.mat'));
control_files={tmp(:).name};
tmp=dir(fullfile(dirs.agg_dir,dirs.nervecut_dir,'*.mat'));
nervecut_files={tmp(1:2:end).name};
all_files=[control_files nervecut_files];
file_idx=[ones(length(control_files),1);ones(length(nervecut_files),1)*2];
file_idx=file_idx==1;
store_data=[];

for i=1:length(all_files)

	if file_idx(i)
		load(fullfile(dirs.agg_dir,dirs.control_dir,all_files{i}),'store');
	else
		load(fullfile(dirs.agg_dir,dirs.nervecut_dir,all_files{i}),'store');
	end

	all_files{i}

	% get the userdata

	if ~isfield(store,'bird_id')
		continue;
	end


	bird_name=store(1).bird_id;
	def_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,['defaults.txt']));
	user_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,[bird_name '.txt']));

	user_names=fieldnames(user_options);

	for i=1:length(user_names)
		def_options.(user_names{i})=user_options.(user_names{i});
	end

	user_options=def_options;

	% match motif, channel, compute stats

	motif_list={store(:).motif_name};
	if file_idx(i)
		motif_idx=strcmp(motif_list,user_options.motif_select)
	else
		motif_idx=strcmp(motif_list,user_options.motif1_select);
	end

	ch_idx=strcmp(lower(store(motif_idx).ch_list),lower(user_options.channel));

	% only use days with a sufficient number of trials

	sz=cellfun(@(x) size(x,2),store(motif_idx).spikes.smooth_rate{ch_idx});
	mu=cellfun(@(x) median(x,2),store(motif_idx).spikes.smooth_rate{ch_idx},'uniformoutput',false);
	mu=cat(2,mu{:});

	dates=store(motif_idx).datenums(ch_idx,:);
	dates(dates==0)=[];
	
	if isfield(user_options,'exclude')
		sz(user_options.exclude)=[];
		mu(:,user_options.exclude)=[];
		dates(user_options.exclude)=[];
	end

	sz_include=sz>=user_options.trial_limit;
	
	if ~any(sz_include)
		continue;
	end

	mu=mu(:,sz_include);
	corrmat=corr(zscore(mu(100:end-100,:)));
	dates=dates(:,sz_include);

	dates_diff=dates-dates(1);
	store_data=[store_data;[ dates_diff(:) corrmat(1,:)' ]];

end


tmp=dir(fullfile(dirs.agg_dir,dirs.nervecut_dir,'*.mat'));
nervecut_files={tmp(:).name};

figure();plot(store_data(:,1),store_data(:,2),'k.');
