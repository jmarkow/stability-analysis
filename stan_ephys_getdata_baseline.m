function BASELINE_DATA=stan_ephys_getdata_baseline()
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

tmp=dir(fullfile(dirs.agg_dir,dirs.control_dir,'*.mat'));
control_files={tmp(:).name};
tmp=dir(fullfile(dirs.agg_dir,dirs.nervecut_dir,'*.mat'));
nervecut_files={tmp(1:2:end).name};
all_files=[control_files nervecut_files];
file_idx=[ones(length(control_files),1);ones(length(nervecut_files),1)*2];
file_idx=file_idx==1;

BASELINE_DATA.dates={};
BASELINE_DATA.days_since={};
BASELINE_DATA.ntrials={};
BASELINE_DATA.spike_rate={};
BASELINE_DATA.spike_threshold={};

BASELINE_DATA.rms={};

padding_smps=round([options.padding-.1]*options.spike_fs);
rem_bird=[];

for i=1:length(all_files)

	data=struct();

	if file_idx(i)
		load(fullfile(dirs.agg_dir,dirs.control_dir,all_files{i}),'store');
	else
		load(fullfile(dirs.agg_dir,dirs.nervecut_dir,all_files{i}),'store');
	end

	% get the userdata

	bird_name=store(1).bird_id;	
	disp([bird_name])

	def_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,['defaults.txt']));
	user_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,[bird_name '.txt']));

	user_names=fieldnames(user_options);

	for j=1:length(user_names)
		def_options.(user_names{j})=user_options.(user_names{j});
	end

	user_options=def_options;

	% match motif, channel, compute stats

	motif_list={store(:).motif_name};
	
	if file_idx(i)
		motif_idx=strcmp(motif_list,user_options.motif_select);
	else
		motif_idx=strcmp(motif_list,user_options.motif_select1);
	end

	ch_idx=strcmp(lower(store(motif_idx).ch_list),lower(user_options.channel));

	% only use days with a sufficient number of trials

	data.dates=store(motif_idx).datenums(ch_idx,:);
	
	data.ntrials=cellfun(@(x) size(x,2),store(motif_idx).spikes.smooth_rate{ch_idx});
	data.spike_rate=store(motif_idx).spikes.smooth_rate{ch_idx};
	data.spike_threshold=cell(size(data.spike_rate));

	for j=1:length(data.spike_rate)
		[nsamples,ntrials]=size(data.spike_rate{j});
		data.spike_threshold{j}=zeros(1,ntrials);

		cur_spikes=store(motif_idx).spikes.threshold{ch_idx}{j};
		cur_trials=store(motif_idx).spikes.trial{ch_idx}{j};

		for k=1:ntrials
			data.spike_threshold{j}(k)=mean(cur_spikes(cur_trials==k));
		end
	end
	
	% get rms values

	data.rms=store(motif_idx).rms.data{ch_idx};

	% now collect spike threshold stats

	data.dates=store(motif_idx).datenums(ch_idx,:);
	data.dates(data.dates==0)=[];

	data_types=fieldnames(data);

	if isfield(user_options,'exclude')
		user_options.exclude(user_options.exclude>length(data.ntrials))=[];

		for j=1:length(data_types)
			tmp=size(data.(data_types{j}),2);
			if tmp==1
				data.(data_types{j})(user_options.exclude)=[];
			else
				data.(data_types{j})(:,user_options.exclude)=[];
			end
		end
	end

	sz_include=data.ntrials>=options.ephys_trial_limit;

	if ~any(sz_include)
		continue;
	end

	% remove pad from ephys data

	for j=1:length(data_types)		
		data.(data_types{j})=data.(data_types{j})(sz_include);

		if iscell(data.(data_types{j})) 
			[m,n]=size(data.(data_types{j}){1});
			if	(m>1&n>1)
				for k=1:length(data.(data_types{j}))
					data.(data_types{j}){k}=data.(data_types{j}){k}(padding_smps(1):end-padding_smps(2),:);
				end
			end
		end
	end

	params=store(motif_idx).spikes.parameters{ch_idx}{1}; 
	data.days_since=data.dates-data.dates(1);
	data_types=fieldnames(data); % update data types

	% store data

	for j=1:length(data_types)
		BASELINE_DATA.(data_types{j}){i}=data.(data_types{j});
	end

end

save(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_baseline_data.mat']),'BASELINE_DATA')
