function store_data=stan_compute_baseline()
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

tmp=dir(fullfile(dirs.agg_dir,dirs.nervecut_dir,'*.mat'));
nervecut_files_pre={tmp(1:2:end).name};
nervecut_files_post={tmp(2:2:end).name};
template_key=stan_read_templates;

store_data.rms_corr_mu=[];
store_data.rms_corr_ci=[];
store_data.spikes_corr_mu=[];
store_data.spikes_corr_ci=[];
store_data.days_since=[];
store_data.birdid=[];

for i=1:length(nervecut_files_pre)

	load(fullfile(dirs.agg_dir,dirs.nervecut_dir,nervecut_files_pre{i}),'store');
	
	nervecut_files_pre{i}

	% get the userdata

	if ~isfield(store,'bird_id')
		continue;
	end

	bird_name=store(1).bird_id

	def_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,['defaults.txt']));
	user_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,[bird_name '.txt']));

	user_names=fieldnames(user_options);

	for j=1:length(user_names)
		def_options.(user_names{j})=user_options.(user_names{j});
	end

	user_options=def_options;

	% match motif, channel, compute stats

	motif_list={store(:).motif_name};	
	motif_idx=strcmp(motif_list,user_options.motif_select1);
	ch_idx=strcmp(lower(store(motif_idx).ch_list),lower(user_options.channel));

	% need motif1 name and motif2 name, grab from template directory to get offset

	% only use days with a sufficient number of trials

	sz=cellfun(@(x) size(x,2),store(motif_idx).spikes.smooth_rate{ch_idx});
	mu=cellfun(@(x) median(x,2),store(motif_idx).spikes.smooth_rate{ch_idx},'uniformoutput',false);
	mu=cat(2,mu{:});

	% get final day from the pre-data, use to compare with the post-cut

	dates=store(motif_idx).datenums(ch_idx,:);
	dates(dates==0)=[];

	if isfield(user_options,'exclude')
		user_options.exclude(user_options.exclude>length(sz))=[];
		sz(user_options.exclude)=[];
		mu(:,user_options.exclude)=[];
		dates(user_options.exclude)=[];
	end

	sz_include=sz>=options.ephys_trial_limit;
	
	if ~any(sz_include)
		continue;
    end
   
	% now take final day
	% take an extra 100 ms 
    
	params=store(motif_idx).spikes.parameters{ch_idx}{1}; 
	padding_smps=round((options.padding-[.1 .1])*params.smooth_fs)

	mu=mu(:,sz_include);
	dates=dates(sz_include);
	template_data=mu(:,end);
	template_date=dates(end);

	% now get post-cut data, compare with template

	clearvars store;

	nervecut_files_post{i}
	load(fullfile(dirs.agg_dir,dirs.nervecut_dir,nervecut_files_post{i}),'store');

	motif_list={store(:).motif_name}
	motif_idx=strcmp(motif_list,user_options.motif_select2)
	ch_idx=strcmp(lower(store(motif_idx).ch_list),lower(user_options.channel))

	% only use days with a sufficient number of trials

	sz=cellfun(@(x) size(x,2),store(motif_idx).spikes.smooth_rate{ch_idx});
	mu=cellfun(@(x) median(x,2),store(motif_idx).spikes.smooth_rate{ch_idx},'uniformoutput',false);
	mu=cat(2,mu{:});

	% get final day from the pre-data, use to compare with the post-cut

	idx=1:length(store(motif_idx).spikes.smooth_rate{ch_idx});
	dates=store(motif_idx).datenums(ch_idx,:);
	dates(dates==0)=[];

	if isfield(user_options,'exclude2')
		user_options.exclude2(user_options.exclude2>length(sz))=[];
		user_options.exclude2
		sz(user_options.exclude2)=[];
		mu(:,user_options.exclude2)=[];
		dates(user_options.exclude2)=[];
		idx(user_options.exclude2)=[];
	end

	sz_include=sz>=options.ephys_trial_limit;
	
	if ~any(sz_include)
		continue;
    end

	% compare mu to template mu

	mu=mu(:,sz_include);
	dates=dates(sz_include);
	idx=idx(sz_include);

	corr=zeros(size(mu,2),1);
	ci=zeros(size(mu,2),2);
	
	for j=1:size(mu,2)
		x1=zscore(mu(padding_smps(1):end-padding_smps(2),j));
		x2=zscore(template_data(padding_smps(1):end-padding_smps(2)));
		norm_fact=sqrt(sum(x1.^2)*sum(x2.^2));
		corr(j)=max(xcorr(x1,x2))./norm_fact; % convert to corr coefficient

		% bootstrap the correlation
		bootval=zeros(1,options.nbootstraps);
		bootdata=store(motif_idx).spikes.smooth_rate{ch_idx}{idx(j)};
		ntrials=size(bootdata,2);
		trial_pool=1:ntrials;

		for k=1:options.nbootstraps
			new_trials=randsample(trial_pool,ntrials,true);
			tmp=median(bootdata(:,new_trials)');
			tmp=zscore(tmp(padding_smps(1):end-padding_smps(2)));
			norm_fact=sqrt(sum(tmp.^2)*sum(x2.^2));
			bootval(k)=max(xcorr(tmp,x2))./norm_fact;
		end

		ci(j,:)=prctile(bootval,[99.5 .5]);

	end

	%figure();plot(dates,corr,'k.--');
	%datetick('x');
	%pause();
	
	dates=dates-dates(1);

	corr=corr';
	ci=ci'

	store_data.corr_mu=[store_data.corr_mu corr];
	store_data.corr_ci=[store_data.corr_ci ci];
	store_data.days_since=[store_data.days_since dates];
	store_data.birdid=[store_data.birdid ones(size(corr))*i];

end

