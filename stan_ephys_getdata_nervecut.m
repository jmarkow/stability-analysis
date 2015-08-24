function NERVECUT_STATS=stan_compute_baseline()
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;
tmp=dir(fullfile(dirs.agg_dir,dirs.nervecut_dir,'*.mat'));
nervecut_files_pre={tmp(1:2:end).name};
nervecut_files_post={tmp(2:2:end).name};
template_key=stan_read_templates;

NERVECUT_STATS.rms=[];
NERVECUT_STATS.rms_ci=[];
NERVECUT_STATS.spikes=[];
NERVECUT_STATS.spikes_ci=[];
NERVECUT_STATS.days_since=[];
NERVECUT_STATS.birdid=[];

bootstrap_ci=[ 1-options.bootstrap_alpha/2 options.bootstrap_alpha/2  ]

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


	% align audio

	load(fullfile(dirs.agg_dir,dirs.template_dir,[ bird_name '_' user_options.motif_select1 '.mat']),'template');
	precut_template=template;

	load(fullfile(dirs.agg_dir,dirs.template_dir,[ bird_name '_' user_options.motif_select2 '.mat']),'template');
	postcut_template=template;

	[shift,shift_id]=stan_get_offset(precut_template.data,postcut_template.data,'fs',precut_template.fs,'audio_proc',1,'rms_tau',.05);
	shift_template=postcut_template.data;

	precut_len=length(precut_template.data);
	postcut_len=length(postcut_template.data);

	% round shift to spike_fs
	
	shift_t=shift/precut_template.fs;
	shift=round(shift_t*options.spike_fs);
	

	% match motif, channel, compute stats

	motif_list={store(:).motif_name};	
	motif_idx=strcmp(motif_list,user_options.motif_select1);
	ch_idx=strcmp(lower(store(motif_idx).ch_list),lower(user_options.channel));

	% need motif1 name and motif2 name, grab from template directory to get offset

	% only use days with a sufficient number of trials

	sz=cellfun(@(x) size(x,2),store(motif_idx).spikes.smooth_rate{ch_idx});
	template_spikes_mu=cellfun(@(x) median(x,2),store(motif_idx).spikes.smooth_rate{ch_idx},'uniformoutput',false);
	template_rms_mu=cellfun(@(x) median(x,2),store(motif_idx).rms.data{ch_idx},'uniformoutput',false);
	
	template_spikes_mu=cat(2,template_spikes_mu{:});
	template_rms_mu=cat(2,template_rms_mu{:});

	% get final day from the pre-data, use to compare with the post-cut

	dates=store(motif_idx).datenums(ch_idx,:);
	dates(dates==0)=[];

	if isfield(user_options,'exclude')
		user_options.exclude(user_options.exclude>length(sz))=[];
		sz(user_options.exclude)=[];
		template_spikes_mu(:,user_options.exclude)=[];
		template_rms_mu(:,user_options.exclude)=[];
		dates(user_options.exclude)=[];
	end

	sz_include=sz>=options.ephys_trial_limit;
	
	if ~any(sz_include)
		continue;
    end
   
	% now take final day
	% take an extra 100 ms 
    
	params=store(motif_idx).spikes.parameters{ch_idx}{1}; 
	padding_smps=round((options.padding*params.smooth_fs))

	template_spikes_mu=template_spikes_mu(:,sz_include);
	template_rms_mu=template_rms_mu(:,sz_include);

	dates=dates(sz_include);
	
	template_spikes=template_spikes_mu(:,end);
	template_rms=template_rms_mu(:,end);
	template_date=dates(end);	

	% now get post-cut data, compare with template

	clearvars store template_spikes_mu template_rms_mu sz dates;

	load(fullfile(dirs.agg_dir,dirs.nervecut_dir,nervecut_files_post{i}),'store');

	motif_list={store(:).motif_name};
	motif_idx=strcmp(motif_list,user_options.motif_select2);
	ch_idx=strcmp(lower(store(motif_idx).ch_list),lower(user_options.channel));

	% only use days with a sufficient number of trials

	sz=cellfun(@(x) size(x,2),store(motif_idx).spikes.smooth_rate{ch_idx});
	
	spikes_mu=cellfun(@(x) median(x,2),store(motif_idx).spikes.smooth_rate{ch_idx},'uniformoutput',false);
	rms_mu=cellfun(@(x) median(x,2),store(motif_idx).rms.data{ch_idx},'uniformoutput',false);
	
	spikes_mu=cat(2,spikes_mu{:});
	rms_mu=cat(2,rms_mu{:});

	% get final day from the pre-data, use to compare with the post-cut

	idx=1:length(store(motif_idx).spikes.smooth_rate{ch_idx});
	dates=store(motif_idx).datenums(ch_idx,:);
	dates(dates==0)=[];

	if isfield(user_options,'exclude2')
		user_options.exclude2(user_options.exclude2>length(sz))=[];
		user_options.exclude2
		sz(user_options.exclude2)=[];
		spikes_mu(:,user_options.exclude2)=[];
		rms_mu(:,user_options.exclude2)=[];
		dates(user_options.exclude2)=[];
		idx(user_options.exclude2)=[];
	end

	sz_include=sz>=options.ephys_trial_limit;
	
	if ~any(sz_include)
		continue;
    end

	% compare mu to template mu

	spikes_mu=spikes_mu(:,sz_include);
	rms_mu=rms_mu(:,sz_include);

	dates=dates(sz_include);
	idx=idx(sz_include);

	spikes_corr=zeros(size(spikes_mu,2),1);
	spikes_ci=zeros(size(spikes_mu,2),2);

	template_spikes=template_spikes(padding_smps(1):end-padding_smps(2));
	spikes_mu=spikes_mu(padding_smps(1):end-padding_smps(2),:);
	template_rms=template_rms(padding_smps(1):end-padding_smps(2));
	rms_mu=rms_mu(padding_smps(1):end-padding_smps(2),:);

	% shift

	if shift==0
		shift=1;
	end

	if shift_id==1 
		template_spikes=template_spikes(shift:shift+size(spikes_mu,1)-1);
		template_rms=template_rms(shift:shift+size(rms_mu,1)-1);
	else
		spikes_mu=spikes_mu(shift:shift+length(template_spikes)-1,:);
		rms_mu=rms_mu(shift:shift+length(template_rms)-1,:);
	end

	for j=1:size(spikes_mu,2)

		x1=zscore(spikes_mu(:,j));
		x2=zscore(template_spikes);

		norm_fact=sqrt(sum(x1.^2)*sum(x2.^2));
		spikes_corr(j)=max(xcorr(x1,x2))./norm_fact; % convert to corr coefficient

		% bootstrap the correlation
		bootval=zeros(1,options.nbootstraps);
		bootdata=store(motif_idx).spikes.smooth_rate{ch_idx}{idx(j)};

		bootdata=bootdata(padding_smps(1):end-padding_smps(2),:);
		
		if shift_id==2
			bootdata=bootdata(shift:shift+length(template_spikes-1),:);
		end

		ntrials=size(bootdata,2);
		trial_pool=1:ntrials;

		for k=1:options.nbootstraps
			new_trials=randsample(trial_pool,ntrials,true);
			tmp=zscore(median(bootdata(:,new_trials)'));
			norm_fact=sqrt(sum(tmp.^2)*sum(x2.^2));
			bootval(k)=max(xcorr(tmp,x2))./norm_fact;
		end

		spikes_ci(j,:)=quantile(bootval,bootstrap_ci);

	end

	rms_corr=zeros(size(rms_mu,2),1);
	rms_ci=zeros(size(rms_mu,2),2);

	for j=1:size(rms_mu,2)

		x1=zscore(rms_mu(:,j));
		x2=zscore(template_rms);

		norm_fact=sqrt(sum(x1.^2)*sum(x2.^2));
		rms_corr(j)=max(xcorr(x1,x2))./norm_fact; % convert to corr coefficient

		% bootstrap the correlation
		bootval=zeros(1,options.nbootstraps);
		bootdata=store(motif_idx).rms.data{ch_idx}{idx(j)};

		bootdata=bootdata(padding_smps(1):end-padding_smps(2),:);
		
		if shift_id==2
			bootdata=bootdata(shift:shift+length(template_rms-1),:);
		end


		ntrials=size(bootdata,2);
		trial_pool=1:ntrials;

		for k=1:options.nbootstraps
			new_trials=randsample(trial_pool,ntrials,true);
			tmp=zscore(median(bootdata(:,new_trials)'));
			norm_fact=sqrt(sum(tmp.^2)*sum(x2.^2));
			bootval(k)=max(xcorr(tmp,x2))./norm_fact;
		end

		rms_ci(j,:)=quantile(bootval,bootstrap_ci);

	end

	dates=dates-dates(1);

	NERVECUT_STATS.spikes=[NERVECUT_STATS.spikes spikes_corr'];
	NERVECUT_STATS.spikes_ci=[NERVECUT_STATS.spikes_ci spikes_ci'];
	
	NERVECUT_STATS.rms=[NERVECUT_STATS.rms rms_corr'];
	NERVECUT_STATS.rms_ci=[NERVECUT_STATS.rms_ci rms_ci'];
	
	NERVECUT_STATS.days_since=[NERVECUT_STATS.days_since dates];
	NERVECUT_STATS.birdid=[NERVECUT_STATS.birdid ones(size(spikes_corr'))*i];

end

save(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_nervecut_stats.mat']),'NERVECUT_STATS')
