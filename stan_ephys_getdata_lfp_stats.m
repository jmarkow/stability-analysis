function LFP_DATA=stan_ephys_getdata_lfp_stats()
%
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;
nboots=1e3;

tmp=dir(fullfile(dirs.agg_dir,dirs.lfp_dir,'*.mat'));
lfp_files={tmp(:).name};

LFP_DATA.dates=[];
LFP_DATA.days_since=[];
LFP_DATA.ang_diff=[];
LFP_DATA.channel_id=[];
LFP_DATA.bird_id=[];
counter=1;
padding_smps=round([options.padding_lfp+.1]*options.lfp_fs);

% filter setting

%[n,Wn,beta,ftype]=kaiserord([10 25 35 50],[0 1 0],[.01 .05 .01],options.lfp_fs);
%b=fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
%a=1;

%[b,a]=ellip(4,.2,40,[25 35]/(options.lfp_fs/2),'bandpass');
[b,a]=sfield_filt_coeffs(options.lfp_fs,2);

for i=1:length(lfp_files)

	data=struct();

	load(fullfile(dirs.agg_dir,dirs.lfp_dir,lfp_files{i}),'store');

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
	motif_idx=strcmp(motif_list,user_options.motif_select);

	% only use days with a sufficient number of trials

	ch_list=store(motif_idx).ch_list;
	dist={};

	for j=1:length(ch_list)

		% agg by bird

		data.dates=store(motif_idx).datenums(j,:);	
		data.ntrials=cellfun(@(x) size(x,2),store(motif_idx).lfp.data{j});
		data.lfp=store(motif_idx).lfp.data{j};

		% now collect spike threshold stats

		data.dates=store(motif_idx).datenums(j,:);
		data.dates(data.dates==0)=[];
		data_types=fieldnames(data);

		sz_include=data.ntrials>=options.lfp_trial_limit;

		if ~any(sz_include)
			continue;
		end

		t0=min(find(sz_include));

		data.dates=data.dates(sz_include);
		days_since=data.dates-min(data.dates);

		filt_data=cellfun(@(x) filtfilt(b,a,x),data.lfp(sz_include),'uniformoutput',0);
		ang_data=cellfun(@(x) angle(hilbert(x)),filt_data,'uniformoutput',0);

		%template=angle(mean(exp(1j.*ang_data{1}),2));
		template=ang_data{1};
		dist=cellfun(@(x) stan_angdist(template,x),ang_data,'uniformoutput',0)

		idx=find(triu(ones(size(dist{1})),1));
		tmp1=dist{1}(idx);

		for k=1:length(dist)
			idx=find(triu(ones(size(dist{k})),1));
			tmp2=dist{k}(idx);
			pop=(tmp2-mean(tmp1))/std(tmp1);
			ci=bootci(1e3,{@mean,pop},'type','per')
			pause();
		end

	end

	% remove pad from ephys data

end


