function stan_analyze_nervecut_audio()
%
%
%

% get SAP features, compare pre/post for best match motif
%


save_dir='features';

% get options

[options,dirs]=stan_preflight;

key=stan_read_nervecut_audio;
bird_names={key(:).bird_name};
birds=unique(bird_names);

save_dir=fullfile(dirs.agg_dir,dirs.nervecut_audio_dir,save_dir);
if ~exist(save_dir,'dir'), mkdir(save_dir); end

for i=1:length(birds)

	idx=strcmp(bird_names,birds{i});

	cur_key=key(idx);
	motifs={cur_key(:).motif_name_short};

	% first hit is pre
	
	motif_idx=find(strcmp(cur_key(1).motif_name_short,motifs))

	if sum(motif_idx)<2
		error('Not enough motif hits!');
	end

	pre_key=cur_key(motif_idx(1));
	
	savefile=fullfile(save_dir,[ pre_key.bird_name '_audiofeatures.mat']);
	
	%if exist(savefile,'file'),continue; end

	post_key=cur_key(motif_idx(2));

	load(fullfile(dirs.agg_dir,dirs.template_dir,[pre_key.bird_name '_' pre_key.motif_name '.mat']),'template');
	precut_template=template;
	
	load(fullfile(dirs.agg_dir,dirs.template_dir,[pre_key.bird_name '_' post_key.motif_name '.mat']),'template');
	postcut_template=template;

	clearvars template;

	[shift,id]=stan_get_offset(precut_template.data,postcut_template.data,'fs',precut_template.fs,'audio_proc',1,'rms_tau',.05);
	shift_template=postcut_template.data;

	precut_len=length(precut_template.data);
	postcut_len=length(postcut_template.data);

	load(pre_key.filename,'audio');
	[nsamples,ntrials]=size(audio.data)

	pad_smps=round(options.padding*audio.fs);

	% trim pads

	audio.data=audio.data(pad_smps(1):end-pad_smps(2),:);

	% align

	if id==1
		audio.data=audio.data(shift:shift+postcut_len-1,:);
	end

	for j=1:ntrials
		[pre_tmp(j)]=zftftb_sap_score(audio.data(:,j),audio.fs);
	end

	feature_names=fieldnames(pre_tmp(1));
	feature_names(strcmp(feature_names,'spec_deriv'))=[];

	load(post_key.filename,'audio');
	[nsamples,ntrials]=size(audio.data)

	% trim pads
	
	audio.data=audio.data(pad_smps(1):end-pad_smps(2),:);

	% align

	if id==2
		audio.data=audio.data(shift:shift+precut_len-1,:);
	end

	for j=1:ntrials
		[post_tmp(j)]=zftftb_sap_score(audio.data(:,j),audio.fs);
	end

	for j=1:length(feature_names)
		pre_features.(feature_names{j})=cat(2,pre_tmp(:).(feature_names{j}));
		post_features.(feature_names{j})=cat(2,post_tmp(:).(feature_names{j}));
	end

	save(savefile,'pre_features','post_features');
	clearvars pre_tmp post_tmp pre_features post_features;

end
