function LFP_DATA=stan_ephys_getdata_lfp()
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
		pli_data=cellfun(@(x) abs(mean(exp(1j.*x(:,1:options.lfp_trial_limit)),2)),ang_data,'uniformoutput',0);

		pli_data=cat(2,pli_data{:});	
		[nsamples,ntrials]=size(pli_data);

		thresh=1e-3/nsamples;

		% cat the pli data

		n=options.lfp_trial_limit;
		r=pli_data.*n;
		z=(r.^2)./n;
		pval=exp(sqrt(1+4*n+4.*(n^2-r.^2))-(1+2*n));

		% get bootstrap distribution
	
		fft_base=fft(filt_data{1}(:,1:options.lfp_trial_limit),nsamples);
		fft_mag=repmat(abs(fft_base),[1 1 nboots]);
		%fft_ang=angle(fft_base);

		pli_boot=zeros(nsamples,nboots);

		fft_rand=fft(randn(nsamples,options.lfp_trial_limit,nboots));
		fft_ang_rand=angle(fft_rand);
		synth_sig=real(ifft(fft_mag.*exp(1j.*fft_ang_rand)));

		for	k=1:nboots
			pli_boot(:,k)=abs(mean(exp(1j.*angle(hilbert(synth_sig(:,:,k)))),2));
		end

		pval=mean(repmat(pli_data(:,1),[1 nboots])<=repmat(max(pli_boot),[nsamples 1]),2);
		pli_t0=pval(:,1);
		pli_mask=pli_t0<=thresh;
		
		%pli_mask=pli_data(:,1)>=.5;
		pli_idx=find(pli_mask);
		pli_idx(pli_idx<padding_smps(1))=[];
		pli_idx(pli_idx>nsamples-padding_smps(2))=[];

		% get only contiguous regions???

		if isempty(pli_idx)
			continue;
		end

		[nsamples,ntrials]=size(pli_data);

		ang_diff=zeros(1,ntrials);
	
		%template=angle(mean(exp(1j.*ang_data{1}),2));
		template=angle(hilbert(mean(zscore(filt_data{1}),2)));
		template=template(pli_idx);

		for k=1:ntrials
		
			compare=mean(zscore(filt_data{k}),2);
			compare=angle(hilbert(compare));
			
			%compare=angle(mean(exp(1j.*ang_data{k}),2));
			compare=compare(pli_idx);	
			ang_diff(k)=median(abs(angle(exp(1j.*(template-compare)))));

		end
	
		%ang_diff

		LFP_DATA.ang_diff=[LFP_DATA.ang_diff ang_diff(2:end)];
		LFP_DATA.days_since=[LFP_DATA.days_since days_since(2:end)];
		LFP_DATA.bird_id=[LFP_DATA.bird_id ones(size(ang_diff(2:end)))*i];
		LFP_DATA.channel_id=[LFP_DATA.channel_id ones(size(ang_diff(2:end)))*counter];	
		counter=counter+1;

	end

	% remove pad from ephys data

end

save(fullfile(dirs.agg_dir,dirs.fig_dir,['ephys_baseline_lfp_data.mat']),'LFP_DATA','options')
