function stan_control_raster(COLORS)
% searches for mat files and associated log files, processes and copies to target directory
%


if nargin<1
	COLORS='jet';
end
time_order=1e-2;
% get options

[options,dirs]=stan_preflight;

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
	
	ctrl(i).precut_path=fullfile(dirs.agg_dir,ctrl(i).precut_path);
	ctrl(i).postcut_path=fullfile(dirs.agg_dir,ctrl(i).postcut_path);
	
	% load data

	load(ctrl(i).precut_path,'store');
	precut_store=store;
	clear store;

	load(ctrl(i).postcut_path,'store');
	postcut_store=store;
	clear store;

	% get templates
	

	motif_list={precut_store(:).motif_name};
	precut_idx=find(strcmp(ctrl(i).precut_motif_name,motif_list));

	motif_list={postcut_store(:).motif_name};
	postcut_idx=find(strcmp(ctrl(i).postcut_motif_name,motif_list));

	load(fullfile(dirs.agg_dir,dirs.template_dir,...
		[precut_store(precut_idx).bird_id '_' precut_store(precut_idx).motif_name '.mat']),'template');
	pad_smps=round(template.fs*options.padding);
	template.data=[repmat(template.data(1),[pad_smps(1) 1]);template.data(:);repmat(template.data(end),[pad_smps(2) 1])]
	precut_template=template;

	clear template;

	load(fullfile(dirs.agg_dir,dirs.template_dir,...
		[postcut_store(postcut_idx).bird_id '_' postcut_store(postcut_idx).motif_name '.mat']),'template');
	pad_smps=round(template.fs*options.padding);
	template.data=[repmat(template.data(1),[pad_smps(1) 1]);template.data(:);repmat(template.data(end),[pad_smps(2) 1])]

	postcut_template=template;
	clear template;

	% get offset between templates

	shift=stan_get_offset(precut_template.data,postcut_template.data,'fs',precut_template.fs,'audio_proc',1);

	% smooth rate, raster

	precut_ch_idx=find(strcmp(precut_store(precut_idx).ch_list,ctrl(i).channel));
	postcut_ch_idx=find(strcmp(postcut_store(postcut_idx).ch_list,ctrl(i).channel));
	precut_spikes=struct();

	[spect.s,spect.f,spect.t]=zftftb_pretty_sonogram(precut_template.data,...
		precut_template.fs,'filtering',300,'clipping',[-3 2],'len',70,'overlap',69.5,'zeropad',0);

	plot_spikes=struct();

	for j=1:length(ctrl(i).precut_days)
		plot_spikes(j).times=precut_store(precut_idx).spikes.times{precut_ch_idx}{ctrl(i).precut_days(j)};
		plot_spikes(j).trial=precut_store(precut_idx).spikes.trial{precut_ch_idx}{ctrl(i).precut_days(j)};
		plot_spikes(j).threshold=precut_store(precut_idx).spikes.parameters{precut_ch_idx}{ctrl(i).precut_days(j)};
		plot_spikes(j).fs=precut_store(precut_idx).spikes.parameters{precut_ch_idx}{ctrl(i).precut_days(j)}.spike_fs;
		plot_spikes(j).smooth_rate=precut_store(precut_idx).spikes.smooth_rate{precut_ch_idx}{ctrl(i).precut_days(j)};
		r=corr(zscore(plot_spikes(j).smooth_rate(100:end-100,:)));

		% find good case, mean r > .4
	
		r_vec=mean(r,2);
		[~,good_trial]=max(r_vec);
		r_check=r(good_trial,:);
		outliers=find(r_check<(mean(r_check)-1.5*(std(r_check))));
		fr_check=mean(plot_spikes(j).smooth_rate);
		outliers2=find(fr_check>(mean(fr_check)+2*(std(fr_check))));
		outliers3=find(fr_check<(mean(fr_check)-2*(std(fr_check))));
		outliers=unique([outliers(:);outliers2(:);outliers3(:)]);

		outlier_idx=[];

		for k=1:length(outliers)
			outlier_idx=[outlier_idx find(plot_spikes(j).trial==outliers(k))];
		end

		plot_spikes(j).trial(outlier_idx)=[];
		plot_spikes(j).times(outlier_idx)=[];

		uniq_trials=unique(plot_spikes(j).trial);

		for k=1:length(uniq_trials)
			% correct trials 

			correction=sum(outliers<uniq_trials(k));
			new_idx=(plot_spikes(j).trial==uniq_trials(k));
			plot_spikes(j).trial(new_idx)=plot_spikes(j).trial(new_idx)-correction;

		end


	end

	plot_spect=repmat(spect,[length(ctrl(i).precut_days) 1]);
	plot_trials=repmat(ctrl(i).precut_trials,[length(ctrl(i).precut_days) 1]);
	
	shift_template=postcut_template.data;

	%if shift<0
	%	shift_template(1:-shift)=[];
	%elseif shift>0
	%	shift_template(end-shift:end)=[];
	%end

	[spect.s,spect.f,spect.t]=zftftb_pretty_sonogram(shift_template,...
		precut_template.fs,'filtering',300,'clipping',[-3 2],'len',70,'overlap',69.5,'zeropad',0);
	spect.t
	offset=length(plot_spikes);

	for j=1:length(ctrl(i).postcut_days)

		plot_spikes(j+offset).times=postcut_store(postcut_idx).spikes.times{postcut_ch_idx}{ctrl(i).postcut_days(j)};
		plot_spikes(j+offset).trial=postcut_store(postcut_idx).spikes.trial{postcut_ch_idx}{ctrl(i).postcut_days(j)};
		plot_spikes(j+offset).threshold=postcut_store(postcut_idx).spikes.parameters{postcut_ch_idx}{ctrl(i).postcut_days(j)};
		plot_spikes(j+offset).fs=postcut_store(postcut_idx).spikes.parameters{postcut_ch_idx}{ctrl(i).postcut_days(j)}.spike_fs;
		plot_spikes(j+offset).smooth_rate=postcut_store(postcut_idx).spikes.smooth_rate{postcut_ch_idx}{ctrl(i).postcut_days(j)};
		r=corr(zscore(plot_spikes(j+offset).smooth_rate(100:end-100,:)));

		% find good case, mean r > .4
	
		r_vec=mean(r,2);
		[~,good_trial]=max(r_vec);
		r_check=r(good_trial,:);
		outliers=find(r_check<(mean(r_check)-1.5*(std(r_check))));

		fr_check=mean(plot_spikes(j+offset).smooth_rate);
		
		outliers2=find(fr_check>(mean(fr_check)+2*(std(fr_check))));
		outliers3=find(fr_check<(mean(fr_check)-2*(std(fr_check))));
		outliers=unique([outliers(:);outliers2(:);outliers3(:)])

		outlier_idx=[];

		for k=1:length(outliers)
			outlier_idx=[outlier_idx find(plot_spikes(j+offset).trial==outliers(k))];
		end

		plot_spikes(j+offset).trial(outlier_idx)=[];
		plot_spikes(j+offset).times(outlier_idx)=[];

		uniq_trials=unique(plot_spikes(j+offset).trial);

		for k=1:length(uniq_trials)
			% correct trials 

			correction=sum(outliers<uniq_trials(k));
			new_idx=(plot_spikes(j+offset).trial==uniq_trials(k));
			plot_spikes(j+offset).trial(new_idx)=plot_spikes(j+offset).trial(new_idx)-correction;

		end


	end

	plot_spect=[plot_spect;repmat(spect,[length(ctrl(i).postcut_days) 1])];
	plot_trials=[plot_trials;repmat(ctrl(i).postcut_trials,[length(ctrl(i).postcut_days) 1])];
	
	fig=figure();
	ax=stan_plot_raster_horizontal(plot_spect,plot_spikes,'plot_trials',plot_trials,'colors',COLORS,'spike_width',.01,'spike_height',.5);
	
	set(gca,'xtick',[]);
	xlabel('');
	h=line([.25 .45],[170 170],'linewidth',1.5,'color','k');
	set(h,'clipping','off');
	xlimits=xlim();;
	xlim([xlimits(1)+.25 xlimits(end)-.25]);
	
	time_range=range(xlim());

	set(fig,'units','centimeters','position',[2 2 (time_range/1.5)*18 8],'paperpositionmode','auto');
	% rms, raster

	markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),['nervecut_raster_' precut_store(precut_idx).bird_id ],'eps,png,fig,pdf','renderer','painters');


end
