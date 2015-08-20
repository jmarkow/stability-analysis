function fig=stan_control_raster(COLORS)
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

idx=regexp(option_names,'control_raster(\d+)_(\w+)','tokens'); 

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
	
	ctrl(i).path=fullfile(dirs.agg_dir,ctrl(i).path);
	
	load(ctrl(i).path,'store');

	% get template
	
	motif_list={store(:).motif_name};
	idx=find(strcmp(ctrl(i).motif_name,motif_list));

	load(fullfile(dirs.agg_dir,dirs.template_dir,...
		[ store(idx).bird_id '_' store(idx).motif_name ]),'template','parameters');

	% parse bird and motif_name from 

	% pad out the template
	
	pad_smps=round(template.fs*options.padding);

	template.data=[repmat(template.data(1),[pad_smps(1) 1]);template.data(:);repmat(template.data(end),[pad_smps(2) 1])]

	[spect.s,spect.f,spect.t]=zftftb_pretty_sonogram(template.data,template.fs,'filtering',300,'clipping',[-3 2],...
		'len',70,'overlap',69.5,'zeropad',0);
	
	% collect spike rasters, throw together into figure

	ch_list=store(idx).ch_list;
	ch_idx=find(strcmp(lower(ctrl(i).channel),lower(ch_list)));

	plot_spikes=struct();
	plot_trials=ctrl(i).trials;
	
	for j=1:length(ctrl(i).days)

		plot_spikes(j).trial=store(idx).spikes.trial{ch_idx}{ctrl(i).days(j)};
		plot_spikes(j).times=store(idx).spikes.times{ch_idx}{ctrl(i).days(j)};
		plot_spikes(j).threshold=store(idx).spikes.threshold{ch_idx}{ctrl(i).days(j)};
		plot_spikes(j).fs=store(idx).spikes.parameters{ch_idx}{ctrl(i).days(j)}.spike_fs;
		plot_spikes(j).smooth_rate=store(idx).spikes.smooth_rate{ch_idx}{ctrl(i).days(j)};

		r=corr(zscore(plot_spikes(j).smooth_rate(100:end-100,:)));

		% find good case, mean r > .4
	
		r_vec=mean(r,2);
		[~,good_trial]=max(r_vec);
		r_check=r(good_trial,:);
		outliers=find(r_check<(mean(r_check)-.3*(std(r_check))));
		fr_check=mean(plot_spikes(j).smooth_rate);
		outliers2=find(fr_check>(mean(fr_check)+1*(std(fr_check))));
		outliers3=find(fr_check<(mean(fr_check)-1*(std(fr_check))));
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

	% spike raster


	% remove outliers
	
	fig.(store(idx).bird_id)=figure();
	ax=stan_plot_raster(spect,plot_spikes,'spike_height',.5,'spike_width',.01,'plot_trials',plot_trials,'colors',COLORS);
	
	if isfield(ctrl(i),'xlim') & ~isempty(ctrl(i).xlim)
		xlimits=ctrl(i).xlim;
		xlim([ctrl(i).xlim]);
	end



	xlimits=get(ax(end),'xlim');
	new_xlimits=[ round(xlimits/time_order)*time_order ];
	span=range(new_xlimits);
	xlim([new_xlimits]);
	set(ax(3),'xtick',[]);
	set(ax(2),'ylim',[1 150]);
	set(ax(3),'ylim',[1 150]);
	h=line([new_xlimits(1) new_xlimits(1)+.2],[170 170],'linewidth',1.5,'color','k');
	set(h,'clipping','off');
	
	for i=1:length(ax)
		pos=get(ax(i),'position');
		set(ax(i),'position',[ .025 pos(2) .95 pos(4)]);
	end	
	%set(ax(end),'xtick',new_xlimits','xticklabel',new_xlimits-new_xlimits(1));
	%set(fig,'units','inches','position',[2 2 (span/2)*4.5 4],'paperpositionmode','auto')
	%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),['control_raster_' store(idx).bird_id ],'eps,png,fig,pdf','renderer','painters');

end
