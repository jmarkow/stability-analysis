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

	for j=1:length(ctrl(i).days)

		plot_spikes(j).trial=store(idx).spikes.trial{ch_idx}{ctrl(i).days(j)};
		plot_spikes(j).times=store(idx).spikes.times{ch_idx}{ctrl(i).days(j)};
		plot_spikes(j).threshold=store(idx).spikes.threshold{ch_idx}{ctrl(i).days(j)};
		plot_spikes(j).fs=store(idx).spikes.parameters{ch_idx}{ctrl(i).days(j)}.spike_fs;

	end

	% spike raster

	plot_trials=ctrl(i).trials;


	fig=figure('units','inches','paperpositionmode','auto','position',[2 2 3 4]);
	ax=stan_songalign_raster(spect,plot_spikes,'spike_height',.5,'spike_width',1,'plot_trials',plot_trials);
	
	if isfield(ctrl(i),'xlim') & ~isempty(ctrl(i).xlim)
		xlimits=ctrl(i).xlim;
		xlim([ctrl(i).xlim]);
	end

	xlimits=get(ax(end),'xlim')
	new_xlimits=[ round(xlimits/time_order)*time_order ]
	xlim([new_xlimits]);

	set(ax(end),'xtick',new_xlimits','xticklabel',new_xlimits-new_xlimits(1));


	% set ticks for pub

	% smooth rate, raster




	% rms, raster


end
