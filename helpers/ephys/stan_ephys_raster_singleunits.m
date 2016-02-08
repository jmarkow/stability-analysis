function fig=stan_plot_singleunits(COLORS)
%
%
%
%

[options,dirs]=stan_preflight;

if nargin<1 | isempty(COLORS)
	COLORS='jet';
end

% get listing of single units
%

tmp=dir(fullfile(dirs.agg_dir,dirs.su_plot_dir));
unit_dirs={};
unit_names={};

for	i=1:length(tmp)
	if tmp(i).isdir & tmp(i).name(1)~='.'
		unit_dirs{end+1}=fullfile(dirs.agg_dir,dirs.su_plot_dir,tmp(i).name);
		unit_names{end+1}=tmp(i).name;
	end
end

for i=1:length(unit_dirs)

	listing=dir(fullfile(unit_dirs{i},'*.mat'));
	nplots=length(listing);

	plot_spikes=struct();

	tmp=regexp(listing(1).name,'^(\w+)\_','tokens');
	bird_id=tmp{1}{1};
	date_num=[];

	for j=1:length(listing)

		% get cluster number

		tmp=regexp(listing(j).name,'\_(\d+-\d+-\d+)\_','tokens');
		date_num(j)=datenum(tmp{1}{1});

		tmp=regexp(listing(j).name,'\_cl(\d+)\.mat','tokens');
		clust_num=str2num(tmp{1}{1});
		listing(j).name

		load(fullfile(unit_dirs{i},listing(j).name),'cluster');
		plot_spikes(j).times=cluster.times{clust_num};
		plot_spikes(j).trial=cluster.trials{clust_num};
		plot_spikes(j).windows=cluster.windows{clust_num};
		plot_spikes(j).fs=cluster.parameters.fs;

	end

	% grab the appropriate template, correct for padding

	load(fullfile(dirs.agg_dir,dirs.template_dir,[ bird_id '_motif1_padding.mat' ]),'template');

	pad_smps=round(template.fs*options.padding_su);
	template.data=[repmat(template.data(1),[pad_smps(1) 1]);...
		template.data(:);...
		repmat(template.data(end),[pad_smps(2) 1])];

	[spect.s,spect.f,spect.t]=zftftb_pretty_sonogram(template.data,...
		template.fs,'filtering',300,'clipping',[-3 2],'len',70,'overlap',69.5,'zeropad',0);

	% put raster in leftmost column, expand

	fig.(unit_names{i})=figure();
	cols=4;
	ax=stan_ephys_plot_raster(spect,plot_spikes,...
		'colors',COLORS,'spike_width',.01,'spike_height',.5,'columns',cols);

	for j=1:length(ax)
		pos=get(ax(j),'position');
		set(ax(j),'position',[pos(1:2) pos(3)*(cols-1) pos(4)]);
	end

	% add day labels

	for j=2:length(ax)
		days_since=daysdif(date_num(1),date_num(j-1));
		ylabel(ax(j),sprintf('Day %i',days_since));
		set(ax(j),'ylim',[0 100]);
	end

	% add cluster waveforms

	length(ax)


	counter=cols*2;
	ax2=[];
	tvec=[-.0005:1/cluster.parameters.interpolate_fs:.0005];
	for j=1:length(plot_spikes)

		length(ax)*2-1
		ax2(j)=subplot(length(ax)*2-1,cols,[counter counter+cols]);
		counter=counter+cols*2;
		mu=mean(plot_spikes(j).windows');
		sig=std(plot_spikes(j).windows');
		markolab_shadeplot(tvec*1e3,[mu-sig;mu+sig],[1 0 0],[0 0 0]);
		hold on;
		plot([tvec*1e3],mean(plot_spikes(j).windows'),'k-');
		set(ax2(j),'ytick',[],'xtick',[]);


	end

	linkaxes(ax2,'xy');
	xlim([-.5 .5]);
	ylim([-250 250]);

	ylimits=ylim();
	%set(gca,'YTick',ylimits,'XTick',[-.5 0 .5],'TickLength',[ 0 0 ],'TickDir','out');
	ylabel(ax2(end),'500 \muV');
	xlabel(ax2(end),'1 ms');
	%h=line([-.6 -.6],[ylimits(1) ylimits(1)+50]);
	%set(h,'clipping','off');




end
