function fig=stan_sonograms_nervecut(COLORS)
%
%
%
%

if nargin<1
	COLORS='jet';
end

[options,dirs]=stan_preflight;

tmp=dir(fullfile(dirs.agg_dir,dirs.nervecut_dir,'*.mat'));
nervecut_listing={tmp(:).name};
nervecut_listing=nervecut_listing(1:2:end);

for i=1:length(nervecut_listing)

	% get bird name, extract the relevant templates, make a nice align pair of stacked sonograms

	birdid=stan_read_filename(nervecut_listing{i})

	% get this bird's userdata

	def_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,[ 'defaults.txt' ]));	
	user_options=stan_read_options(fullfile(dirs.agg_dir,dirs.user_dir,[ birdid '.txt']));

	user_names=fieldnames(user_options);

	for j=1:length(user_names)
		def_options.(user_names{j})=user_options.(user_names{j});
	end

	user_options=def_options;

	load(fullfile(dirs.agg_dir,dirs.template_dir,...
		[ birdid '_' user_options.motif_select1 '.mat' ]),'template');

	precut_template=template;
	clearvars template;

	load(fullfile(dirs.agg_dir,dirs.template_dir,...
		[ birdid '_' user_options.motif_select2 '.mat' ]),'template');

	postcut_template=template;
	clearvars template;

	[shift,id]=stan_get_offset(precut_template.data,postcut_template.data,'fs',precut_template.fs,'audio_proc',1,'rms_tau',.05);
	shift_template=postcut_template.data;

	precut_len=length(precut_template.data);
	postcut_len=length(postcut_template.data);

	if id==1
		precut_template.data=precut_template.data(shift:shift+length(postcut_template.data)-1);
	else
		postcut_template.data=postcut_template.data(shift:shift+length(precut_template.data)-1);
	end

	[precut.s,precut.f,precut.t]=zftftb_pretty_sonogram(precut_template.data,precut_template.fs,'filtering',300,...
		'clipping',[-3 2],'len',70,'overlap',69.5,'zeropad',0);

	[postcut.s,postcut.f,postcut.t]=zftftb_pretty_sonogram(postcut_template.data,postcut_template.fs,'filtering',300,...
		'clipping',[-3 2],'len',70,'overlap',69.5,'zeropad',0);

	fig.(birdid)=figure();
	ax(1)=subplot(2,1,1);
	imagesc(precut.t,precut.f/1e3,precut.s)
	%ylabel('FS (kHz)');
	box off;
	set(gca,'YTick',[],'XTick',[],'TickLength',[0 0]);
	axis xy;
	ax(2)=subplot(2,1,2);
	pos=get(ax(2),'position');
	%set(ax(2),'position',[ pos(1) pos(2)+.1 pos(3:4)]);
	imagesc(postcut.t,postcut.f/1e3,postcut.s);
	box off;
	set(gca,'XTick',[],'YTick',[]);
	linkaxes(ax,'xy');
	ylim([0 10]);
	axis xy;
	colormap(COLORS);
	h=line([0 .2],[-1.5 -1.5],'linewidth',1.5,'color','k');
	set(h,'clipping','off');

	%set(fig,'position',[300 300 500*(diff(xlim())/2) 200],'paperpositionmode','auto');
	%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),[ birdid '_nervecut_sonograms'],'eps,png,fig,pdf');

end
