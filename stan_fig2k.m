function stan_fig2k(SINGLEUNIT_DATA)
% Generates Fig 2k, example single unit rasters
%
% note:  requires spikoclust_raster...

scaling_fun=@(x) (x/1.5)*3;
bird_id=fieldnames(SINGLEUNIT_DATA);

for i=1:length(bird_id)

	fig.(bird_id{i})=figure();
	cols=4;
	ax=stan_ephys_plot_raster(SINGLEUNIT_DATA.(bird_id{i}).spect,SINGLEUNIT_DATA.(bird_id{i}).spikes,...
		'spike_width',.01,'spike_height',.5,'columns',cols);
	interpolate_fs=SINGLEUNIT_DATA.(bird_id{i}).spikes(1).fs*8;

	for j=1:length(ax)
		pos=get(ax(j),'position');
		set(ax(j),'position',[pos(1:2) pos(3)*(cols-1) pos(4)]);
	end

	% add day labels

	for j=2:length(ax)
		days_since=daysdif(SINGLEUNIT_DATA.(bird_id{i}).date_num(1),...
			SINGLEUNIT_DATA.(bird_id{i}).date_num(j-1));
		ylabel(ax(j),sprintf('Day %i',days_since));
		set(ax(j),'ylim',[0 100]);
	end

	% add cluster waveforms

	length(ax)

	counter=cols*2;
	ax2=[];
	tvec=[-.0005:1/interpolate_fs:.0005];

	for j=1:length(SINGLEUNIT_DATA.(bird_id{i}).spikes)

		length(ax)*2-1
		ax2(j)=subplot(length(ax)*2-1,cols,[counter counter+cols]);
		counter=counter+cols*2;
		mu=mean(SINGLEUNIT_DATA.(bird_id{i}).spikes(j).windows');
		sig=std(SINGLEUNIT_DATA.(bird_id{i}).spikes(j).windows');
		markolab_shadeplot(tvec*1e3,[mu-sig;mu+sig],[1 0 0],[0 0 0]);
		hold on;
		plot([tvec*1e3],mean(SINGLEUNIT_DATA.(bird_id{i}).spikes(j).windows'),'k-');
		set(ax2(j),'ytick',[],'xtick',[]);

	end

	linkaxes(ax2,'xy');
	xlim([-.5 .5]);
	ylim([-250 250]);

	ylimits=ylim();
	ylabel(ax2(end),'500 \muV');
	xlabel(ax2(end),'1 ms');
end

cell_names=fieldnames(fig);

for i=1:length(cell_names)
	set(fig.(cell_names{i}),'units','centimeters','position',[3 3 8 7],'paperpositionmode','auto');

	% axes 3,4 are spikes, 5 sonogram

	ax=findall(fig.(cell_names{i}),'type','axes');

	n=ceil(length(ax)/2);

	for j=n:length(ax)
		set(ax(j),'units','centimeters');
		pos=get(ax(j),'position');
		xrange=range(get(ax(j),'xlim'));
		newpos=scaling_fun(xrange);
		storepos=newpos+pos(1);
		set(ax(j),'position',[ pos(1:2) newpos pos(4)]);

	end

	ylimits=get(ax(n),'ylim');
	h=line([0 .2],[ylimits(2)+10 ylimits(2)+10],'parent',ax(n));
	set(h,'clipping','off');
	set(ax(n),'xtick',[]);

	for j=1:n-1
		axis(ax(j),'off');
		set(ax(j),'units','centimeters');
		pos=get(ax(j),'position')
		set(ax(j),'position',[ storepos+.15 pos(2) .75 pos(4) ]);
	end

	ylimits=get(ax(1),'ylim');
	xlimits=get(ax(1),'xlim');

	h=line([xlimits(1)-.1 xlimits(1)-.1],[ylimits(1)-21 ylimits(1)+79],'parent',ax(1));
	set(h,'clipping','off');
	h2=line([xlimits(1)-.1 xlimits(1)+.4],[ylimits(1)-21 ylimits(1)-21],'parent',ax(1));
	set(h2,'clipping','off');

end
