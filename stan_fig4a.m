function stan_fig6a
%
%
%
%

[options,dirs]=stan_preflight;
load custom_colormaps;
fig=stan_raster_singleunits(fee_map);
scaling_fun=@(x) (x/1.5)*3;

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

	markolab_multi_fig_save(fig.(cell_names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),['figure_7a_' cell_names{i}],'eps,png,fig,pdf');
end

%ax=findall(fig,'type','axes');
%set(ax(:),'xtick',[]);
%xlimits=get(ax(1),'xlim');

%h=line([xlimits(1) xlimits(1)+.2],[-10 -10]);
%set(h,'clipping','off');

