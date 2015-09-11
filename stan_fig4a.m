function stan_fig6a
%
%
%
%

[options,dirs]=stan_preflight;
load custom_colormaps;
fig=stan_raster_singleunits(fee_map);
scaling_fun=@(x) (x/1.5)*4;

cell_names=fieldnames(fig);

for i=1:length(cell_names)
	set(fig.(cell_names{i}),'units','centimeters','position',[3 3 8 7],'paperpositionmode','auto');

	% axes 3,4 are spikes, 5 sonogram
	
	ax=findall(fig.(cell_names{i}),'type','axes');

	for j=3:5
		set(ax(j),'units','centimeters');
		pos=get(ax(j),'position');
		xrange=range(get(ax(j),'xlim'));
		set(ax(j),'position',[ pos(1:2) scaling_fun(xrange) pos(4)]);

	end

	ylimits=get(ax(3),'ylim');
	h=line([0 .2],[ylimits(2)+5 ylimits(2)+5],'parent',ax(3));
	set(h,'clipping','off');
	set(ax(3),'xtick',[]);

	markolab_multi_fig_save(fig.(cell_names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),['figure_7a_' cell_names{i}],'eps,png,fig,pdf');
end

%ax=findall(fig,'type','axes');
%set(ax(:),'xtick',[]);
%xlimits=get(ax(1),'xlim');

%h=line([xlimits(1) xlimits(1)+.2],[-10 -10]);
%set(h,'clipping','off');

