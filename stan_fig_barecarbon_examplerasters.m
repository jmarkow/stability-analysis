function stan_fig1c
%
%
%
%

[options,dirs]=stan_preflight;
load custom_colormaps;
figs=stan_ephys_raster_control(fee_map);
scaling_fun=@(x) (x/1.69)*3.6;

%tightfig(figs.y273);
%tightfig(figs.lpur72);

% scale in a manner that keeps time consistent

ax=get(figs.y273,'currentAxes');
xrange1=range(get(ax,'xlim'));

ax=get(figs.lpur72,'currentAxes');
xrange2=range(get(ax,'xlim'));

% rescale axes (can't make figure small enough here)

set(figs.y273,'units','centimeters','position',[3 3 6 6.5],'paperpositionmode','auto');
set(figs.lpur72,'units','centimeters','position',[3 3 6 6.5],'paperpositionmode','auto');

ax=findall(figs.y273,'type','axes');

for i=1:length(ax)

	set(ax(i),'units','centimeters');
	pos=get(ax(i),'position');

	% width change

	new_width=scaling_fun(xrange1);
	width_change=new_width-pos(3);

	set(ax(i),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
	%set(ax(i),'xcolor',get(gcf,'color'),'ycolor',get(gcf,'color'))
end


ax=findall(figs.lpur72,'type','axes');

for i=1:length(ax)

	set(ax(i),'units','centimeters');
	pos=get(ax(i),'position');

	% width change

	new_width=scaling_fun(xrange2);
	width_change=new_width-pos(3);

	set(ax(i),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
	%set(ax(i),'xcolor',get(gcf,'color'),'ycolor',get(gcf,'color'))
end

markolab_multi_fig_save(figs.y273,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_2a_y273'],'eps,png,fig,pdf');
markolab_multi_fig_save(figs.lpur72,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_2a_lpur72'],'eps,png,fig,pdf');
