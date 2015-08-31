function stan_fig2a
%
%
%
%

[options,dirs]=stan_preflight;
load custom_colormaps;
figs=stan_ephys_plot_lfp_example(fee_map);

scaling_fun=@(x) (x/1.5)*5;

%tightfig(figs.y273);
%tightfig(figs.lpur72);

% scale in a manner that keeps time consistent

ax=get(figs.lhp33,'currentAxes');
xrange1=range(get(ax,'xlim'));

ax=get(figs.rm7,'currentAxes');
xrange2=range(get(ax,'xlim'));

% rescale axes (can't make figure small enough here)

set(figs.lhp33,'units','centimeters','position',[3 3 6 6.5],'paperpositionmode','auto');
set(figs.rm7,'units','centimeters','position',[3 3 6 6.5],'paperpositionmode','auto');

ax=findall(figs.lhp33,'type','axes');

for i=1:length(ax)

	set(ax(i),'units','centimeters');
	pos=get(ax(i),'position');

	% width change
	
	new_width=scaling_fun(xrange1);
	width_change=new_width-pos(3);

	set(ax(i),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
	%set(ax(i),'xcolor',get(gcf,'color'),'ycolor',get(gcf,'color'))
end

set(ax(1),'ylim',[-5 23],'xtick',[],'ytick',[]);
xlimits=get(ax(1),'xlim');
h=line([xlimits(1) xlimits(1)+.2],[-7 -7],'parent',ax(1));
h2=line([xlimits(1)-.1*xrange1 xlimits(1)-.1*xrange1],[-5 0],'parent',ax(1));
set(h,'clipping','off');
set(h2,'clipping','off');

ax=findall(figs.rm7,'type','axes');

for i=1:length(ax)

	set(ax(i),'units','centimeters');
	pos=get(ax(i),'position');

	% width change
	
	new_width=scaling_fun(xrange2);
	width_change=new_width-pos(3);

	set(ax(i),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
	%set(ax(i),'xcolor',get(gcf,'color'),'ycolor',get(gcf,'color'))
end

set(ax(1),'ylim',[-5 23],'xtick',[],'ytick',[]);
xlimits=get(ax(1),'xlim');
h=line([xlimits(1) xlimits(1)+.2],[-7 -7],'parent',ax(1));
h2=line([xlimits(1)-.1*xrange2 xlimits(1)-.1*xrange2],[-5 0],'parent',ax(1));
set(h,'clipping','off');
set(h2,'clipping','off');


markolab_multi_fig_save(figs.lhp33,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_5a_lfp33'],'eps,png,fig,pdf');
markolab_multi_fig_save(figs.rm7,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_5a_rm7'],'eps,png,fig,pdf');



