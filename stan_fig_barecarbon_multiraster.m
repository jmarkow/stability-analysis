function stan_fig6a
%
%
%
%

[options,dirs]=stan_preflight;
load custom_colormaps;
%fig=stan_raster_multi(fee_map);
fig=stan_ephys_raster_multi_raw(fee_map);
set(fig,'units','centimeters','position',[3 3 14 16],'paperpositionmode','auto');

ax=findall(fig,'type','axes');
set(ax(:),'xtick',[]);
xlimits=get(ax(1),'xlim');

h=line([xlimits(1) xlimits(1)+.2],[-10 -10]);
set(h,'clipping','off');

markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),['figure_s1'],'eps,png,fig,pdf');
