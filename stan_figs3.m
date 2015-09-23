function stan_fig3a
%
%
%
%

[options,dirs]=stan_preflight;
scaling_fun=@(x) (x/1)*3;

load custom_colormaps;
fig=stan_sonograms_nervecut(fee_map);
bird_names=fieldnames(fig);

for i=1:length(bird_names)

	ax=get(fig.(bird_names{i}),'CurrentAxes');
	xlimits=get(ax,'xlim');
	xrange=range(xlimits);
	
	set(fig.(bird_names{i}),'units','centimeters','position',[4 4 6 3.5],'paperpositionmode','auto');

	ax=findall(fig.(bird_names{i}),'type','axes');

	for j=1:length(ax)	
		set(ax(j),'units','centimeters');
		pos=get(ax(j),'position')
		new_width=scaling_fun(xrange)
		width_change=new_width-pos(3)
		set(ax(j),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
	end



	markolab_multi_fig_save(fig.(bird_names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),['figure_3a_' bird_names{i} ],'eps,png,fig,pdf');

end



