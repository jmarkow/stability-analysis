% REQUIRES plotSpread for beeswarm plot
%
%
%

[options,dirs]=stan_preflight;
load custom_colormaps;

figs=stan_raster_nervecut(fee_map);
scaling_fun=@(x) (x/1.69)*5;

%tightfig(figs.y273);
%tightfig(figs.lpur72);

% scale in a manner that keeps time consistent
bird_names=fieldnames(figs);

for i=1:length(bird_names)

	ax=get(figs.(bird_names{i}),'currentAxes');
	xrange1=range(get(ax,'xlim'));

	set(figs.(bird_names{i}),'units','centimeters','position',[3 3 10 6.5],'paperpositionmode','auto');
	ax=findall(figs.(bird_names{i}),'type','axes');

	for j=1:length(ax)

		set(ax(j),'units','centimeters');
		pos=get(ax(j),'position');

		% width change

		new_width=scaling_fun(xrange1);
		width_change=new_width-pos(3);

		set(ax(j),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
	end

	markolab_multi_fig_save(figs.(bird_names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),['figure_4a_' bird_names{i}],'eps,png,fig,pdf');

end

