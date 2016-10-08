function stan_fig2h(SPIKEDATA)
% Generates Fig 2h, multi-unit spike raster
%

bird_id=fieldnames(SPIKEDATA);
for i=1:length(bird_id)
	figs.(bird_id{i})=figure();
	stan_ephys_plot_raster(SPIKEDATA.(bird_id{i}).spect,SPIKEDATA.(bird_id{i}).spikes);
end

scaling_fun=@(x) (x/1.69)*3.6;

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
end


ax=findall(figs.lpur72,'type','axes');

for i=1:length(ax)

	set(ax(i),'units','centimeters');
	pos=get(ax(i),'position');

	% width change

	new_width=scaling_fun(xrange2);
	width_change=new_width-pos(3);

	set(ax(i),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
end
