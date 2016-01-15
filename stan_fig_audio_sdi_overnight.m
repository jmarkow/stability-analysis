function stan_fig_audio_overnight();
%%%
[options,dirs]=stan_preflight;
[fig,stats]=stan_audio_overnight_sdi;
scaling_fun=@(x) (x/1)*3;

names=fieldnames(fig);

for i=1:length(names)
	if length(strfind(names{i},'zoom'))>0
		set(fig.(names{i}),'units','centimeters','position',[10 10 3.7 6],'paperpositionmode','auto');
		ax=findall(fig.(names{i}),'type','axes')
		for j=1:length(ax)
			set(ax(j),'XTick',[],'YTick',[]);
		end

		% scale axis
    %set(gca,'units','centimeters');
    %pos=get(gca,'position');
    %set(gca,'position',[pos(1:2) scaling_fun(diff(xlim())) pos(4)]);

	elseif strcmp(names{i},'sdi')
		set(fig.(names{i}),'units','centimeters','position',[10 10 14 7],'paperpositionmode','auto');
	end

	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'sdiovernight_' names{i} ],'eps,png,fig',...
		'renderer','painters');

end
