function stan_cadata_schnitzplot()
%
%
%
%
%

[options,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.fig_dir,'AwesomeNess.mat'));

newdata=stan_cadata_format(data1,data2,data3,data4,data5);

fig.localnorm=figure();
stan_cadata_sortmat(newdata,'scaling','l','sort_day',1,'fig_nrows',3,'fig_row',1,'chk_day',1);
stan_cadata_sortmat(newdata,'scaling','l','sort_day',3,'fig_nrows',3,'fig_row',2,'chk_day',1);
stan_cadata_sortmat(newdata,'scaling','l','sort_day',5,'fig_nrows',3,'fig_row',3,'chk_day',1);
ax=findall(gcf,'type','axes')
linkaxes(ax,'x')
xlim([1 2.4])
title(ax(end-2),'Within day and roi scaling');
ylabel(ax(end),'Day 1');
ylabel(ax(end-5),'Day 3');
ylabel(ax(end-10),'Day 5');
set(ax(1),'xtick',xlim(),'TickLength',[0 0]);
xlabel(ax(1),'Time (s)');

fig.roinorm=figure();
stan_cadata_sortmat(newdata,'scaling','r','sort_day',1,'fig_nrows',3,'fig_row',1,'chk_day',1);
stan_cadata_sortmat(newdata,'scaling','r','sort_day',3,'fig_nrows',3,'fig_row',2,'chk_day',1);
stan_cadata_sortmat(newdata,'scaling','r','sort_day',5,'fig_nrows',3,'fig_row',3,'chk_day',1);
ax=findall(gcf,'type','axes')
linkaxes(ax,'x')
xlim([1 2.4])
title(ax(end-2),'Within roi across day scaling');
ylabel(ax(end),'Day 1');
ylabel(ax(end-5),'Day 3');
ylabel(ax(end-10),'Day 5');
set(ax(1),'xtick',xlim(),'TickLength',[0 0]);
xlabel(ax(1),'Time (s)');

fig.sortnorm=figure();
stan_cadata_sortmat(newdata,'scaling','s','sort_day',1,'fig_nrows',3,'fig_row',1,'chk_day',1);
stan_cadata_sortmat(newdata,'scaling','s','sort_day',3,'fig_nrows',3,'fig_row',2,'chk_day',1);
stan_cadata_sortmat(newdata,'scaling','s','sort_day',5,'fig_nrows',3,'fig_row',3,'chk_day',1);
ax=findall(gcf,'type','axes')
linkaxes(ax,'x')
xlim([1 2.4])
title(ax(end-2),'Within roi sort day');
ylabel(ax(end),'Day 1');
ylabel(ax(end-5),'Day 3');
ylabel(ax(end-10),'Day 5');
set(ax(1),'xtick',xlim(),'TickLength',[0 0]);
xlabel(ax(1),'Time (s)');

names=fieldnames(fig);

for i=1:length(names)
	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'schnitzplot_' names{i} ],'eps,png,fig,pdf');
end
