
stan_cadata_drift_analyze_pertime_ave;

%%
[~,idx]=sort(all_time);
fig=figure();stan_plot_regress(all_time(idx)*24,all_ca(idx),all_id(idx))
xlim([7 17])
ylabel('Correlation to n-1 ave (r)')
xlabel('Time (hr)')
set(gca,'TickLength',[0 0],'YTick',[.2 .9])
title('')
ylim([.2 .9])
set(fig,'position',[300 700 230 210]);

%%
[~,idx]=sort(all_trial);
fig2=figure();stan_plot_regress(all_trial(idx),all_ca(idx),all_id(idx))
xlim([0 400])
ylabel('Correlation to n-1 ave (r)')
xlabel('Trial num.')
set(gca,'TickLength',[0 0])
title('')
ylim([.2 .9])

%%
[r,p]=corrcoef(all_time,all_ca);
[r2,p2]=corrcoef(all_trial,all_ca);
