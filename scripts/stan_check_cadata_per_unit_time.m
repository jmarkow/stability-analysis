
stan_cadata_drift_analyze_pertime_ave;

[~,idx]=sort(all_time);
fig=figure();stan_plot_regress(all_time(idx)*24,all_ca(idx),all_id(idx))
xlim([7 17])
ylabel('Correlation to n-1 ave (r)')
xlabel('Time (hr)')
set(gca,'TickLength',[0 0])
title('')
ylim([.2 .9])

[r,p]=corrcoef(all_time,all_ca);
