function plot_data=stan_ephys_plot(STATS)
% ephys stats plots
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collect spike and rms correlation

ci_inv=.01; % alpha level (/2 to get two-tailed)
ylimit_rounding=1e-1;
r_cutoff=.4; % below this value typically due to equipment failure
save_name='baseline_regression';

% first column is days since 1, second is regression value, third is bird ID
% remove points where x=0 (by definition == 1, artifacts appear to have r<.4)

plot_data.rms=[];
plot_data.spikes=[];
plot_data.birdid=[];
plot_data.days_since=[];

for i=1:length(STATS.rms_corr)
	npoints=size(STATS.spikes_corr{i},2)-1;
	plot_data.spikes=[plot_data.spikes STATS.spikes_corr{i}(1,2:end)];
	plot_data.rms=[plot_data.rms STATS.rms_corr{i}(1,2:end)];
	plot_data.days_since=[plot_data.days_since STATS.days_since{i}(2:end)];
	plot_data.birdid=[plot_data.birdid ones(1,npoints)*i];
end

plot_data

artifacts=plot_data.spikes<=r_cutoff;
artifacts2=isnan(plot_data.spikes);

to_del=find(artifacts|artifacts2)

plot_data.spikes(to_del)=[];
plot_data.days_since(to_del)=[];
plot_data.birdid(to_del)=[];
plot_data.rms(to_del)=[];

% spike and rms plot

fig=figure();

ax(1)=subplot(2,1,1);
stan_plot_regress(plot_data.days_since(:),plot_data.spikes(:),plot_data.birdid(:),'markersize',20);
ylabel('FR correlation (R)');
set(gca,'XTick',[],'TickLength',[0 0],'YTick',[.5:.25:1]);

ax(2)=subplot(2,1,2);
stan_plot_regress(plot_data.days_since(:),plot_data.rms(:),plot_data.birdid(:),'markersize',20);
ylabel('RMS correlation (R)');
xlabel('Days');
set(gca,'TickLength',[0 0],'YTick',[.5:.25:1]);
ylim([.5 1]);
linkaxes(ax,'xy');

set(fig,'position',[200 200 150 300],'paperpositionmode','auto');
markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),[save_name],'eps,fig,png,pdf','renderer','painters');
save(fullfile(dirs.agg_dir,dirs.fig_dir,[ save_name '.mat']),'plot_data');

