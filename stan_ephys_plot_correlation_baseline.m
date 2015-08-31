function fig=stan_ephys_plot(STATS)
% ephys stats plots
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collect spike and rms correlation

ci_inv=.01; % alpha level (/2 to get two-tailed)
ylimit_rounding=1e-1;
r_cutoff=.4; % below this value typically due to equipment failure
save_name='baseline_regression';
markersize=5;

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


fontsize=7;

% spike and rms plot

fig=figure();

ax(1)=subplot(2,1,1);
stan_plot_regress(plot_data.days_since(:),plot_data.spikes(:),plot_data.birdid(:),'markersize',markersize);
yh=ylabel('FR');
xlim([-.5 120])
set(yh,'position',get(yh,'position')+[.2 -.2 0])
set(gca,'XTick',[],'TickLength',[0 0],'YTick',[0 1],'fontsize',fontsize);

pos=get(ax(1),'position')
asp_ratio=pos(3)/pos(4);
width=.33;
new_width=width/asp_ratio;
h_offset=.28;
v_offset=.1;

new_axis(1)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);
stan_plot_regress(plot_data.days_since(:),plot_data.spikes(:),plot_data.birdid(:),'markersize',markersize);
xlim([-.2 20.2]);
ylim([0 1]);
set(gca,'XTick',[0 20],'YTick',[],'TickDir','out','TickLength',[ 0 0 ],'fontsize',fontsize);

ax(2)=subplot(2,1,2);
stan_plot_regress(plot_data.days_since(:),plot_data.rms(:),plot_data.birdid(:),'markersize',markersize);
yh=ylabel('RMS');
xlim([-.5 120]);
set(yh,'position',get(yh,'position')+[.2 -.2 0])
xlabel('Days');
set(gca,'TickLength',[0 0],'YTick',[0 1],'fontsize',fontsize);
ylim([0 1]);

pos=get(ax(2),'position')

new_axis(2)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);
stan_plot_regress(plot_data.days_since(:),plot_data.rms(:),plot_data.birdid(:),'markersize',markersize);
xlim([-.2 20.2]);
ylim([0 1]);
set(gca,'XTick',[0 20],'YTick',[],'TickDir','out','TickLength',[ 0 0 ],'fontsize',fontsize);

linkaxes(ax,'xy');
linkaxes(new_axis,'xy');
linkaxes([ax(:);new_axis(:)],'y');
ylim([0 1]);


%set(fig,'units','centimeters','position',[3 3 8 13],'paperpositionmode','auto');
%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),[save_name],'eps,fig,png,pdf','renderer','painters');
%save(fullfile(dirs.agg_dir,dirs.fig_dir,[ save_name '.mat']),'plot_data');

