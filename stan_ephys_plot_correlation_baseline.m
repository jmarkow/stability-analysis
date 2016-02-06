function fig=stan_ephys_plot_correlation_baseline(STATS)
% ephys stats plots
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collect spike and rms correlation

ci_inv=.01; % alpha level (/2 to get two-tailed)
ylimit_rounding=1e-1;
r_cutoff=.4; % due to equipment failure
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

b(1,:)=stan_plot_regress(plot_data.days_since(:),plot_data.spikes(:),plot_data.birdid(:),'markersize',markersize);
yh=ylabel('FR');
set(yh,'position',get(yh,'position')+[5 -.2 0])
set(gca,'XTick',[],'TickLength',[0 0],'YTick',[0 1],'fontsize',fontsize);

pos=get(ax(1),'position')
asp_ratio=pos(3)/pos(4);
width=.3;
new_width=width/asp_ratio;
h_offset=.25;
v_offset=.27;

new_axis(1)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);
%set(ax(1),'position',[ pos(1) pos(2) pos(3)-new_width pos(4)-v_offset-.05])
stan_plot_regress(plot_data.days_since(:),plot_data.spikes(:),plot_data.birdid(:),'markersize',markersize);
ylim([0 1]);
set(gca,'XTick',[],'YTick',[],'TickDir','out','TickLength',[ 0 0 ],'fontsize',fontsize);

ax(2)=subplot(2,1,2);
b(2,:)=stan_plot_regress(plot_data.days_since(:),plot_data.rms(:),plot_data.birdid(:),'markersize',markersize);
yh=ylabel('RMS');
set(yh,'position',get(yh,'position')+[5 -.2 0])
set(gca,'TickLength',[0 0],'YTick',[0 1],'fontsize',fontsize);
ylim([0 1]);

pos=get(ax(2),'position')

new_axis(2)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);
%set(ax(2),'position',[ pos(1) pos(2) pos(3)-new_width pos(4)-v_offset-.05])
stan_plot_regress(plot_data.days_since(:),plot_data.rms(:),plot_data.birdid(:),'markersize',markersize);
ylim([0 1]);
set(gca,'XTick',[],'YTick',[],'TickDir','out','TickLength',[ 0 0 ],'fontsize',fontsize);

linkaxes(ax,'xy');
linkaxes(new_axis,'xy');
linkaxes([ax(:);new_axis(:)],'y');
ylim([0 1]);

linkaxes(ax,'x');
set(ax(1),'xlim',[0 40]);
xh=xlabel(ax(end),'Days');
set(ax(end),'xtick',get(ax(end),'xlim'));
set(xh,'position',get(xh,'position')+[0 .01 0]);

linkaxes(new_axis,'x');
set(new_axis(1),'xlim',[0 120]);
set(new_axis(end),'xtick',get(new_axis(end),'xlim'));

fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,'fig1_regression.txt'),'w+');
fprintf(fid,'FR R=%g\nRMS R=%g\nN(points)=%g\nN(birds)=%g',b(1,1),b(2,1),length(plot_data.days_since(:)),length(unique(plot_data.birdid(:))));
fclose(fid);
