function [figs,figs_stats]=stan_plot_barecarbon_ca_timecourse()
%
%
%

% upsample?
interp_factor=5;
ndays=4;
swarm_offset=1;
nboots=1e2;
filewrite=false;

[options,dirs]=stan_preflight;

peakstats=load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats_peaktime.mat'),'stats');
peakstats=peakstats.stats;

% plot time-courses from cadata

if exist('parula')>0
  cmap=parula(length(peakstats));
else
  cmap=paruly(length(peakstats));
end

figs.frac_unstable_peaktime=figure();
frac_peaktime=zeros(length(peakstats),ndays+1);

for i=1:length(peakstats)

  peakstats(i).peak_stable=peakstats(i).peak_stable(:,all(peakstats(i).peak_ispeak));
  nrois=size(peakstats(i).peak_stable,2)
  unstable=zeros(ndays+1,nrois);

  for j=1:ndays
    unstable(j+1,:)=(peakstats(i).peak_stable(1,:)==1&peakstats(i).peak_stable(j+1,:)==0);
  end

  count=cumsum(unstable);
  n=sum(count>0,2);
  frac_peaktime(i,:)=n./nrois;
  plot([0:4],frac_peaktime(i,:),'ko-','color',cmap(i,:),'markersize',8,'markerfacecolor',[1 1 1]);
  hold on;

end

figs_stats.frac_peaktime=frac_peaktime;

ylim([0 1])
xlim([-.5 4.5])
set(gca,'TickLength',[0 0],'YTick',[0:.5:1],'XTick',[0:4],'FontSize',7)

if filewrite
    fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,'fig5_catimecourse.txt'),'w+');
    fprintf(fid,'Multi-unit stability vs calcium: p=%e z=%g\n',figs_stats.mu_v_ca_all.pval,figs_stats.mu_v_ca_all.zval);
    fprintf(fid,'Within day v between day variability: p=%e r=%g\n',figs_stats.drift.var_v_change.p,figs_stats.drift.var_v_change.r);
    fprintf(fid,'N(multi-unit): %i\nN(ROIS,1): %i\nN(ROIS,2): %i\nN(ROIS,3): %i',length(plotpoints{1}),length(plotpoints{2}),length(plotpoints{3}),length(plotpoints{4}));
    fclose(fid);
end
