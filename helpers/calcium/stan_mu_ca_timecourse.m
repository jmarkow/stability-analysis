function [figs,figs_stats]=stan_plot_barecarbon_ca_timecourse()
%
%
%

% upsample?
interp_factor=5;
ndays=4;
nboots=1e2;
filewrite=false;

[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats_new.mat'),'stats');
peakstats=load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats_peaktime.mat'),'stats');
peakstats=peakstats.stats;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'mu_baseline_stability.mat'),'teststats');

% get baseline from ephys data

% lag correlation time-course

if exist('parula')>0
    colors=parula(length(stats));
else
    colors=paruly(length(stats));
end

figs.catimecourse=stan_plot_cacorr_timecourse(stats,colors,nboots);

% plotspread version??

% collect points for each day

figs.beeswarm=stan_plot_mu_ca_compare(stats,teststats);

% monte carlo permutation test for fraction of unstable cells

figs_stats.drift.pval=cell(1,length(stats));

for i=1:length(stats)

  % get pvals for each lag

  nrois=size(stats(i).rmat_mu.lag.all{1},2);

  cur_data=stats(i).rmat_mu.lag.all;
  boot_data=stats(i).rmat_mu.bootstrap.lag.all;
  x=find(cellfun(@length,cur_data)>0);
  figs_stats.drift.pval{i}=cell(1,length(x));

  for j=1:length(x)
    cur_datum=stats(i).rmat_mu.lag.all{x(j)};
    boot_datum=stats(i).rmat_mu.bootstrap.lag.all{x(j)};
    boot_datum=mean(cat(3,boot_datum{:}),3);

    if size(cur_datum,1)>1
      cur_datum=mean(cur_datum);
    end

    cur_datum=repmat(cur_datum,[size(boot_datum,1) 1]);
    pval=mean(cur_datum>boot_datum)+1/size(boot_datum,1);
    figs_stats.drift.pval{i}{x(j)}=pval;

  end

end

% test for whether within day variability accounts for change
variability=[];
change=[];

for i=1:length(stats)
  x=find(cellfun(@length,stats(i).rmat_mu.lag.day)>0);
  variability=[variability mean(stats(i).rmat_mu.lag.day{x(1)})];

  if size(stats(i).rmat_mu.lag.day{x(2)},1)>0
    comparison=mean(stats(i).rmat_mu.lag.day{x(2)});
  else
    comparison=stats(i).rmat_mu.lag.day{x(2)};
  end
  change=[change mean(stats(i).rmat_mu.lag.day{x(1)})-comparison];
end

figs.var_v_change=figure();
scatter(variability,change)
[r3,p3]=corr(variability(:),change(:),'type','pearson');
figs_stats.drift.var_v_change.p=p3;
figs_stats.drift.var_v_change.r=r3;

% plot cum fraction of unstable cells for each bird (MC correction for number of ROIs)

figs.frac_unstable=figure();
%frac=zeros(length(stats),ndays+1);
frac=cell(1,length(stats));
figs_stats.surv_time=cell(1,length(stats));

for i=1:length(stats)

  nrois=size(stats(i).rmat_mu.lag.day{1},2);

  x=find(cellfun(@length,stats(i).rmat_mu.lag.day)>0);
  unstable=nan(length(x),nrois);

  for j=1:length(x)
    unstable(j,:)=markolab_bonf_holm(figs_stats.drift.pval{i}{x(j)},.05)<0.01;
    %unstable(j,:)=figs_stats.drift.pval{i}{x(j)}<=.001;
  end

  tmp=[];
  for j=1:size(unstable,2)
    tmp2=x(min(find(unstable(:,j))));
    if ~isempty(tmp2)
      tmp=[tmp tmp2];
    end
  end

  count=cumsum(unstable);
  n=sum(count>0,2);

  frac{i}=n./nrois;

  figs_stats.surv_time{i}=tmp;
  figs_stats.drift.unstable_n(i)=size(stats(i).rmat_mu.lag.all{1},2);
  plot(x-1,frac{i},'ko-','color',colors(i,:),'markersize',8,'markerfacecolor',[1 1 1]);
  hold on;

  frac{i}=unstable;
end

figs_stats.drift.unstable=frac;
ylim([0 1])
xlim([-.5 4.5])
set(gca,'TickLength',[0 0],'YTick',[0:.5:1],'XTick',[0:4],'FontSize',7)

figs.frac_unstable_peaktime=figure();
frac_peaktime=zeros(length(peakstats),ndays+1);

for i=1:length(peakstats)
  nrois=size(peakstats(i).peak_stable,2);
  unstable=zeros(ndays+1,nrois);
  for j=1:ndays
    unstable(j+1,:)=(peakstats(i).peak_stable(1,:)==1&peakstats(i).peak_stable(j+1,:)==0);
  end

  count=cumsum(unstable);
  n=sum(count>0,2);
  frac_peaktime(i,:)=n./nrois;
  plot([0:4],frac_peaktime(i,:),'ko-','color',colors(i,:),'markersize',8,'markerfacecolor',[1 1 1]);
  hold on;
end

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
