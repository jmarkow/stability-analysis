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

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats_new.mat'),'stats');
peakstats=load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats_peaktime.mat'),'stats');
peakstats=peakstats.stats;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'mu_baseline_stability.mat'),'teststats');

% get baseline from ephys data

% plot time-courses from cadata


figs.catime=figure();
if exist('parula')>0
    colors=parula(3);
else
    colors=paruly(3);
end
x=0:ndays;

for i=1:length(stats)

  mu=ones(1,ndays+1);
  mu_ci=ones(2,ndays+1);
  cur_data=stats(i).rmat_mu.lag.all;

  % day one from bootstrap

  for j=1:ndays
    mu_ci(:,j+1)=bootci(nboots,{@mean,cur_data{j+1}(:)},'type','cper');
    mu(:,j+1)=mean(cur_data{j+1}(:));
  end

  xx=0:1/interp_factor:ndays;

  mu_interp=interp1(x(:),mu(:),xx(:),'spline');
  mu_ci_interp=interp1(x(:),mu_ci',xx(:),'spline');

  markolab_shadeplot(xx,mu_ci_interp',colors(i,:),'k');
  hold on;
  plot(xx,mu_interp,'k--');
  plot(x,mu,'ko','markerfacecolor','w');

end

yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0]);
ylim([.5 1]);
ylimits=ylim();
set(gca,'YTick',ylimits,'TickLength',[0 0],'FontSize',7);

% plotspread version??

% collect points for each day

allcapoints=[];
idx=teststats.days_since<5;

bird=teststats.birdid(idx);
val=teststats.val_mu(idx);
[uniq_bird]=unique(bird);

plotpoints{1}=[];
for i=1:length(uniq_bird)
  plotpoints{1}(end+1)=mean(val(bird==uniq_bird(i)));
end

figs_stats.mu_v_ca.pval=ones(1,length(stats))*NaN;

for i=1:length(stats)
  tmp=cat(1,stats(i).rmat_mu.lag.all{2:5});
  tmp=mean(tmp);
  plotpoints{i+1}=tmp(:);
  figs_stats.mu_v_ca.pval(i)=ranksum(plotpoints{1},plotpoints{i+1},'tail','right')
end

[figs_stats.mu_v_ca_all.pval,~,tmp_stats]=ranksum(plotpoints{1},cat(1,plotpoints{2:end}));
figs_stats.mu_v_ca_all.zval=tmp_stats.zval;
% compare within a day
% for statistical comparison?

npre=1;
pos=[1 ones(1,npre)+swarm_offset ones(1,length(stats)-npre)];
pos=cumsum(pos);

if exist('parula')>0
    cmap=parula(length(plotpoints)-npre);
else
    cmap=paruly(length(plotpoints)-npre);
end
swarm_colors=[repmat([.7 .7 .7],[npre 1]);cmap];

% pairwise ranksum, Holm-Bonferonni stepdown

length(plotpoints{1})
length(plotpoints{2})
length(plotpoints{3})
length(plotpoints{4})

figs.beeswarm=figure();
h=plotSpread(plotpoints,'xValues',pos,'binWidth',.6,'distributionColors',swarm_colors,'spreadFcn',{'xp',[]});
for i=1:length(h)-1
  set(h{i},'markersize',8);
end
ylim([0 1]);
ylimits=ylim();
set(gca,'Xtick',pos([1 3]),'XTickLabel',{'Multi-unit','PNs'},'TickLength',[0 0],'YTick',[ylimits(1) ylimits(2)],...
  'FontSize',7,'outerPosition',[.05 0 1 1])
yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0]);

figs.box=figure();
markolab_boxplot(plotpoints,[],'box_color',swarm_colors,'feature_idx',[4:-1:1],'med_color',repmat([0 0 0],[4 1]));
ylim([.2 1])
set(gca,'YTick',get(gca,'YLim'),'outerPosition',[.05 0 1 1])
yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0],'FontSize',7);

% monte carlo permutation test for fraction of unstable cells

figs_stats.drift.pval=cell(1,length(stats));
for i=1:length(stats)

  % get pvals for each lag
  nrois=size(stats(i).rmat_mu.lag.all{1},2);
  figs_stats.drift.pval{i}=zeros(ndays,nrois);

  for j=1:ndays
    cur_data=stats(i).rmat_mu.lag.all{j+1};
    boot_data=stats(i).rmat_mu.bootstrap.lag.all{j+1};

    boot_data=mean(cat(3,boot_data{:}),3);

    if size(cur_data,1)>1
      cur_data=mean(cur_data);
    end

    cur_data=repmat(cur_data,[size(boot_data,1) 1]);
    pval=mean(cur_data>boot_data)+1/size(boot_data,1);
    figs_stats.drift.pval{i}(j,:)=pval;
  end

end

% test for whether within day variability accounts for change
variability=[];
change=[];

for i=1:length(stats)
  variability=[variability mean(stats(i).rmat_mu.lag.day{1})];
  change=[change mean(stats(i).rmat_mu.lag.day{4})-mean(stats(i).rmat_mu.lag.day{1})];
end

figs.var_v_change=figure();
scatter(variability,change)
[r,p]=corr(variability(:),change(:),'type','spearman')
[r2,p2]=corrcoef(variability,change)
[r3,p3]=corr(variability(:),change(:),'type','pearson')
figs_stats.drift.var_v_change.p=p3;
figs_stats.drift.var_v_change.r=r3;

% plot cum fraction of unstable cells for each bird (MC correction for number of ROIs)

figs.frac_unstable=figure();
frac=zeros(length(stats),ndays+1);
figs_stats.surv_time=cell(1,length(stats));
for i=1:length(stats)
  nrois=size(stats(i).rmat_mu.lag.day{1},2);
  unstable=zeros(ndays+1,nrois);
  for j=1:ndays
    unstable(j+1,:)=(markolab_bonf_holm(figs_stats.drift.pval{i}(j,:),.05)<.01);
  end

  tmp=[];
  for j=1:size(unstable,2)
    tmp2=min(find(unstable(:,j)))
    if ~isempty(tmp2)
      tmp=[tmp tmp2];
    end
  end

  count=cumsum(unstable);
  n=sum(count>0,2);
  frac(i,:)=n./nrois;

  figs_stats.surv_time{i}=tmp;
  figs_stats.drift.unstable_n(i)=size(stats(i).rmat_mu.lag.all{1},2);
  plot([0:4],frac(i,:),'ko-','color',cmap(i,:),'markersize',8,'markerfacecolor',[1 1 1]);
  hold on;
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
  plot([0:4],frac_peaktime(i,:),'ko-','color',cmap(i,:),'markersize',8,'markerfacecolor',[1 1 1]);
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
