function [figs,fig_stats]=stan_plot_barecarbon_ca_timecourse()
%
%
%

% upsample?
interp_factor=5;
ndays=4;
swarm_offset=1;
nboots=1e3;

[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats.mat'),'stats');
load(fullfile(dirs.agg_dir,dirs.datastore_dir,'mu_baseline_stability.mat'),'teststats');

% get baseline from ephys data

% plot time-courses from cadata


figs.catime=figure();
colors=paruly(3);
x=0:ndays;

%baseline_ci=bootci(nboots,{@mean,teststats.val_mu(teststats.days_since<5)},'type','cper');
%markolab_shadeplot(x,repmat([1;baseline_ci(2)],[1 length(x)]),[.7 .7 .7],'k');
%hold on;

for i=1:length(stats)

  mu=ones(1,ndays+1);
  mu_ci=ones(2,ndays+1);
  cur_data=stats(i).rmat_mu.lag.all;

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
ylim([.6 1]);
ylimits=ylim();
set(gca,'YTick',ylimits,'TickLength',[0 0],'FontSize',9);

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

%plotpoints{1}=teststats.val_mu(teststats.days_since<5);
pval=ones(1,length(stats))*NaN;
pval2=ones(1,length(stats))*NaN;

% plotpoints{2}=[];
% inc_rois={};
% for i=1:length(stats)
%   tmp=(stats(i).rmat_mu_withinday(1,:));
%   %inc_rois{i}=find(tmp>.6);
%   %tmp(tmp<.6)=[];
%   plotpoints{2}=[plotpoints{2};tmp(:)];
% end

for i=1:length(stats)
  %tmp=stats(i).rmat_mu(2:5,inc_rois{i});
  %tmp=stats(i).rmat_mu.acrossday(2:5,:);
  tmp=cat(1,stats(i).rmat_mu.lag.all{2:5});
  tmp=mean(tmp);
  plotpoints{i+1}=tmp(:);
  pval(i)=ranksum(plotpoints{1},plotpoints{i+1},'tail','right');
  %compare=stats(i).rmat_mu_withinday;
  %pval2(i)=ranksum(compare(:),plotpoints{i+1},'tail','right');
end

% compare within a day
% for statistical comparison?

npre=1;
pos=[1 ones(1,npre)+swarm_offset ones(1,length(stats)-npre)];
pos=cumsum(pos);
cmap=paruly(length(plotpoints)-npre);
swarm_colors=[repmat([.7 .7 .7],[npre 1]);cmap];

% pairwise ranksum, Holm-Bonferonni stepdown

figs.beeswarm=figure();
pval_cor=markolab_bonf_holm(pval,.05)
h=plotSpread(plotpoints,'xValues',pos,'binWidth',.4,'distributionColors',swarm_colors);
for i=1:length(h)-1
  set(h{i},'markersize',8);
end
ylim([0 1]);
ylimits=ylim();
set(gca,'Xtick',pos([1 3]),'XTickLabel',{'Multi-unit','PNs'},'TickLength',[0 0],'YTick',[ylimits(1) ylimits(2)],...
  'FontSize',9,'outerPosition',[.05 0 1 1])
yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0]);


swarm_colors

figs.box=figure();
markolab_boxplot(plotpoints,[],'box_color',swarm_colors,'feature_idx',[4:-1:1],'med_color',repmat([0 0 0],[4 1]));
ylim([.2 1])
set(gca,'YTick',get(gca,'YLim'),'outerPosition',[.05 0 1 1])
yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0]);
