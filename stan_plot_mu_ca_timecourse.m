function figs=stan_plot_barecarbon_ca_timecourse()
%
%
%

% upsample?
interp_factor=5;
ndays=4;
swarm_offset=1;
nboots=10e3;

[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats.mat'),'stats');
load(fullfile(dirs.agg_dir,dirs.datastore_dir,'mu_baseline_stability.mat'),'teststats');

% get baseline from ephys data

% plot time-courses from cadata


figs.catime=figure();
colors=paruly(3);
x=0:ndays;

baseline_ci=bootci(nboots,{@mean,teststats.val_mu(teststats.days_since<5)},'type','cper');
markolab_shadeplot(x,repmat([1;baseline_ci(2)],[1 length(x)]),[.7 .7 .7],'k');
hold on;

for i=1:length(stats)
  mu_ci=bootci(nboots,{@mean,stats(i).rmat_mu'},'type','cper');
  mu=mean(stats(i).rmat_mu');
  xx=0:1/interp_factor:ndays;

  mu_interp=interp1(x(:),mu(:),xx(:),'spline');
  mu_ci_interp=interp1(x(:),mu_ci',xx(:),'spline');

  markolab_shadeplot(xx,mu_ci_interp',colors(i,:),'k');
  plot(xx,mu_interp,'k--');
  plot(x,mu,'ko','markerfacecolor','w');
end

yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0]);
ylim([.6 1]);
ylimits=ylim();
set(gca,'YTick',ylimits,'TickLength',[0 0],'FontSize',9);

% plotspread version??

figs.cabeeswarm=figure();

% collect points for each day

allcapoints=[];
plotpoints{1}=teststats.val_mu(teststats.days_since<5);
pval=ones(1,length(stats))*NaN;

for i=1:length(stats)
  tmp=stats(i).rmat_mu(3:5,:);
  plotpoints{i+1}=tmp(:);
  [pval(i)]=ranksum(plotpoints{1},plotpoints{i+1},'tail','right');
end

% for statistical comparison?

pos=[1 1+swarm_offset ones(1,length(stats)-1)];
pos=cumsum(pos)
cmap=paruly(length(plotpoints)-1);
swarm_colors=[.7 .7 .7;cmap];

% pairwise ranksum, Holm-Bonferonni stepdown

pval_cor=markolab_bonf_holm(pval,.05)
plotSpread(plotpoints,'xValues',pos,'binWidth',.15,'distributionColors',swarm_colors);
ylim([-1 1]);
ylimits=ylim();
set(gca,'Xtick',pos([1 3]),'XTickLabel',{'Multi-unit','PNs'},'TickLength',[0 0],'YTick',[ylimits(1) ylimits(2)],...
  'FontSize',9);
yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0]);
