function [FIG_SWARM,FIG_STATS,plotpoints]=stan_fig_mu_ca_compare(CADATA,MUDATA)
%
%

swarm_offset=1;
idx=MUDATA.days_since<5;

bird=MUDATA.birdid(idx);
val=MUDATA.val_mu(idx);
[uniq_bird]=unique(bird);

plotpoints{1}=[];
for i=1:length(uniq_bird)
  plotpoints{1}(end+1)=mean(val(bird==uniq_bird(i)));
end

FIG_STATS.pval=nan(1,length(CADATA));

for i=1:length(CADATA)
  tot_lags=length(CADATA(i).rmat_mu.lag.all);
  tmp=cat(1,CADATA(i).rmat_mu.lag.all{2:min(tot_lags,5)});
  if size(tmp,1)>1
    tmp=mean(tmp);
  end
  plotpoints{i+1}=tmp(:);
  FIG_STATS.pval(i)=ranksum(plotpoints{1},plotpoints{i+1},'tail','right')
end

[FIG_STATS.all.pval,~,tmp_stats]=ranksum(plotpoints{1},cat(1,plotpoints{2:end}));
FIG_STATS.all.zval=tmp_stats.zval;

npre=1;
pos=[1 ones(1,npre)+swarm_offset ones(1,length(CADATA)-npre)];
pos=cumsum(pos);

if exist('parula')>0
    cmap=parula(length(plotpoints)-npre);
else
    cmap=paruly(length(plotpoints)-npre);
end

swarm_colors=[repmat([.7 .7 .7],[npre 1]);cmap];

% pairwise ranksum, Holm-Bonferonni stepdown

FIG_SWARM=figure();
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

FIG_HIST=figure();
bins=[0:.05:1];
n=ksdensity(plotpoints{1},bins);
n2=ksdensity(cat(1,plotpoints{2:end}),bins);

plot(bins,cumsum(n)./sum(n));
hold on;
plot(bins,cumsum(n2)./sum(n2));

figure();
markolab_stairplot(n2./sum(n2),bins,'facecolor','r','method','a');
hold on;
markolab_stairplot(n./sum(n),bins,'facecolor','b','method','a');

%
% figs.box=figure();
% markolab_boxplot(plotpoints,[],'box_color',swarm_colors,'feature_idx',[4:-1:1],'med_color',repmat([0 0 0],[4 1]));
% ylim([.2 1])
% set(gca,'YTick',get(gca,'YLim'),'outerPosition',[.05 0 1 1])
% yh=ylabel('Correlation (R)')
% set(yh,'position',get(yh,'position')+[.16 0 0],'FontSize',7);
