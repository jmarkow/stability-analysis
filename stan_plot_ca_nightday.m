function [figs figs_stats]=stan_plot_ca_nightdat(CASTATS)
%
%
%

% stitch everything together appropriately, alternating first/last

nlags=length(CASTATS(1).rmat_mu.lag.day)-1; % first point is lag 0
npoints=nlags*2; % nlags morning/evening
nboots=1e3;

agg=cell(1,npoints);
counter=1;

for i=1:nlags

  tmp=[];
  tmp2=[];

  for j=1:length(CASTATS)
    tmp=[tmp;CASTATS(j).rmat_mu.lag.day{i+1}(:)];
    tmp2=[tmp2;CASTATS(j).rmat_mu.lag.night{i+1}(:)];
  end

  agg{counter}=tmp;
  agg{counter+1}=tmp2;
  counter=counter+2;

end

mu=cellfun(@mean,agg);
mu_ci=zeros(2,length(agg));

for i=1:length(agg)
  mu_ci(:,i)=bootci(nboots,{@mean,agg{i}},'type','cper');
end

x=[1:length(agg)];

figs.daynightcompare=figure();

counter=1;
for i=1:nlags
  markolab_shadeplot([x(counter)-.5 x(counter)+.5],[0 0;1 1],'c','none');
  hold on;
  markolab_shadeplot([x(counter+1)-.5 x(counter+1)+.5],[0 0;1 1],'b','none');
  hold on;
  plot([x(counter+1)+.5 x(counter+1)+.5],[0 1],'k-');
  counter=counter+2;
end

hold on;
plot(x,mu,'k-');
hold on;
stan_plot_dot_error(x,mu,mu_ci,repmat([1 2],[1 4]),'colors',[0 0 0;0 0 0]);
hold on;
xlim([x(1)-.5 x(end)+.5]);
ylim([.45 .85]);
ylabel('Correlation (R)');
xlabel('Lag (day)');
set(gca,'XTick',x(1)+.5:2:x(end),'XTickLabel',[1:4],'Ytick',get(gca,'ylim'));
set(gca,'layer','top');

% group together night and day (histogram of day-night????)

days=agg(1:2:npoints);
nights=agg(2:2:npoints);

plot_mu=zeros(1,length(days)*2);
plot_mu_ci=zeros(2,length(days)*2);
fig_stats.pairtest=nan(1,length(days));
counter=1;

% form a paired plot

for i=1:length(days)

  mu_day=mean(days{i});
  std_day=std(days{i});

  norm_day{i}=(days{i}-mu_day)./(std_day);
  norm_night{i}=(nights{i}-mu_day)./(std_day);

  plot_mu(counter)=mean(norm_day{i});
  plot_mu(counter+1)=mean(norm_night{i});

  plot_mu_ci(:,counter)=bootci(nboots,{@mean,norm_day{i}},'type','cper');
  plot_mu_ci(:,counter+1)=bootci(nboots,{@mean,norm_night{i}},'type','cper');

  figs_stats.lagtest(i)=ranksum(norm_night{i},norm_night{1});
  figs_stats.pairtest(i)=signrank(days{i},nights{i});
  counter=counter+2;

end

cmap=winter(4);
offset=.05;
figs.testing=figure();
counter=1;
x=[];
for i=1:2:length(days)*2
  x=[x 1+offset*counter 2+offset*counter];
  h(counter)=plot([1+offset*counter 2+offset*counter],plot_mu([i i+1]),'k-','color',cmap(counter,:));
  hold on;
  counter=counter+1;
end

stan_plot_dot_error(x,plot_mu,plot_mu_ci,[1 1 2 2 3 3 4 4],'colors',cmap);
L=legend(h,{'1 day','2 days','3 days','4 days'})
legend boxoff;
set(L,'location','northwest','FontSize',6)
xlim([.5 2.5])
set(gca,'XTick',[1 2],'XTickLabel',{'First','Last'})

% for each animal, check shift within each ROI (between relative to within)

pool1=[];
pool2=[];
pool3=[];
pool4=[];

dprime=[];

for i=1:length(CASTATS)

  mu=mean(CASTATS(i).rmat_mu.lag.day{1});
  %idx=mu>prctile(mu,75);
  idx=mu>0

  tmp1=mean(CASTATS(i).rmat_mu.lag.day{4}(:,idx));
  tmp2=mean(CASTATS(i).rmat_mu.lag.day{1}(:,idx));
  tmp3=std(CASTATS(i).rmat_mu.lag.day{1}(:,idx));
  tmp4=std(CASTATS(i).rmat_mu.lag.day{2}(:,idx));
  nrois=length(tmp1);

  %for j=1:nrois
  %  pval(j)=ranksum(CASTATS(i).rmat_mu.lag.day{1}(:,j),CASTATS(i).rmat_mu.lag.day{2}(:,j),'tail','right');
  %end

  pool1=[pool1 tmp1];
  pool2=[pool2 tmp2];
  pool3=[pool3 tmp3];
  pool4=[pool4 tmp4];
  %pool4=[pool4 pval];

  %dprime=[dprime tmp3];
  %dprime=[dprime (pool1-pool2)./pool2];

end

d{1}=pool1;
d{2}=pool2;
figure();markolab_boxplot(d);
pause();
e{1}=(pool1-pool2)./sqrt(pool3.*pool3);
%figure();hist(e{1},20)

figs.overnight_swarm=figure();
plotSpread(e,'binWidth',.01);
ylim([-5 5]);
set(gca,'YTick',[-5 0 5],'TickLength',[0 0],'XTick',[],'outerposition',[.05 0 1 1]);
hold on
plot(get(gca,'xlim'),[0 0],'k-');
yh=ylabel('\DeltaCorr. (Z)');
%set(yh,'position',get(yh,'position')+[.0005 0 0])

figs.overnight_box=figure();
markolab_boxplot(e);
ylim([-3 3]);
set(gca,'YTick',[-3 0 3],'TickLength',[0 0],'XTick',[],'outerposition',[.05 0 1 1]);
plot(get(gca,'xlim'),[0 0],'k-');
ylabel('\DeltaCorr. (Z)');

%bootci(1e3,{@median,(pool1-pool2)},'type','cper')
