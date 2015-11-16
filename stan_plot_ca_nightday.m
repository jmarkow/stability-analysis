function [figs figs_stats]=stan_plot_ca_nightdat(CASTATS)
%
%
%

% stitch everything together appropriately, alternating first/last

nlags=4;
npoints=nlags*2; % nlags morning/evening
nboots=1e3;

agg=cell(1,npoints);
counter=1;
for i=1:nlags
  tmp=[];
  tmp2=[];
  for j=1:length(CASTATS)
    tmp=[tmp;CASTATS(j).rmat_mu.lag.day{i}(:)];
    tmp2=[tmp2;CASTATS(j).rmat_mu.lag.night{i}(:)];
  end

  agg{counter}=tmp;
  agg{counter+1}=tmp2;
  counter=counter+2;

end

% get mean

mu=cellfun(@mean,agg);

% get confidence interfval

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
  % block off each day
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

% figs.daycompare=figure();
% plotSpread(days([1 length(days)]),'binwidth',.0075,'distributionColors',{'c','c'});
% hold on;
% markolab_boxplot(days([1 length(days)]),[],'box_color','none','med_color',[0 0 0])
% set(gca,'XTick',[1:2],'XTickLabel',[1 4],'TickLength',[0 0]);
% ylim([-1 1.1])
% set(gca,'YTick',[-1 0 1]);
% ylabel('Correlation (R)');
% xlabel('Lag (days)');

% figs.nightcompare=figure();
% plotSpread(nights([1 length(nights)]),'binwidth',.0075,'distributionColors',{'b','b'});
% hold on;
% markolab_boxplot(nights([1 length(nights)]),[],'box_color','none','med_color',[0 0 0]);
% set(gca,'XTick',[1:2],'XTickLabel',[1 4],'TickLength',[0 0]);
% ylim([-1 1.1])
% set(gca,'YTick',[-1 0 1]);
% ylabel('Correlation (R)');
% xlabel('Lag (days)');
%
% figs.firstlast=figure();
% e{1}=nights{1};
% e{2}=days{end};
% plotSpread(e,'binwidth',.0075,'distributionColors',{'b','c'});
% hold on;
% markolab_boxplot(e,[],'box_color','none','med_color',[0 0 0]);
% set(gca,'XTick',[1:2],'XTickLabel',[1 4],'TickLength',[0 0]);
% ylim([-1 1.1])
% set(gca,'YTick',[-1 0 1]);
% ylabel('Correlation (R)');
% xlabel('Lag (days)');

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
  figs_stats.pairtest(i)=signrank(norm_day{i},norm_night{i});

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

% for i=1:length(days)
%   days{i}=atanh(days{i});
%   nights{i}=atanh(nights{i});
%   tmp_mu=mean(days{i});
%   tmp_std=std(days{i});
%   days{i}=(days{i}-tmp_mu)./tmp_std;
%   nights{i}=(nights{i}-tmp_mu)./tmp_std;
% end
%
% % normalize, is day-night change significant when pooled?
%
% nights=cat(1,nights{:});
% days=cat(1,days{:});
% e{1}=days;
% e{2}=nights;
%
% figs.daynightz=figure();
% plotSpread(e,'binwidth',.0075,'distributionColors',{'c','b'});
% hold on;
% markolab_boxplot(e,[],'box_color','none','med_color',[0 0 0]);
% set(gca,'TickLength',[0 0])
% set(gca,'YTick',[-4 0 4]);
% set(gca,'XTick',[1:2],'XTickLabel',{'Day','Night'});
% ylabel('Night-day corr. (Z)')
