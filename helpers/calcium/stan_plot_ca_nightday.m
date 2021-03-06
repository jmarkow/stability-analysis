function [figs figs_stats]=stan_plot_ca_nightdat(EXT)
%
%
%

% stitch everything together appropriately, alternating first/last

if nargin<1
  EXT='con';
end

[options,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-' EXT '.mat']),'stats');

%nlags=length(stats(1).rmat_mu.lag.day)-1; % first point is lag 0
%npoints=nlags*2; % nlags morning/evening
nlags=4;
npoints=nlags*2;
nboots=1e3;
filewrite=true;
agg={};
counter=1;

for i=1:nlags

  % TODO: normalize

  tmp=[];
  tmp2=[];

  for j=1:length(stats)
    if length(stats(j).rmat_mu.lag.all)>=i+1
      if ~isempty(stats(j).rmat_mu.lag.all{i+1})
        for k=1:size(stats(j).rmat_mu.lag.day{i+1},1)
            tmp=[tmp stats(j).rmat_mu.lag.day{i+1}(k,:)-mean(stats(j).rmat_mu.bootstrap.lag.day{1}{1})];
        end

        %tmp=[tmp;stats(j).rmat_mu.lag.day{i+1}(:)];

        for k=1:size(stats(j).rmat_mu.lag.night{i+1},1)
            tmp2=[tmp2 stats(j).rmat_mu.lag.night{i+1}(k,:)-mean(stats(j).rmat_mu.bootstrap.lag.day{1}{1})];
        end


        %tmp2=[tmp2;stats(j).rmat_mu.lag.night{i+1}(:)];
      end
    end
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

% figs.daynightcompare=figure();
% 
% counter=1;
% for i=1:nlags
%   markolab_shadeplot([x(counter)-.5 x(counter)+.5],[0 0;1 1],'c','none');
%   hold on;
%   markolab_shadeplot([x(counter+1)-.5 x(counter+1)+.5],[0 0;1 1],'b','none');
%   hold on;
%   plot([x(counter+1)+.5 x(counter+1)+.5],[0 1],'k-');
%   counter=counter+2;
% end
% 
% hold on;
% plot(x,mu,'k-');
% hold on;
% stan_plot_dot_error(x,mu,mu_ci,repmat([1 2],[1 4]),'colors',[0 0 0;0 0 0]);
% hold on;
% xlim([x(1)-.5 x(end)+.5]);
% ylim([.45 .85]);
% ylabel('Correlation (R)');
% set(gca,'XTick',x(1)+.5:2:x(end),'XTickLabel',[1:4],'Ytick',get(gca,'ylim'));
% set(gca,'layer','top','FontSize',7);

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

  %figs_stats.circadian.lagtest(i)=ranksum(norm_night{i},norm_night{1});
  %figs_stats.circadian.pairtest(i)=signrank(days{i},nights{i});
  counter=counter+2;

end

cmap=winter(4);
offset=.05;
% figs.testing=figure();
% counter=1;
% x=[];
% for i=1:2:length(days)*2
%   x=[x 1+offset*counter 2+offset*counter];
%   h(counter)=plot([1+offset*counter 2+offset*counter],plot_mu([i i+1]),'k-','color',cmap(counter,:));
%   hold on;
%   counter=counter+1;
% end
% 
% stan_plot_dot_error(x,plot_mu,plot_mu_ci,[1 1 2 2 3 3 4 4],'colors',cmap);
% L=legend(h,{'1 day','2 days','3 days','4 days'})
% legend boxoff;
% set(L,'location','northwest','FontSize',6)
% xlim([.5 2.5])
% set(gca,'XTick',[1 2],'XTickLabel',{'First','Last'},'FontSize',7)

% for each animal, check shift within each ROI (between relative to within)

pool_pval1_left=[];
pool_pval1_right=[];
pool_pval2_left=[];
pool_pval2_right=[];

poolz1=[];
poolz2=[];
pool1=[];
pool2=[];

dprime=[];
%
% for i=1:length(stats)
%
%   poolz1=[];
%   poolz2=[];
%
%   pool1=[];
%   pool2=[];
%
%   bootmu=mean(cat(3,stats(i).rmat_mu.bootstrap.lag.day{1}{:}),3);
%   bootmu_mu=mean(bootmu);
%   bootmu_var=std(bootmu);
%
%   for j=1:size(stats(i).rmat_mu.lag.day{1},1)
%     tmp1=stats(i).rmat_mu.lag.day{1}(j,:);
%     pool1=[pool1;tmp1];
%     poolz1=[poolz1;(tmp1-bootmu_mu)./bootmu_var];
%   end
%
%   for j=1:size(stats(i).rmat_mu.lag.day{2},1)
%     tmp2=stats(i).rmat_mu.lag.day{2}(j,:);
%     pool2=[pool2;tmp2];
%     poolz2=[poolz2;(tmp2-bootmu_mu)./bootmu_var];
%   end
%
%   mean(poolz1)
%   mean(poolz2)
%
% end

for i=1:length(stats)

  if size(stats(i).rmat_mu.lag.day{2},1)>1
    mu2=mean(stats(i).rmat_mu.lag.day{2});
  else
    mu2=stats(i).rmat_mu.lag.day{2};
  end

  if size(stats(i).rmat_mu.lag.day{1},1)>1
    mu1=mean(stats(i).rmat_mu.lag.day{1});
  else
    mu1=stats(i).rmat_mu.lag.day{1};
  end

  tmp1=stats(i).rmat_mu.lag.day{1};
  tmp2=stats(i).rmat_mu.lag.day{2};

  pool1=[pool1 mu1];
  pool2=[pool2 mu2];

  bootmu=mean(cat(3,stats(i).rmat_mu.bootstrap.lag.day{1}{:}),3);

  pval1_right=mean(repmat(mu1,[size(bootmu,1) 1])>bootmu)+(1/size(bootmu,1));
  pval1_left=mean(repmat(mu1,[size(bootmu,1) 1])<bootmu)+(1/size(bootmu,1));

  pval2_right=mean(repmat(mu2,[size(bootmu,1) 1])>bootmu)+(1/size(bootmu,1));
  pval2_left=mean(repmat(mu2,[size(bootmu,1) 1])<bootmu)+(1/size(bootmu,1));

  poolz1=[poolz1 (mu1-mean(bootmu))./std(bootmu)];
  poolz2=[poolz2 (mu2-mean(bootmu))./std(bootmu)];

  pool_pval1_left=[pool_pval1_left pval1_left];
  pool_pval1_right=[pool_pval1_right pval1_right];

  pool_pval2_left=[pool_pval2_left pval2_left];
  pool_pval2_right=[pool_pval2_right pval2_right];

end

%significant pool

figs_stats.overnight.roi_pval1.left=[pool_pval1_left];
figs_stats.overnight.roi_pval1.right=[pool_pval1_right];
[figs_stats.overnight.all_pval1,~,tmpstats]=signrank(poolz1,0,'tail','left')
figs_stats.overnight.all_zval1=tmpstats.zval;

figs_stats.overnight.roi_pval2.left=[pool_pval2_left];
figs_stats.overnight.roi_pval2.right=[pool_pval2_right];
[figs_stats.overnight.all_pval2,~,tmpstats]=signrank(poolz2,poolz1,'tail','left')
figs_stats.overnight.all_zval2=tmpstats.zval;
figs_stats.overnight.all_pval2_raw=signrank(pool2,pool1,'tail','left');


bins=[-25:1:15];
est1=histc(poolz1,bins);
est2=histc(poolz1(pool_pval1_right<.05),bins);
est3=histc(poolz1(pool_pval1_left<.05),bins);

maxval=150;
%
% figs.roihist1=figure();
% ax(1)=markolab_stairplot(est1,bins,'facecolor',[.5 .5 .5],'edgecolor','k','method','p');
% hold on;
% ax(2)=markolab_stairplot(est2,bins,'facecolor',[0 0 1],'edgecolor','k','method','p');
% ax(3)=markolab_stairplot(est3,bins,'facecolor',[1 0 0],'edgecolor','k','method','p');
% plot([0 0],[0 maxval],'k--');
% xlim([-15 5]);
% set(gca,'YTick',[0:20:maxval],'xtick',[-15:5:5],'TickLength',[0 0],'FontSize',7,'layer','top');

est12=histc(poolz2,bins);
est22=histc(poolz2(pool_pval2_right<.05),bins);
est32=histc(poolz2(pool_pval2_left<.05),bins);

% figs.roihist1=figure();
% ax(1)=markolab_stairplot(est12,bins,'facecolor',[.5 .5 .5],'edgecolor','k','method','p');
% hold on;
% ax(2)=markolab_stairplot(est22,bins,'facecolor',[0 0 1],'edgecolor','k','method','p');
% ax(3)=markolab_stairplot(est32,bins,'facecolor',[1 0 0],'edgecolor','k','method','p');
% plot([0 0],[0 maxval],'k--');
% xlim([-15 5]);
% set(gca,'YTick',[0:20:maxval],'xtick',[-15:5:5],'TickLength',[0 0],'FontSize',7,'layer','top');

ax=[];
figs.comparehist=figure();
ax(1)=markolab_stairplot(est1,bins,'facecolor',[1 0 0],'edgecolor','k','method','p');
hold on;
ax(2)=markolab_stairplot(est12,bins,'facecolor',[0 0 1],'edgecolor','k','method','p');
plot([0 0],[0 maxval],'k--');
xlim([-15 5]);
set(gca,'YTick',[0:50:maxval],'xtick',[-15:5:10],'TickLength',[0 0],'FontSize',7,'layer','top');
L=legend(ax,{'Day','Night'});
legend boxoff;
set(L,'FontSize',5,'Location','Northwest')

%
% bins_raw=[-1:.01:1];
% est_raw=histc(pool2-pool1,bins_raw);
%
% figs.rawhist=figure();
% ax(1)=markolab_stairplot(est_raw,bins_raw,'facecolor',[.3 .3 .3],'edgecolor','k','method','p');
% hold on;
% plot([0 0],[0 maxval],'k--');
% xlim([-1 1]);

if filewrite
  fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,['fig6_caovernight-' EXT '.txt']),'w+');
  fprintf(fid,'Within day v overnight corr: p=%e z=%g',figs_stats.overnight.all_pval2,figs_stats.overnight.all_zval2);
  fclose(fid);
end

% figs.compareswarm=figure();
% d{1}=poolz1;
% d{2}=poolz2;
% plotSpread(d,'binwidth',.25)
% hold on;
% markolab_boxplot(d,[],'box_color','none','med_color',repmat([0 0 0],[4 1]))
