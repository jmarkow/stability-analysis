
[~,dirs]=stan_preflight;

ca.con=load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-con.mat']),'stats');
ca.lib=load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-lib.mat']),'stats');
load(fullfile(dirs.agg_dir,dirs.datastore_dir,'mu_baseline_stability.mat'),'teststats');

%%


idx=teststats.days_since<5;
bird=teststats.birdid(idx);
val=teststats.val_mu(idx);
uniq_bird=unique(bird);

mupoints=[];
for i=1:length(uniq_bird)
  mupoints(i)=mean(val(bird==uniq_bird(i)));
end

catypes={'con','lib'};

capoints=[];
for i=1:length(catypes)
  for j=1:length(ca.(catypes{i}).stats)
    tot_lags=length(ca.(catypes{i}).stats(j).rmat_mu.lag.all);
    tmp=cat(1,ca.(catypes{i}).stats(j).rmat_mu.lag.all{2:min(tot_lags,5)});
    if size(tmp,1)>1
        tmp=mean(tmp);
    end
    capoints.(catypes{i}){j}=tmp(:);
  end
end


FIG_HIST=figure();
bins=[0:.005:1];
n=ksdensity(mupoints,bins);
n2=ksdensity(cat(1,capoints.con{2:end}),bins);
n3=ksdensity(cat(1,capoints.lib{2:end}),bins);
h=[];
h(1)=plot(bins,cumsum(n)./sum(n),'k-','color',[.7 .7 .7]);
hold on;
h(2)=plot(bins,cumsum(n2)./sum(n2),'k-','color',[1 0 0]);
h(3)=plot(bins,cumsum(n3)./sum(n3),'k-','color',[0 0 1]);
box off;

L=legend(h,{'MU','PNs (all)','PNs (any)'});
legend boxoff;
set(L,'location','NorthWest','FontSize',7);
set(gca,'FontSize',7,'layer','top','ticklength',[0 0])
set(FIG_HIST,'units','centimeters','position',[10 10 4.5 2],'paperpositionmode','auto');
markolab_multi_fig_save(FIG_HIST,fullfile(dirs.agg_dir,dirs.fig_dir),'driftanalysis_revision-cdfs','eps,png,fig',...
 		'renderer','painters');