%
%
%
%

[options,dirs]=stan_preflight;
%load(fullfile(dirs.agg_dir,dirs.fig_dir,'AwesomeNess.mat'));
load(fullfile(dirs.agg_dir,dirs.fig_dir,'combined_rois.mat'));

newdata=stan_cadata_format(data1,data2,data3,data4,data5);
[corrvals,comparevals,pmat,zmat]=stan_cadata_drift_analyze(newdata,'smoothing',0.05,'padding',1);

test_alpha=.05;
[ndays,nrois]=size(zmat);

tmp=pmat([ 2 3 4 5],:);
tmp=tmp(:);
pmat_cor=markolab_bonf_holm(tmp,.05);
pmat_cor=reshape(pmat_cor',size(pmat(2:end,:)))
pmat_idx=pmat_cor<test_alpha;

nunstable=zeros(ndays-1,1);
nunstable(1)=sum(pmat_idx(1,:));

for i=2:ndays-1
	nunstable(i)=sum(any(pmat_idx(1:i,:)));
end


% shade plots

plot_idx=any(pmat_idx);

ci_stable=bootci(1e3,{@mean,zmat(:,~plot_idx)'},'type','cper');
ci_unstable=bootci(1e3,{@mean,zmat(:,plot_idx)'},'type','cper');
mu_stable=mean(zmat(:,~plot_idx)');
mu_unstable=mean(zmat(:,plot_idx)');

% spline smooth

xdata=1:5;
newxdata=1:1/10:5;

ci_stable_interp=interp1(xdata(:),ci_stable',newxdata(:),'spline');
ci_unstable_interp=interp1(xdata(:),ci_unstable',newxdata(:),'spline');
mu_stable_interp=interp1(xdata(:),mu_stable',newxdata(:),'spline');
mu_unstable_interp=interp1(xdata(:),mu_unstable',newxdata(:),'spline');


fig.drift_plot=figure();
h(1)=markolab_shadeplot(newxdata,ci_stable_interp',[0 0 1],'b');
hold on;
h(2)=markolab_shadeplot(newxdata,ci_unstable_interp',[1 0 0],'r');
plot(newxdata,mu_stable_interp,'b-','markerfacecolor','b');
plot(newxdata,mu_unstable_interp,'r-','markerfacecolor','r');
plot(xdata,mu_stable,'bo','markerfacecolor','b','markersize',6);
plot(xdata,mu_unstable,'ro','markerfacecolor','r','markersize',6);
alpha(.25);
L=legend(h,{'Stable','Unstable'});
legend boxoff;
set(L,'location','SouthWest');
set(fig.drift_plot,'position',[500 500 400 350],'PaperPositionMode','auto');
xlim([1 5.1]);
set(gca,'Ticklength',[0 0],'XTick',xdata,'FontSize',12);
ylim([-2 .5]);
ylabel('Correlation (rel to baseline, Z)');
xlabel('Days');
box off;
set(gca,'TickDir','out')
set(gca,'ticklength',[.025 .025])
set(gca,'FontSize',14)

fig.per_bars=figure();
plot([xdata(2:end)],(nunstable./nrois)*1e2,'ko-','markerfacecolor','k','markersize',15,'markerfacecolor','w');
ylim([30 50]);
set(fig.per_bars,'position',[500 500 400 350],'PaperPositionMode','auto');
xlim([1.9 5.1])
set(gca,'TickLength',[0 0],'xtick',[xdata(2:end)],'ytick',[30:10:50],'FontSize',14)
ylabel('Percent cells unstable');
xlabel('Days');
box off;

names=fieldnames(fig);

for i=1:length(names)
	%set(fig.(names{i}),'units','centimeters','position',[10 10 16 5],'paperpositionmode','auto');
    set(fig.(names{i}),'paperpositionmode','auto');
	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'driftanalysis_' names{i} ],'eps,png,fig',...
		'renderer','painters');
end
