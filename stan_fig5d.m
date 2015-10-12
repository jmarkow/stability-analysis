%
%
%
%

[options,dirs]=stan_preflight;
%load(fullfile(dirs.agg_dir,dirs.fig_dir,'AwesomeNess.mat'));
load(fullfile(dirs.agg_dir,dirs.fig_dir,'combined_rois.mat'));

newdata=stan_cadata_format(data1,data2,data3,data4,data5);
[corrvals,comparevals,pmat,zmat]=stan_cadata_drift_analyze(newdata,'smoothing',0.05,'padding',.6);

test_alpha=.05;
[ndays,nrois]=size(zmat);

tmp=pmat(2:end,:);
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

xdata=0:4;
newxdata=0:1/5:4;

ci_stable_interp=interp1(xdata(:),ci_stable',newxdata(:),'spline');
ci_unstable_interp=interp1(xdata(:),ci_unstable',newxdata(:),'spline');
mu_stable_interp=interp1(xdata(:),mu_stable',newxdata(:),'spline');
mu_unstable_interp=interp1(xdata(:),mu_unstable',newxdata(:),'spline');


fig.drift_plot=figure();
h(1)=markolab_shadeplot(newxdata,ci_stable_interp',[0 0 1],'none');
hold on;
h(2)=markolab_shadeplot(newxdata,ci_unstable_interp',[1 0 0],'none');
plot(newxdata,mu_stable_interp,'b-','markerfacecolor','b');
plot(newxdata,mu_unstable_interp,'r-','markerfacecolor','r');
plot(xdata,mu_stable,'bo','markerfacecolor','b');
plot(xdata,mu_unstable,'ro','markerfacecolor','r');
alpha(.25);
L=legend(h,{'Unstable','Stable'});
legend boxoff;
set(L,'location','SouthWest');
set(fig.drift_plot,'position',[500 500 400 350],'PaperPositionMode','auto');
xlim([-.1 4]);
set(gca,'Ticklength',[0 0],'XTick',0:4,'FontSize',12);
ylim([-2 .5]);
ylabel('Correlation (rel to baseline, Z)');
xlabel('Days');
box on;

fig.per_bars=figure();
plot([1:4],nunstable./nrois,'ro-','markerfacecolor','r');
ylim([.3 .5]);
set(fig.per_bars,'position',[500 500 400 350],'PaperPositionMode','auto');
xlim([.9 4.1])
set(gca,'TickLength',[0 0],'xtick',[1:4],'ytick',[.3:.1:.5],'FontSize',12)
ylabel('Percent cells unstable');
xlabel('Days');

names=fieldnames(fig);

for i=1:length(names)
	%set(fig.(names{i}),'units','centimeters','position',[10 10 16 5],'paperpositionmode','auto');
	markolab_multi_fig_save(fig.(names{i}),fullfile(dirs.agg_dir,dirs.fig_dir),[ 'driftanalysis_' names{i} ],'eps,png,fig',...
		'renderer','painters');
end
