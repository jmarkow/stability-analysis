function fig=stan_raster_multi(COLORS)
%
%

if nargin<1 | isempty(COLORS)
	COLORS='jet';
end

load('/Volumes/steel_age/intan_data/lpur95/templates/motif1_padding/template_data.mat')

[s,f,t]=zftftb_pretty_sonogram(TEMPLATE,25e3,'len',70,'overlap',69.5,'zeropad',0);

fig=figure();

t=t+.2;
ax(1)=subplot(4,2,1);imagesc(t,f/1e3,s);
axis xy
colormap(COLORS)
set(gca,'xtick',[],'ytick',[]);
ylim([0 10])

cat_times=[];
cat_trials=[];

load('/Volumes/bronze_age/figures/hvc_fields/data/spikefield_stats/int2_mu2/int2/lpur95/hvc/cell1/spikedata_ch13_cl1.mat')

cat_times=[cat_times cluster.times{1}];
cat_trials=[cat_trials cluster.trials{1}];

%
ax(2)=subplot(4,2,3);
spikoclust_raster(cluster.times{1}/cluster.parameters.fs,cluster.trials{1});
set(gca,'ytick',[],'xtick',[]);
ylim([0 200]);

load('/Volumes/bronze_age/figures/hvc_fields/data/spikefield_stats/int2_mu2/mu2/lpur95/hvc/cell3/spikedata_ch15_cl1.mat')

cat_times=[cat_times cluster.times{1}];
cat_trials=[cat_trials cluster.trials{1}];

ax(3)=subplot(4,2,5);
spikoclust_raster(cluster.times{1}/cluster.parameters.fs,cluster.trials{1});
set(gca,'ytick',[],'xtick',[]);
ylim([0 200]);

load('/Volumes/bronze_age/figures/hvc_fields/data/spikefield_stats/int2_mu2/mu2/lpur95/hvc/cell1/spikedata_ch1_cl1.mat')

cat_times=[cat_times cluster.times{1}];
cat_trials=[cat_trials cluster.trials{1}];

ax(4)=subplot(4,2,7);
spikoclust_raster(cluster.times{1}/cluster.parameters.fs,cluster.trials{1});
set(gca,'ytick',[]);
ylim([0 200]);

%ax(1)=subplot(4,1,1);imagesc(t,f/1e3,s);
%axis xy
%colormap(COLORS)
%set(gca,'xtick',[],'ytick',[]);
%ylim([0 10])

ax(5)=subplot(4,2,[4:2:8]);
spikoclust_raster(cat_times/cluster.parameters.fs,cat_trials);
%ax(5)=get(fig.combined,'CurrentAxes');
ylim([0 200]);
linkaxes(ax,'x');
set(gca,'ytick',[]);


