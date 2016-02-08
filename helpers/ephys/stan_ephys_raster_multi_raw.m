function fig=stan_raster_multi(COLORS)
%
%

[options,dirs]=stan_preflight;

channels=[1 6 9 13 15];

if nargin<1 | isempty(COLORS)
	COLORS='jet';
end

load(fullfile(dirs.agg_dir,dirs.template_dir,'lpur95_motif1_padding.mat'));
load(fullfile(dirs.agg_dir,'multi_example','lpur95_2012-07-16_motif1padding_ephys1.mat'));

nplots=length(channels)+1;
[s,f,t]=zftftb_pretty_sonogram(template.data,template.fs,'len',70,'overlap',69.5,'zeropad',0);

fig=figure();

t=t+.2;
ax(1)=subplot(nplots,2,1);imagesc(t,f/1e3,s);
axis xy
colormap(COLORS)
set(gca,'xtick',[],'ytick',[]);
ylim([0 10])

idx=zeros(1,length(channels));

for i=1:length(channels)
	idx(i)=find(CHANNELS==channels(i));
end

counter=3;

[b,a]=ellip(3,.2,40,[600 4e3]/(25e3/2),'bandpass');

for i=1:length(channels)
	spike_data=filtfilt(b,a,double(EPHYS_DATA(:,:,idx(i))));
	spike_threshold=options.sigma_t*median(abs(spike_data)/.6745);
	spikes=spikoclust_spike_detect_mu(spike_data,spike_threshold,25e3,'visualize','n','method','b');
	ax(i+1)=subplot(nplots,2,counter);
	spikoclust_raster(spikes.times,spikes.trial,'fs',25e3);
	set(gca,'YTick',[],'FontSize',7);
	title(['Channel ' num2str(i)]);
	counter=counter+2;
end

spike_data=filtfilt(b,a,sum(double(EPHYS_DATA(:,:,idx)),3));
spike_threshold=options.sigma_t*median(abs(spike_data)/.6745);
spikes=spikoclust_spike_detect_mu(spike_data,spike_threshold,25e3,'visualize','n','method','b');
ax(end+1)=subplot(nplots,2,4:2:nplots*2);
spikoclust_raster(spikes.times,spikes.trial,'fs',25e3);
set(gca,'FontSize',7);
title('Combined')
linkaxes(ax,'x');
linkaxes(ax(2:end),'xy');
ylim([0 200]);
set(gca,'ytick',[]);


