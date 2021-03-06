function ax=stan_songalign_raster(SPECT,SPIKES,varargin)
%
%
%
%
%


colors='jet';
nparams=length(varargin);
name='';
fs=[];
spike_width=1;
spike_height=.3;
plot_trials=[];
time_bar=.2;
columns=1;

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'colors'
			colors=varargin{i+1};
		case 'name'
			name=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'spike_width'
			spike_width=varargin{i+1};
		case 'spike_height'
			spike_height=varargin{i+1};
		case 'plot_trials'
			plot_trials=varargin{i+1};
		case 'time_bar'
			time_bar=varargin{i+1};
		case 'columns'
			columns=varargin{i+1};
	end
end


if length(SPIKES)>1
	nplots=2*length(SPIKES)+1;
else
	nplots=4;
end

ax(1)=subplot(nplots,columns,1);
imagesc(SPECT.t,SPECT.f/1e3,SPECT.s);
colormap(colors);
axis xy;
box off;
set(gca,'TickDir','out','TickLength',[0 0]);
ylim([0 10]);
set(gca,'YTick',[],'XTick',[]);
title([name]);
%ylabel('Fs (kHz)');


nplots
columns

plotcounter=columns+1;

if length(SPIKES)>1
	for i=1:length(SPIKES)

		idx=[plotcounter plotcounter+columns];
		plotcounter=plotcounter+columns*2;

		if isempty(fs)
			spike_fs=SPIKES(i).fs;
		else
			spike_fs=fs;
		end

		ax(i+1)=subplot(nplots,columns,idx);
		spikoclust_raster(SPIKES(i).times/spike_fs,SPIKES(i).trial,...
            'spike_width',spike_width,'spike_height',spike_height);
		
		if isempty(plot_trials) | size(plot_trials,1)<i
			ylimits=[min(SPIKES(i).trial) max(SPIKES(i).trial)];
		else
			ylimits=[plot_trials(i,1) plot_trials(i,2)];
		end

		if diff(ylimits)>0
			ylim(ylimits);
		else
			ylimits=ylim();
		end
		box off;
		set(gca,'ydir','rev','ytick',[],'yticklabel',ylimits-(ylimits(1)-1),'tickdir','out','ticklength',[0 0]);

		if i<length(SPIKES)
			set(gca,'xtick',[]);
		end
	end
else

		if isempty(fs)
			spike_fs=SPIKES.fs;
		else
			spike_fs=fs;
		end

		ax(2)=subplot(nplots,1,2:nplots);
		spikoclust_raster(SPIKES.times/spike_fs,SPIKES.trial,'spike_width',spike_width,...
            'spike_height',spike_height);

		if isempty(plot_trials)
			ylimits=[min(SPIKES.trial) max(SPIKES.trial)];
		else
			ylimits=[plot_trials(1) plot_trials(2)];
		end

		if diff(ylimits)>0
			ylim(ylimits);
		else
			ylimits=ylim();
		end
		box off;
		set(gca,'ydir','rev','ytick',[],'tickdir','out','ticklength',[0 0]);
		%xlabel('Time (s)');
		linkaxes(ax,'x');

end

ylimits=ylim();
ylim_range=range(ylimits);
ylimits(2)+ylim_range*.01;
linkaxes(ax,'x');
xlim([SPECT.t(1) SPECT.t(end)]);
