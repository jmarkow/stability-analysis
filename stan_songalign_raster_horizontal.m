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
y_label=1;

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
		case 'y_label'
			y_label=varargin{i+1};
	end
end

nrows=length(SPECT);
ax=[];

for i=1:length(SPECT)

	idx=((i-1)*2)+1;
	ax(idx)=subplot(nrows,2,idx);
	imagesc(SPECT(i).t,SPECT(i).f/1e3,SPECT(i).s);
	colormap(colors);
	axis xy;
	box off;
	set(gca,'TickDir','out','TickLength',[0 0]);
	ylim([0 9]);
	set(gca,'YTick',[0 9],'XTick',[]);
	title([name]);
	
	if y_label
		ylabel('kHz');
	else
		set(gca,'YTick',[]);
	end

	if isempty(fs)
		spike_fs=SPIKES(i).fs;
	else
		spike_fs=fs;
	end

	ax(idx+1)=subplot(nrows,2,idx+1);
	spikoclust_raster(SPIKES(i).times/spike_fs,SPIKES(i).trial,'spike_width',spike_width,...
		'spike_height',spike_height);	

	if isempty(plot_trials)
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

	set(gca,'ydir','rev','ytick',ylimits,'tickdir','out','ticklength',[0 0]);
	linkaxes(ax,'x');

	if ~y_label
		set(gca,'ytick',[]);
	end

	if i<length(SPECT)
		set(gca,'xtick',[]);
	else
		xlabel('Time (s)');
	end

end

linkaxes(ax,'x');
%xlim([SPECT.t(1) SPECT.t(end)]);

%pos=get(gcf,'position')
%set(gcf,'position',[ pos(1) pos(2) pos(3) pos(4) ]); % scale height
