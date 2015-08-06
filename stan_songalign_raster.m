function stan_songalign_raster(SPECT,SPIKES,varargin)
%
%
%
%
%


colors='jet';
nparams=length(varargin);
name='';

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'colors'
			colors=varargin{i+1};
		case 'name'
			name=varargin{i+1};
	end
end

ax(1)=subplot(4,1,1);
imagesc(SPECT.t,SPECT.f,SPECT.s);
colormap(colors);
axis xy;
box off;
set(gca,'TickDir','out','TickLength',[0 0]);
ylim([0 9]);
set(gca,'YTick',[0 9],'XTick',[]);
title([name]);
ylabel('Fs (kHz)');

ax(2)=subplot(4,1,2:4);
spikoclust_raster(SPIKES.times/SPIKES.fs,SPIKES.trial);
ylimits=[min(SPIKES.trial) max(SPIKES.trial)];

if diff(ylimits)>0
	ylim(ylimits);
else
	ylimits=ylim();
end
box off;
set(gca,'ydir','rev','ytick',ylimits,'tickdir','out','ticklength',[0 0]);
xlabel('Time (s)');
linkaxes(ax,'x');
xlim([SPECT.t(1) SPECT.t(end)]);

%pos=get(gcf,'position')
%set(gcf,'position',[ pos(1) pos(2) pos(3) pos(4) ]); % scale height
