function stan_songalign_raster(SPECT,RMS_DATA,RMS_FS,varargin)
%
%
%
%
%


spect_colors='jet';
rms_colors='parula';
nparams=length(varargin);
name='';
clipping=[2.5 97.5];
clim_order=[];
cbar_dist=.006;
cbar_width = .025; 
cbar_label='';
offset=0;

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'spect_colors'
			spect_colors=varargin{i+1};	
		case 'rms_colors'
			rms_colors=varargin{i+1};
		case 'name'
			name=varargin{i+1};
		case 'spect_colors'
			spect_colors=varargin{i+1};
		case 'clipping'
			clipping=varargin{i+1};
		case 'clim_order'
			clim_order=varargin{i+1};
		case 'cbar_width'
			cbar_width=varargin{i+1};
		case 'cbar_label'
			cbar_label=varargin{i+1};
		case 'offset'
			offset=varargin{i+1};
	end
end


ax(1)=subplot(4,1,1);
imagesc(SPECT.t,SPECT.f,SPECT.s);
colormap(spect_colors);
axis xy;
box off;
set(gca,'TickDir','out','TickLength',[0 0]);
ylim([0 9]);
set(gca,'YTick',[0 9],'XTick',[]);
title([name]);
ylabel('Fs (kHz)');
freezeColors();

[nsamples,ntrials]=size(RMS_DATA);
t_vec=[1:nsamples]/RMS_FS;


ax(2)=subplot(4,1,2:4);
imagesc(t_vec,[],RMS_DATA');
colormap(rms_colors);
ylimits=[1 ntrials];
ylim(ylimits);
box off;
set(gca,'ydir','rev','ytick',ylimits,'tickdir','out','ticklength',[0 0]);
xlabel('Time (s)');

clips=prctile(RMS_DATA(:),clipping);
caxis([clips(1) clips(2)]);

linkaxes(ax,'x');
xlim([SPECT.t(1) SPECT.t(end)]);

if offset>0
	xlim([offset SPECT.t(end)]);
end

pos=get(ax(2),'position');
clims=caxis();

if isempty(clim_order)
	clims=[0 0];
else
	clim_order=1/clim_order;
	clims(1)=floor(clims(1)*clim_order)/clim_order;
	clims(2)=ceil(clims(2)*clim_order)/clim_order;
	caxis(clims);
end

if diff(clims)<=0
	clims=caxis();
	order=floor(log(abs(clims))./log(10));
	order(order==-inf)=[];
	[~,loc]=max(abs(order));
	fprintf('Detected order %g\n',order(loc));
	clim_order=1/(10^order(loc));
	clims=round(clims*clim_order)/clim_order;
end



for i=1:length(ax)
	pos=get(ax(i),'position');
	set(ax(i),'position',[pos(1) pos(2) pos(3)-cbar_width-cbar_dist pos(4)]);
end

for i=1:length(clims)
	xticklabels{i}=sprintf('%g',clims(i));
end

newpos=get(ax(2),'position');

hc = colorbar('location','eastoutside','position',...
	[newpos(1)+newpos(3)+cbar_width+cbar_dist newpos(2) cbar_width newpos(4)],'fontsize',10,'xtick',clims,...
	'xticklabel',xticklabels);

cpos=get(hc,'position');
clab=get(hc,'ylabel');
set(clab,'string',cbar_label);
tpos=get(clab,'position');

set(clab,'units','normalized','position',[1.5 .5 0]);
