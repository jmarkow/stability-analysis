function stan_fig2g(LFP)
% Generates Fig 2g, example song-averaged LFPs over long time-scales
%
% note: sonograms are computed using the zftftb library found at https://github.com/jmarkow/zftftb

scaling_fun=@(x) (x/1.5)*5;

% which options specify control raster

for i=1:length(LFP)

	% grab the template

	[b,a]=ellip(4,.2,40,[25 35]/(LFP(i).lfp_fs/2),'bandpass');

	filt_lfp=cellfun(@(x) filtfilt(b,a,x),LFP(i).lfp_data,'uniformoutput',0);
	mu_lfp=cellfun(@(x) mean(zscore(x),2),filt_lfp,'uniformoutput',0);
	mu_lfp=zscore(cat(2,mu_lfp{:}));

	spacing=6;
	spacing=repmat(spacing,[1 length(LFP(i).days)-1]);
	spacing=[0 cumsum(spacing)];

	[nsamples,ntrials]=size(mu_lfp);

	spacing=repmat(spacing,[nsamples 1]);
	mu_t=[1:nsamples]/LFP(i).lfp_fs;

	fig.(LFP(i).bird_id)=figure();

	[s,f,t]=zftftb_pretty_sonogram(LFP(i).audio_data,LFP(i).audio_fs,...
		'filtering',300,'len',70,'overlap',69.5,'clipping',[-3 2],'zeropad',0);
	ax(1)=subplot(3,1,1);
	imagesc(t+LFP(i).padding(1),f/1e3,s);
	axis xy;
	ylim([0 8]);
	set(gca,'XTick',[],'YTick',[]);

	ax(2)=subplot(3,1,2:3);
	plot(mu_t,mu_lfp+spacing,'b-');

	linkaxes(ax,'x');

	if isfield(LFP(i),'xlim')
		xlim([LFP(i).xlim]);
	end

end

% scale in a manner that keeps time consistent

ax=get(fig.lhp33,'currentAxes');
xrange1=range(get(ax,'xlim'));

ax=get(fig.rm7,'currentAxes');
xrange2=range(get(ax,'xlim'));

% rescale axes (can't make figure small enough here)

set(fig.lhp33,'units','centimeters','position',[3 3 6 6.5],'paperpositionmode','auto');
set(fig.rm7,'units','centimeters','position',[3 3 6 6.5],'paperpositionmode','auto');

ax=findall(fig.lhp33,'type','axes');

for i=1:length(ax)

	set(ax(i),'units','centimeters');
	pos=get(ax(i),'position');

	% width change

	new_width=scaling_fun(xrange1);
	width_change=new_width-pos(3);
	set(ax(i),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);
	
end

set(ax(1),'ylim',[-5 23],'xtick',[],'ytick',[]);
xlimits=get(ax(1),'xlim');
h=line([xlimits(1) xlimits(1)+.2],[-7 -7],'parent',ax(1));
h2=line([xlimits(1)-.1*xrange1 xlimits(1)-.1*xrange1],[-5 0],'parent',ax(1));
set(h,'clipping','off');
set(h2,'clipping','off');

ax=findall(fig.rm7,'type','axes');

for i=1:length(ax)

	set(ax(i),'units','centimeters');
	pos=get(ax(i),'position');

	% width change

	new_width=scaling_fun(xrange2);
	width_change=new_width-pos(3);

	set(ax(i),'position',[pos(1)-width_change/2 pos(2) new_width pos(4)]);

end

set(ax(1),'ylim',[-5 23],'xtick',[],'ytick',[]);
xlimits=get(ax(1),'xlim');
h=line([xlimits(1) xlimits(1)+.2],[-7 -7],'parent',ax(1));
h2=line([xlimits(1)-.1*xrange2 xlimits(1)-.1*xrange2],[-5 0],'parent',ax(1));
set(h,'clipping','off');
set(h2,'clipping','off');
