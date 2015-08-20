function fig=stan_ephys_plot(STATS)
% ephys stats plots
% stability analysis--baseline data
% first take all of the data from control

fontsize=7;

offset1=14;
offset2=72;

break1=32;
break2=41;

markersize=10;

[options,dirs]=stan_preflight;

% collapse data into plot-able vectors

% first spike rate across time

mufun=@(x) median(x);
mufun_ci=@(x) iqr(x);

% first dot/errorbar plot, then regress (?)

% first spikerate

to_plot={'rms_mu','spikes_mu','threshold'};
y_labels={'\langleRMS\rangle (\muV)','\langleFR\rangle (Hz)','Thresh (\muV)'};

spikes_mu=[];
spikes_mu_ci=[];
days_since=[];
birdid=[];

%for i=1:length(STATS.spikes.mu)
%	days_since=[days_since STATS.days_since{i}];

	%birdid=[birdid ones(size(mu))*i];
%end

for i=1:length(to_plot)

	plot_stats.(to_plot{i}).mu=[];
	plot_stats.(to_plot{i}).mu_ci=[];
	
	for j=1:length(STATS.(to_plot{i}))
		mu=cellfun(mufun,STATS.(to_plot{i}){j});
		sig=cellfun(mufun_ci,STATS.(to_plot{i}){j});
		plot_stats.(to_plot{i}).mu=[plot_stats.(to_plot{i}).mu mu];
		plot_stats.(to_plot{i}).mu_ci=[plot_stats.(to_plot{i}).mu_ci [mu+sig;mu-sig]];
	end
end

sz=cellfun(@length,STATS.rms_mu);
birdid=[];
days_since=[];

for i=1:length(sz)
	days_since=[days_since STATS.days_since{i}];
	birdid=[birdid ones(1,sz(i))*i];
end

% collect per trial variables


to_del=find(plot_stats.rms_mu.mu<10)

birdid(to_del)=[];
days_since(to_del)=[];

for i=1:length(to_plot)
	plot_stats.(to_plot{i}).mu(to_del)=[];
	plot_stats.(to_plot{i}).mu_ci(:,to_del)=[];
end


% break up x-axis 

birds=unique(birdid);

tmp=birdid;

for i=1:length(birds)
	tmp(birdid==birds(i))=i;
end

birdid=tmp;
birds=unique(tmp);
cmap=paruly(length(birds));

breaks=[];

% first break is 30-49

idx1=find(days_since>40&days_since<100);
idx2=find(days_since>100);

% translocate 49 to 35



days_since(idx1)=days_since(idx1)-offset1;
days_since(idx2)=days_since(idx2)-offset2;

% make a patch to mark the break

% supplementary stats (mean rms, sigma rms)

fig=figure();
ax=[];
for i=1:length(to_plot)

	i
	ax(i)=subplot(length(to_plot),1,i);

	% any outliers?	

	for j=1:length(birds)
		plot(days_since(birdid==birds(j)),plot_stats.(to_plot{i}).mu(birdid==birds(j)),'-','color',cmap(j,:));
		hold on;
	end

	stan_plot_dot_error(days_since,plot_stats.(to_plot{i}).mu,plot_stats.(to_plot{i}).mu_ci,birdid,'markersize',markersize);
	yh=ylabel(y_labels{i});
	%set(yh,'position',get(yh,'position')+[ 0 0]);
	if i==length(to_plot)
		xlabel('Days');
	else
		set(gca,'XTick',[]);
	end

	ylimits=ylim();
	offset=1;

	% boxes

	%h=patch([ break1 break1+1 break1+1 break1 ],...
	%	[ ylimits(2)+offset ylimits(2)+offset ylimits(1)-offset ylimits(1)-offset ],1,'facecolor',[1 1 1],'edgecolor','none');
	%h2=patch([ break2 break2+1 break2+1 break2 ],...
	%	[ ylimits(2)+offset ylimits(2)+offset ylimits(1)-offset ylimits(1)-offset ],1,'facecolor',[1 1 1],'edgecolor','none');
	%set(h2,'clipping','off');
	%set(h,'clipping','off');
	
	
	%plot([break1 break1],ylimits,'k-');
	%plot([break1+1 break1+1],ylimits,'k-');

	%plot([break2 break2],ylimits,'k-');
	%plot([break2+1 break2+1],ylimits,'k-');
	
	ylim([ylimits]);

	% plot lines on top of white patches to indicate breaks

	% wavy lines

	splitwidth=1;
	theta = linspace(0,2*pi,100);
	amp = splitwidth/2 * 0.9;
    x1= amp * sin(theta) + break1;
	x2= amp * sin(theta) + break1+1;
	ypoints=linspace(ylimits(1)-offset,ylimits(2)+offset,100);

	h=patch([x1 fliplr(x2)],[ypoints fliplr(ypoints)],1,'facecolor',[1 1 1],'edgecolor','none');

	plot(x1,ypoints,'k-');
	plot(x2,ypoints,'k-');

	x1= amp * sin(theta) + break2;
	x2= amp * sin(theta) + break2+1;

	h2=patch([x1 fliplr(x2)],[ypoints fliplr(ypoints)],1,'facecolor',[1 1 1],'edgecolor','none');

	set(h,'clipping','off');
	set(h2,'clipping','off');
	set(gca,'clipping','off');
	plot(x1,ypoints,'k-');
	plot(x2,ypoints,'k-');


	if i==length(to_plot)
		set(gca,'XTick',[0 5 10 15 20 25 33 42],'XTickLabel',[0:5:30 33+offset1 42+offset2]);
	end

	set(gca,'fontsize',fontsize,'ytick',[ylimits(1) range(ylimits)/2+ylimits(1) ylimits(2)],'layer','bottom');
end

linkaxes(ax,'x');
xlim([-1 47]);


%set(fig,'units','centimeters','position',[3 3 8 12],'paperpositionmode','auto');
%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),[ 'signal_timecourse' ],'eps,fig,png,pdf','renderer','painters');
