function stan_ephys_plot(STATS)
% ephys stats plots
% stability analysis--baseline data
% first take all of the data from control

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
cmap=parula(length(birds));

breaks=[];

% first break is 30-49

idx1=find(days_since>40&days_since<100);
idx2=find(days_since>100);

% translocate 49 to 35

offset1=14;
offset2=72;

break1=32;
break2=41;

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

	stan_plot_dot_error(days_since,plot_stats.(to_plot{i}).mu,plot_stats.(to_plot{i}).mu_ci,birdid);
	ylabel(y_labels{i});
	if i==length(to_plot)
		xlabel('Days');
	else
		set(gca,'XTick',[]);
	end

	ylimits=ylim();
	offset=eps;

	h=patch([ break1 break1+1 break1+1 break1 ],...
		[ ylimits(2)+offset ylimits(2)+offset ylimits(1)-offset ylimits(1)-offset ],1,'facecolor',[1 1 1],'edgecolor','k');
	h=patch([ break2 break2+1 break2+1 break2 ],...
		[ ylimits(2)+offset ylimits(2)+offset ylimits(1)-offset ylimits(1)-offset ],1,'facecolor',[1 1 1],'edgecolor','k');
	ylim([ylimits]);

	if i==length(to_plot)
		set(gca,'XTick',[0 5 10 15 20 25 30 33 38 42],'XTickLabel',[0:5:30 [33 38]+offset1 42+offset2]);
	end

end

linkaxes(ax,'x');
xlim([-1 47]);

set(fig,'position',[200 200 560 420],'paperpositionmode','auto');
markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),[ 'signal_timecourse' ],'eps,fig,png,pdf','renderer','painters');
