function stan_plot_lfp_angdiff(LFP_DATA)
%
%
%
%

% mle estimate of log-normal mode (peak of PDF seems most intuitive)

[options,dirs]=stan_preflight;

logmode=@(x) exp(mean(log(x))-std(log(x))^2);
logmedian=@(x) exp(mean(log(x)));
filewrite=true;
time_cuts=[0 30;...
	31 511];...

ntimebins=size(time_cuts,1);

bins=[0:pi/20:pi];

colors=[0 0 1;...
	1 0 0;...
	0 1 1];

ax=[];
labels={};

if filewrite
    fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,'figs1_lfp_drift.txt'),'w+');
    
end


for i=1:ntimebins

	dat=LFP_DATA.ang_diff(LFP_DATA.days_since>time_cuts(i,1)&LFP_DATA.days_since<=time_cuts(i,2));
	%est=ksdensity(dat,bins,'support',[0 pi]);

	est=histc(dat,bins);
	est=est./sum(est);
	mu=logmode(dat)
	ci=bootci(1e3,{logmode,dat},'type','cper','alpha',.01)
	
	%ax(i)=area(bins,est,'facecolor',colors(i,:),'edgecolor','k');
	ax(i)=markolab_stairplot(est,bins,'facecolor',colors(i,:),'edgecolor','k','method','a');
	hold on;
	
	% find the closets estimate to mu and plot from 0 to tip of distribution
	
	[~,mindist]=min(abs(mu-bins));

	% use this y point 
	
	plot([mu mu],[0 est(mindist-1)],'k-')
	%plot(repmat(ci(:)',[2 1]),repmat([0;est(mindist)],[1 2]),'k--')

	labels{end+1}=sprintf('%i-%i',time_cuts(i,1),time_cuts(i,2));

    if filewrite
        fprintf(fid,'modal change %g ci (%g-%g) (%i-%i days)\n',mu,ci(1),ci(2),time_cuts(i,1),time_cuts(i,2));
    end
    
end

if filewrite
    fprintf(fid,'n(channels) %i n(birds) %i',length(unique(LFP_DATA.channel_id)),length(unique(LFP_DATA.bird_id)));
    fclose(fid);
end
L=legend(ax,labels);
legend boxoff;
set(L,'FontSize',6);
set(gca,'FontSize',7,'FontName','Helvetica','TickLength',[0 0]);



%findall(L,'type','area')

