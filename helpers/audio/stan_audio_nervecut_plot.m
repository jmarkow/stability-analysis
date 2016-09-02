function [fig,p]=stan_nervecut_audio_plot(PLOT_FEATURES)
% takes collected stats, plots and performs hypothesis tests
% 
%
%
save_name='nervecut_acoustic_features';

% get options

[options,dirs]=stan_preflight;

plot_features={'AM','FM','entropy','pitch_goodness'};
plot_labels={'AM','FM','Ent.','PG'};
alpha=[.001 .01 .05];
alpha_sign={ '***','**','*' };
feature_names=fieldnames(PLOT_FEATURES.mu);
nfeatures=length(plot_features);

%%%% hypothesis tests
%

nbirds=length(PLOT_FEATURES.mu)

for i=1:nbirds
	for j=1:nfeatures
		tmp=signrank(PLOT_FEATURES.mu(i).(plot_features{j}));
		p.perbird.mu(i).(plot_features{j})=tmp;
		tmp=signrank(PLOT_FEATURES.var(i).(plot_features{j}));
		p.perbird.var(i).(plot_features{j})=tmp;
	end
end

for i=1:nfeatures
	p.allbird.mu.(plot_features{i})=signrank(abs(cat(2,PLOT_FEATURES.mu(:).(plot_features{i}))));
	p.allbird.var.(plot_features{i})=signrank(abs(cat(2,PLOT_FEATURES.var(:).(plot_features{i}))));
end


% get mu/variance and plot appropriately

to_plot_idx=cell(1,nfeatures*2);
to_plot=cell(1,nfeatures*2);
counter=1;
nbirds=length(PLOT_FEATURES.mu);

for i=1:nfeatures
	for j=1:nbirds
		to_plot{counter}=[to_plot{counter} PLOT_FEATURES.mu(j).(plot_features{i})];
		to_plot_idx{counter}=[to_plot_idx{counter} ones(size(PLOT_FEATURES.mu(j).(plot_features{i})))*j];
	end
	counter=counter+1;
	
	for j=1:nbirds	
		to_plot{counter}=[to_plot{counter} PLOT_FEATURES.var(j).(plot_features{i})];
		to_plot_idx{counter}=[to_plot_idx{counter} ones(size(PLOT_FEATURES.var(j).(plot_features{i})))*j];
	end
	counter=counter+1;
end

% sort to make useful for beeswarm plot

swarm_plot_mu=cell(1,nfeatures*nbirds);
counter=1;
for i=1:nfeatures
	for j=1:nbirds
		swarm_plot_mu{counter}=PLOT_FEATURES.mu(j).(plot_features{i});
		counter=counter+1;
	end
end

swarm_plot_var=cell(1,nfeatures*nbirds);
counter=1;
for i=1:nfeatures
	for j=1:nbirds
		swarm_plot_var{counter}=PLOT_FEATURES.var(j).(plot_features{i});
		counter=counter+1;
	end
end

swarm_colors=repmat(parula(nbirds),[nfeatures 1])
offset=2;
pos=ones(1,length(swarm_plot_mu));
pos(1:nbirds:end)=pos(1:nbirds:end)+offset;
pos=cumsum(pos)

ax=[];

fig.swarm=figure();
ax(1)=subplot(2,1,1);
h=plotSpread(swarm_plot_mu,'xValues',pos,'binWidth',.6,'distributionColors',swarm_colors);

% xticks are grouped according to pos

xticks=[];
for i=1:nbirds:length(swarm_plot_mu)
	xticks=[xticks mean(pos(i):pos(i+nbirds-1))];
end

set(h{1},'markersize',2.5);
set(gca,'xtick',[],'TickLength',[0 0],'TickDir','out','FontSize',7);
xlimits=xlim();
hold on;
h=plot([xlimits],[0 0],'k-');
uistack(h,'bottom');
ylim([-10 10]);
ylimits=ylim();
set(gca,'YTick',[ylimits(1) 0 ylimits(2)]);
ylabel('\langleFeature\rangle (Z)')

ax(2)=subplot(2,1,2);
h=plotSpread(swarm_plot_var,'xValues',pos,'binWidth',.6,'distributionColors',swarm_colors);
set(h{1},'markersize',2.5);
set(gca,'xtick',xticks,'xticklabel',plot_labels,'TickLength',[0 0],'TickDir','out','FontSize',7);
hold on;
h=plot([xlimits],[0 0],'k-');
uistack(h,'bottom');
ylabel('\sigma(Feature) (Z)');
ylim([-13 8]);
ylimits=ylim();
set(gca,'YTick',[ylimits(1) 0 ylimits(2)]);
linkaxes(ax,'x');
% rectify

for i=1:length(to_plot)
	to_plot{i}=abs(to_plot{i});
end

fig.box_plot=figure();
repmat(1:2,[1 nfeatures])
[h.box,h.med,h.whisk]=markolab_boxplot(to_plot,[],'feature_idx',repmat(1:2,[1 nfeatures]),'bunching',2,'bunching_offset',3);
xlimits=xlim();
ylimits=ylim();
%plot([xlimits(1) xlimits(2)],[0 0],'k--');
pos=get(gca,'XTick');
new_pos=[];

for i=length(h.box)-1:1
    x1=get(h.box(i),'xdata');
    x2=get(h.box(i+1),'xdata');
    new_pos(end+1)=mean([x1([1 3]) x2([1 3])])
end

for i=1:2:nfeatures*2
	%new_pos(end+1)=mean(pos(i):pos(i+1))
end

% plot mu sig tests

counter=1;

for i=length(h.whisk):-2:(length(h.whisk)/2+1)
	ydata=get(h.whisk(i),'ydata');
	xdata=get(h.whisk(i),'xdata');
	yext=double(ydata(2)+(.01*range(ylimits)));

	p_hit=(p.allbird.mu.(plot_features{counter})*nfeatures)<alpha
	idx=min(find(p_hit));

	if ~isempty(idx) & idx>0
		h_txt=text(xdata(1),yext,alpha_sign{idx});
		pos=get(h_txt,'position');
		ext=get(h_txt,'extent');
		width=ext(3);
		% center over line position	
		newxpos=xdata(1)-1.5*width;
		set(h_txt,'position',[newxpos pos(2:3)]);
	end
	
	counter=counter+1;
end

counter=1;

for i=length(h.whisk)-1:-2:(length(h.whisk)/2+1)
	ydata=get(h.whisk(i),'ydata');
	xdata=get(h.whisk(i),'xdata');
	yext=double(ydata(2)+(.01*range(ylimits)));

	p_hit=(p.allbird.var.(plot_features{counter})*nfeatures)<alpha
	idx=min(find(p_hit));

	if ~isempty(idx) & idx>0
		h_txt=text(xdata(1),yext,alpha_sign{idx});
		pos=get(h_txt,'position');
		ext=get(h_txt,'extent');
		width=ext(3);
		% center over line position	
		newxpos=xdata(1)-1.5*width;
		set(h_txt,'position',[newxpos pos(2:3)]);
	end
	
	counter=counter+1;
end


ylim([-.1 9]);
set(gca,'XTick',new_pos,'XTickLabel',plot_labels,'TickDir','out','TickLength',[0 0]);
ylabel('Abs. change (Z)');
%set(fig,'position',[200 200 180 180],'paperpositionmode','auto');
%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),save_name,'eps,png,fig,pdf','renderer','painters');
%save(fullfile(dirs.agg_dir,dirs.fig_dir,[save_name '.mat']),'p');
%


