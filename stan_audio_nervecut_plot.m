function [to_plot,p]=stan_nervecut_audio_plot(PLOT_FEATURES)
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
	p.allbird.mu.(plot_features{i})=signrank(cat(2,PLOT_FEATURES.mu(:).(plot_features{i})));
	p.allbird.var.(plot_features{i})=signrank(cat(2,PLOT_FEATURES.var(:).(plot_features{i})));
end


% get mu/variance and plot appropriately

counter=1;
for i=1:nfeatures
	to_plot{counter}=cat(2,PLOT_FEATURES.mu(:).(plot_features{i}));
	counter=counter+1;
	to_plot{counter}=cat(2,PLOT_FEATURES.var(:).(plot_features{i}));
	counter=counter+1;
end

fig=figure();
repmat(1:2,[1 nfeatures])
to_plot
[h.box,h.med,h.whisk]=markolab_boxplot(to_plot,[],'feature_idx',repmat(1:2,[1 nfeatures]),'bunching',2,'bunching_offset',3);
xlimits=xlim();
ylimits=ylim();
plot([xlimits(1) xlimits(2)],[0 0],'k--');
pos=get(gca,'XTick');
new_pos=[];
for i=1:2:nfeatures*2
	new_pos(end+1)=mean(pos(i):pos(i+1))
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
		newxpos=xdata(1)-width*2;
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
		newxpos=xdata(1)-width*2;
		set(h_txt,'position',[newxpos pos(2:3)]);
	end
	
	counter=counter+1;
end


set(gca,'XTick',new_pos,'XTickLabel',plot_labels,'TickDir','out','TickLength',[0 0]);
ylabel('Post-cut change in feature (Z)');
set(fig,'position',[200 200 180 180],'paperpositionmode','auto');
markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),save_name,'eps,png,fig,pdf','renderer','painters');
save(fullfile(dirs.agg_dir,dirs.fig_dir,[save_name '.mat']),'p');
