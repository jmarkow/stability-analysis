function fig=stan_ephys_plot_correlation_comparison_smooth(BASELINE,NERVECUT)
% ephys stats plots
% stability analysis--baseline data
% first take all of the data from control

[options,dirs]=stan_preflight;

% collect spike and rms correlation

ci_inv=.01; % alpha level (/2 to get two-tailed)
ylimit_rounding=1e-1;
r_cutoff=.4; % below this value typically due to equipment failure
save_name='baseline_regression';
t_bins=[]; % bins for the baseline data to compute prctile, etc.
h_offset=.28;
v_offset=.095;
width=.35;
shade_color=[.7 .7 1];
save_name=['baseline_comparison_stats'];
% first column is days since 1, second is regression value, third is bird ID
% remove points where x=0 (by definition == 1, artifacts appear to have r<.4)

baseline_data.rms=[];
baseline_data.spikes=[];
baseline_data.birdid=[];
baseline_data.days_since=[];
baseline_data.idx=[];

for i=1:length(BASELINE.rms_corr)
	npoints=size(BASELINE.spikes_corr{i},2)-1;
	baseline_data.spikes=[baseline_data.spikes BASELINE.spikes_corr{i}(1,2:end)];
	baseline_data.rms=[baseline_data.rms BASELINE.rms_corr{i}(1,2:end)];
	baseline_data.days_since=[baseline_data.days_since BASELINE.days_since{i}(2:end)];
	baseline_data.birdid=[baseline_data.birdid ones(1,npoints)*i];
	baseline_data.idx=[baseline_data.idx 2:length(BASELINE.rms_corr{i}(1,:))];
end

baseline_data

artifacts=baseline_data.spikes<=r_cutoff;
artifacts2=isnan(baseline_data.spikes);

to_del=find(artifacts|artifacts2)

baseline_data.spikes(to_del)=[];
baseline_data.days_since(to_del)=[];
baseline_data.birdid(to_del)=[];
baseline_data.rms(to_del)=[];
baseline_data.idx(to_del)=[];

% bin data, create baseline region and plot nervecut data on top using dot and error plot

baseline_data
NERVECUT

%win_steps=[-inf 10:10:30 inf ];

% 16/5 for visualization, 25/0 for testing

win=16;
overlap=5;
win_step=win-overlap;

%win_steps(1)=-inf;
%win_steps(end)=inf;
%[~,bins]=histc(baseline_data.days_since,win_steps);
%[~,bins_nervecut]=histc(NERVECUT.days_since,win_steps);
condition={'spikes','rms'};

steps=0:win_step:25;

% collect p-value for each point

fig=figure();
ax=[];

% hypothesis testing

p.pboot.rms=nan(1,length(NERVECUT.spikes));
p.pboot.spikes=nan(1,length(NERVECUT.spikes));

for ii=1:2
	ax(ii)=subplot(2,1,ii);

	% form x and y data

	xdata=[];
	ydata=[];

	% sliding window?

	for i=1:length(steps)

		left_edge=steps(i)-win/2;
		left_edge(left_edge<0)=0;
		right_edge=steps(i)+win/2;
		right_edge(right_edge>max(steps))=max(steps);

		hits=find(baseline_data.days_since>=left_edge&baseline_data.days_since<right_edge);

		cur_bindata=baseline_data.(condition{ii})(hits);
		cur_bindata(isnan(cur_bindata))=[];

		cur_birdid=baseline_data.birdid(hits);
		cur_idx=baseline_data.idx(hits);

		cur_bootdata=[];

		for j=1:length(cur_birdid)
			cur_bootdata=[cur_bootdata;BASELINE.([condition{ii} '_corr_boot']){cur_birdid(j)}(:,cur_idx(j))];
		end

		% get points that map to this bin

		nervecut_idx=find(NERVECUT.days_since>=left_edge&NERVECUT.days_since<steps(i)+right_edge)

		for j=1:length(nervecut_idx)
			% ncomparisons (how do we want to account for this?)
			p.pboot.(condition{ii})(nervecut_idx(j))=1-mean(NERVECUT.(condition{ii})(nervecut_idx(j))<cur_bootdata);
		end

		% hypothesis testing

		bin_mu=median(cur_bindata);

		bin_conf(i,1)=prctile(cur_bootdata,.5); % bootstrap minimum?
		bin_conf(i,2)=prctile(cur_bootdata,100-.5);
		bin_conf(i,3)=prctile(cur_bootdata,50);

		p.bin_conf.(condition{ii})=bin_conf;
		p.bin_conf.idx(i,:)=[steps(i)-win/2 steps(i)+win/2];

		% how many bootstraps are we < | >?

		%bin_conf(i)=median(cur_bindata)-2*iqr(cur_bindata);

		%bin_conf(i)=min(cur_bootdata);

		xdata=[xdata steps(i)];
		ydata=[ydata bin_conf(i)];

	end

	% plot everything

	xdata(xdata==0)=-1;
	%area(xdata,ydata(2,:),0,'facecolor',shade_color,'edgecolor','none');
	markolab_shadeplot(xdata,bin_conf(:,[1 2])',shade_color,'k');
	hold on;
	plot(xdata,bin_conf(:,3),'k--');
	%plot(xdata,ydata_line,'k-','linewidth',1)
	%plot(xdata_boundary,ydata_boundary,'k-','linewidth',1);

	stan_plot_dot_error(NERVECUT.days_since,NERVECUT.(condition{ii}),NERVECUT.([condition{ii} '_ci']),NERVECUT.birdid,...
		'markersize',10);

	ylim([0 1]);
	set(gca,'YTick',[0 1],'layer','top','FontSize',7);
	yh=ylabel([ condition{ii} ]);

	if ii==2
		xlabel('Days');
	else
		set(gca,'XTick',[]);
	end

	set(yh,'position',get(yh,'position')+[2 0 0])

	%pos=get(ax(ii),'position');
	%asp_ratio=pos(3)/pos(4);
	%new_width=width/asp_ratio;

	%new_axis(ii)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);
	%%markolab_shadeplot(xdata,ydata,shade_color,'none');
	%area(xdata,ydata(2,:),0,'facecolor',shade_color,'edgecolor','none');
	%hold on;
	%plot(xdata,ydata_line,'k-','linewidth',1)
	%plot(xdata_boundary,ydata_boundary,'k-','linewidth',1);
	%stan_plot_dot_error(NERVECUT.days_since,NERVECUT.(condition{ii}),NERVECUT.([condition{ii} '_ci']),NERVECUT.birdid,...
	%	'markersize',10);
	%set(gca,'layer','top','ytick',[])
	%ylim([0 1]);
	%xlim([-1 40]);

	% now make inset


end


win=25;
overlap=0;
win_step=win-overlap;

condition={'spikes','rms'};

steps=0:win_step:25;

for ii=1:2
	ax(ii)=subplot(2,1,ii);


	% sliding window?

	for i=1:length(steps)

		left_edge=steps(i)-win/2;
		left_edge(left_edge<0)=0
		right_edge=steps(i)+win/2;
		right_edge(right_edge>max(steps))=max(steps)

		left_edge
		right_edge

		hits=find(baseline_data.days_since>=left_edge&baseline_data.days_since<right_edge);

		cur_bindata=baseline_data.(condition{ii})(hits);
		cur_bindata(isnan(cur_bindata))=[];

		cur_birdid=baseline_data.birdid(hits);
		cur_idx=baseline_data.idx(hits);

		cur_bootdata=[];

		for j=1:length(cur_birdid)
			cur_bootdata=[cur_bootdata;BASELINE.([condition{ii} '_corr_boot']){cur_birdid(j)}(:,cur_idx(j))];
		end

		% get points that map to this bin

		nervecut_idx=find(NERVECUT.days_since>=left_edge&NERVECUT.days_since<steps(i)+right_edge)

		for j=1:length(nervecut_idx)
			% ncomparisons (how do we want to account for this?)
			p.pboot.(condition{ii})(nervecut_idx(j))=1-mean(NERVECUT.(condition{ii})(nervecut_idx(j))<cur_bootdata);

		end


	end

end


% get all the p-values

fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,'fig3_drift.txt'),'w+');
pvals=sort(markolab_bonf_holm([p.pboot.spikes(NERVECUT.days_since<25) p.pboot.rms(NERVECUT.days_since<25)],.5));
fprintf(fid,'p=');
for i=1:length(pvals)-1
    fprintf(fid,'%g,',pvals(i));
end

fprintf(fid,'%g',pvals(end));
fclose(fid);


set(fig,'units','centimeters','position',[3 3 8.4 12],'paperpositionmode','auto');

linkaxes(ax,'xy');
set(ax(1),'xlim',[-1 25]);
ylim([0 1.01]);
%linkaxes(new_axis,'xy');
%set(new_axis(1),'xlim',[-2 40]);

% do we want to highlight p<.01 with Bonferonni?

%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'baseline_nervecut_comparison','eps,fig,png,pdf','renderer','painters');
save(fullfile(dirs.agg_dir,dirs.datastore_dir,[ save_name '.mat']),'baseline_data','p');
