function [baseline_data,p]=stan_ephys_plot(BASELINE,NERVECUT)
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
shade_color=[.85 .65 .65];
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

win_steps=[-inf 20 40 60 inf ];
%win_steps(1)=-inf;
%win_steps(end)=inf;
[~,bins]=histc(baseline_data.days_since,win_steps);
[~,bins_nervecut]=histc(NERVECUT.days_since,win_steps);
condition={'rms','spikes'};

% collect p-value for each point

fig=figure();
ax=[];

% hypothesis testing

p.binned.rms=nan(1,length(NERVECUT.spikes));
p.binned.spikes=nan(1,length(NERVECUT.spikes));
p.all.rms=nan(size(p.binned.rms));
p.all.spikes=nan(size(p.binned.rms));
p.days_since.rms=nan(size(p.binned.rms));
p.days_since.spikes=nan(size(p.binned.rms));

for ii=1:2
	ax(ii)=subplot(2,1,ii);

	for i=1:length(NERVECUT.(condition{ii}))
		p.all.(condition{ii})(i)=ranksum(baseline_data.(condition{ii}),NERVECUT.(condition{ii})(i),'tail','right');
		
		% get all baseline points with days since <= to current point

		cur_days_since=NERVECUT.days_since(i);
		idx=find(baseline_data.days_since<=cur_days_since)

		if isempty(idx)
			continue;
		end

		p.days_since.(condition{ii})(i)=ranksum(baseline_data.(condition{ii})(idx),NERVECUT.(condition{ii})(i),'tail','right');
	end

	% form x and y data
	
	xdata=[];
	ydata=[];
	ydata_line=[];
	xdata_boundary=[];
	ydata_boundary=[];

	for i=1:length(win_steps)-1

		cur_bin=[win_steps(i) win_steps(i+1)];
		cur_bindata=baseline_data.(condition{ii})(bins==i);
		cur_bindata(isnan(cur_bindata))=[];

		find(bins==i) % get bin_idx

		cur_birdid=baseline_data.birdid(bins==i);
		cur_idx=baseline_data.idx(bins==i);

		cur_bootdata=[];

		for j=1:length(cur_birdid)
			cur_bootdata=[cur_bootdata;BASELINE.([condition{ii} '_corr_boot']){cur_birdid(j)}(:,cur_idx(j))];
		end

		% get points that map to this bin
		
		nervecut_idx=find(bins_nervecut==i)

		for j=1:length(nervecut_idx)

			% ncomparisons (how do we want to account for this?)	
			p.binned.(condition{ii})(nervecut_idx(j))=ranksum(cur_bindata,NERVECUT.(condition{ii})(nervecut_idx(j)),'tail','right');
		
		end

		% hypothesis testing

		bin_mu=median(cur_bindata);
		%ci=bootci(1e3,{@median,cur_bindata});
		bin_conf(i)=prctile(cur_bootdata,1);
		%bin_conf(i)=median(cur_bindata)-2*iqr(cur_bindata);
		cur_bin(cur_bin==-inf)=-5;
		cur_bin(cur_bin==inf)=150;

		xdata=[xdata cur_bin];
		ydata=[ydata [1 1;bin_conf(i) bin_conf(i)]];
		ydata_line=[ydata_line bin_conf(i) bin_conf(i)];

	end

	for i=1:length(win_steps)-2
		xdata_boundary=[xdata_boundary win_steps(i+1) win_steps(i+1)];
		ydata_boundary=[ydata_boundary bin_conf(i) bin_conf(i+1)];
	end

	% plot everything

	area(xdata,ydata(2,:),0,'facecolor',shade_color,'edgecolor','none');
	%markolab_shadeplot(xdata,ydata,shade_color,'none');
	hold on;
	plot(xdata,ydata_line,'k-','linewidth',1)
	plot(xdata_boundary,ydata_boundary,'k-','linewidth',1);

	stan_plot_dot_error(NERVECUT.days_since,NERVECUT.(condition{ii}),NERVECUT.([condition{ii} '_ci']),NERVECUT.birdid,...
		'markersize',10);
	ylim([0 1]);
	ylabel([ condition{ii} ' (R)']);

	if ii==2
		xlabel('Days');
	else
		set(gca,'XTick',[]);
	end

	set(gca,'YTick',[0:.2:1],'layer','top');

	pos=get(ax(ii),'position');
	asp_ratio=pos(3)/pos(4);
	new_width=width/asp_ratio;

	new_axis(ii)=axes('position',[ pos(1)+pos(3)-h_offset pos(2)+pos(4)-v_offset width new_width ]);
	%markolab_shadeplot(xdata,ydata,shade_color,'none');
	area(xdata,ydata(2,:),0,'facecolor',shade_color,'edgecolor','none');
	hold on;
	plot(xdata,ydata_line,'k-','linewidth',1)
	plot(xdata_boundary,ydata_boundary,'k-','linewidth',1);
	stan_plot_dot_error(NERVECUT.days_since,NERVECUT.(condition{ii}),NERVECUT.([condition{ii} '_ci']),NERVECUT.birdid,...
		'markersize',10);
	set(gca,'layer','top','ytick',[])
	ylim([0 1]);
	xlim([-1 40]);

	% now make inset


end

set(fig,'units','centimeters','position',[3 3 8.4 12],'paperpositionmode','auto');

linkaxes(ax,'xy');
set(ax(1),'xlim',[-2 150]);

linkaxes(new_axis,'xy');
set(new_axis(1),'xlim',[-2 40]);

% do we want to highlight p<.01 with Bonferonni?

%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'baseline_nervecut_comparison','eps,fig,png,pdf','renderer','painters');
save(fullfile(dirs.agg_dir,dirs.fig_dir,[ save_name '.mat']),'baseline_data','p');

