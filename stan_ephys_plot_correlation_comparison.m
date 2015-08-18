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
		cur_bin(cur_bin==-inf)=0;
		cur_bin(cur_bin==inf)=150;

		markolab_shadeplot([cur_bin],[1 1;bin_conf(i) bin_conf(i)],[.95 .95 .95],'none');
		hold on;
		plot([cur_bin],[bin_conf(i) bin_conf(i)],'k-','linewidth',1.5)

	end

	for i=1:length(win_steps)-2
		plot([win_steps(i+1) win_steps(i+1)],[bin_conf(i) bin_conf(i+1)],'k-','linewidth',1.5)
	end

	stan_plot_dot_error(NERVECUT.days_since,NERVECUT.(condition{ii}),NERVECUT.([condition{ii} '_ci']),NERVECUT.birdid);
	ylim([0 1]);
	ylabel([ condition{ii} ' correlation (R) ']);
	
	if ii==2
		xlabel('Days');
	else
		set(gca,'XTick',[]);
	end
	set(gca,'YTick',[0:.2:1]);
end

linkaxes(ax,'x');
xlim([0 150]);
set(fig,'position',[200 200 150 300],'paperpositionmode','auto');

% do we want to highlight p<.01 with Bonferonni?

%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'baseline_nervecut_comparison','eps,fig,png,pdf','renderer','painters');
%save(fullfile(dirs.agg_dir,dirs.fig_dir,[ save_name '.mat']),'baseline_data');

