function baseline_data=stan_ephys_plot(BASELINE,NERVECUT)
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

for i=1:length(BASELINE.rms_corr)
	npoints=size(BASELINE.spikes_corr{i},2)-1;
	baseline_data.spikes=[baseline_data.spikes BASELINE.spikes_corr{i}(1,2:end)];
	baseline_data.rms=[baseline_data.rms BASELINE.rms_corr{i}(1,2:end)];
	baseline_data.days_since=[baseline_data.days_since BASELINE.days_since{i}(2:end)];
	baseline_data.birdid=[baseline_data.birdid ones(1,npoints)*i];
end

baseline_data

artifacts=baseline_data.spikes<=r_cutoff;
artifacts2=isnan(baseline_data.spikes);

to_del=find(artifacts|artifacts2)

baseline_data.spikes(to_del)=[];
baseline_data.days_since(to_del)=[];
baseline_data.birdid(to_del)=[];
baseline_data.rms(to_del)=[];

% bin data, create baseline region and plot nervecut data on top using dot and error plot 

baseline_data
NERVECUT

win_steps=[-inf 20 40 60 inf ];
%win_steps(1)=-inf;
%win_steps(end)=inf;
[~,bins]=histc(baseline_data.days_since,win_steps);
condition={'rms','spikes'};

% collect p-value for each point

fig=figure();
for ii=1:2
	subplot(2,1,ii);
	for i=1:length(win_steps)-1

		cur_bin=[win_steps(i) win_steps(i+1)];
		cur_bindata=baseline_data.(condition{ii})(bins==i);
		cur_bindata(isnan(cur_bindata))=[];

		bin_mu=median(cur_bindata);
		%ci=bootci(1e3,{@median,cur_bindata});
		bin_conf(i)=prctile(cur_bindata,.5);
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
	xlabel('Days');
	set(gca,'YTick',[0:.2:1]);
end



%set(fig,'position',[200 200 240 480],'paperpositionmode','auto');
%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),[save_name],'eps,fig,png,pdf');
%save(fullfile(dirs.agg_dir,dirs.fig_dir,[ save_name '.mat']),'baseline_data');

