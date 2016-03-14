function stan_cadata_drift_analyze_all()
%
%
%
%
% batch analysis of all calcium imaging birds

maxlag=.0165;
%maxlag=.1;
[options,dirs]=stan_preflight;

% bird 1

cadata.data(1)=load(fullfile(dirs.agg_dir,dirs.ca_dir,'lw76.mat'),'roi_data'); % inscopix frame rate 22 Hz

% bird 2

cadata.data(2)=load(fullfile(dirs.agg_dir,dirs.ca_dir,'lny13.mat'),'roi_data');

% bird 3

cadata.data(3)=load(fullfile(dirs.agg_dir,dirs.ca_dir,'lny18.mat'),'roi_data');

cadata.fs=[22 100 100];

% analyze data, leave out pads

[stats(1).peak_stable]=stan_cadata_drift_analyze_peaktime(cadata.data(1).roi_data,...
  'padding',[1 1],'movie_fs',cadata.fs(1),'realign',0,'smoothing',0,'dist_thresh',.1);
[stats(2).peak_stable]=stan_cadata_drift_analyze_peaktime(cadata.data(2).roi_data,...
  'padding',[.25 .22],'movie_fs',cadata.fs(2),'realign',0,'smoothing',0,'dist_thresh',.1);
[stats(3).peak_stable]=stan_cadata_drift_analyze_peaktime(cadata.data(3).roi_data,...
  'padding',[0 .18],'movie_fs',cadata.fs(3),'realign',0,'smoothing',0,'dist_thresh',.1);

save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_peaktime.mat']),'stats');