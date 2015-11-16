function stan_cadata_drift_analyze_all()
%
%
%
%
% batch analysis of all calcium imaging birds

maxlag=.033;
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

[stats(1).rmat_mu stats(1).pmat]=...
  stan_cadata_drift_analyze(cadata.data(1).roi_data,'padding',[1 1],'movie_fs',cadata.fs(1),...
  'lag_corr',1,'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag);

[stats(2).rmat_mu stats(2).pmat]=...
  stan_cadata_drift_analyze(cadata.data(2).roi_data,'padding',[.25 .22],'movie_fs',cadata.fs(2),...
  'lag_corr',1,'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag);

[stats(3).rmat_mu stats(3).pmat]=...
  stan_cadata_drift_analyze(cadata.data(3).roi_data,'padding',[0 .18],'movie_fs',cadata.fs(3),...
  'lag_corr',1,'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag);

save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats.mat']),'stats','cadata');
