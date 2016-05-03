function stan_cadata_drift_analyze_all()
%
%
%
%
% batch analysis of all calcium imaging birds

% maxlag=.0165;
% %maxlag=.1;
% [options,dirs]=stan_preflight;
%
% % bird 1
%
% cadata.data(1)=load(fullfile(dirs.agg_dir,dirs.ca_dir,'lw76.mat'),'roi_data'); % inscopix frame rate 22 Hz
%
% % bird 2
%
% cadata.data(2)=load(fullfile(dirs.agg_dir,dirs.ca_dir,'lny13.mat'),'roi_data');
%
% % bird 3
%
% cadata.data(3)=load(fullfile(dirs.agg_dir,dirs.ca_dir,'lny18.mat'),'roi_data');
% cadata.fs=[22 100 100];
%
% % analyze data, leave out pads
%
% [stats(1).rmat_mu stats(1).pmat]=...
%   stan_cadata_drift_analyze(cadata.data(1).roi_data,'padding',[1 1],'movie_fs',cadata.fs(1),...
%   'lag_corr',1,'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag);
%
% [stats(2).rmat_mu stats(2).pmat]=...
%   stan_cadata_drift_analyze(cadata.data(2).roi_data,'padding',[.25 .22],'movie_fs',cadata.fs(2),...
%   'lag_corr',1,'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag);
%
% [stats(3).rmat_mu stats(3).pmat]=...
%   stan_cadata_drift_analyze(cadata.data(3).roi_data,'padding',[0 .18],'movie_fs',cadata.fs(3),...
%   'lag_corr',1,'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag);
%
% save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats.mat']),'stats','cadata');

[options,dirs]=stan_preflight;
motif_select=2;
listing=dir(fullfile(dirs.agg_dir,dirs.ca_dir,'*.mat'));
maxlag=.03;

for i=1:length(listing)

  % use only the selected motif
  listing(i).name
  load(fullfile(dirs.agg_dir,dirs.ca_dir,listing(i).name),'roi_data','roi_motifs','roi_params','roi_dates');

  lag_idx=zeros(1,length(roi_data));

  for j=1:length(roi_data)
    roi_data{j}=roi_data{j}(:,:,roi_motifs{j}==motif_select);
    lag_idx(j)=round(roi_dates{j}(1)-roi_dates{1}(1));
  end

  % easy to assign lag indices, round off day difference between two datenumbers

  [stats(i).rmat_mu stats(i).pmat]=stan_cadata_drift_analyze(...
    roi_data,lag_idx,'padding',roi_params(1).padding,...
    'movie_fs',roi_params(1).fs,'lag_corr',1,...
    'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag,'nboots',1e2);

end

save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new.mat']),'stats');
