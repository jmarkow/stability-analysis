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
%
% cadata.fs=[22 100 100];
%
% % analyze data, leave out pads
%
% [stats(1).peak_stable stats(1).peak_ispeak]=stan_cadata_drift_analyze_peaktime(cadata.data(1).roi_data,...
%   'padding',[1 1],'movie_fs',cadata.fs(1),'realign',0,'smoothing',0,'dist_thresh',.1);
% [stats(2).peak_stable stats(2).peak_ispeak]=stan_cadata_drift_analyze_peaktime(cadata.data(2).roi_data,...
%   'padding',[.25 .22],'movie_fs',cadata.fs(2),'realign',0,'smoothing',0,'dist_thresh',.1);
% [stats(3).peak_stable stats(3).peak_ispeak]=stan_cadata_drift_analyze_peaktime(cadata.data(3).roi_data,...
%   'padding',[0 .18],'movie_fs',cadata.fs(3),'realign',0,'smoothing',0,'dist_thresh',.1);
%
% save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_peaktime.mat']),'stats');


[options,dirs]=stan_preflight;
motif_select=2;
listing=dir(fullfile(dirs.agg_dir,dirs.ca_dir,'*.mat'));
maxlag=.02;

for i=1:length(listing)

  % use only the selected motif
  disp([listing(i).name]);
  cur=load(fullfile(dirs.agg_dir,dirs.ca_dir,listing(i).name),...
    'roi_data','roi_motifs','roi_params','roi_dates');

  lag_idx=zeros(1,length(cur.roi_data));

  len=cellfun(@length,cur.roi_data);

  to_del=len==0;
  cur.roi_data(to_del)=[];
  cur.roi_motifs(to_del)=[];
  cur.roi_params(to_del)=[];
  cur.roi_dates(to_del)=[];

  if strcmp(listing(i).name,'lw76.mat')
    tmp_motif_select=1;
    lag_corr=0;
  else
    tmp_motif_select=motif_select;
    lag_corr=0;
  end

  if strcmp(listing(i).name,'lny13.mat')
    cur.roi_params(1).padding=[.7 .5];
  end

  if strcmp(listing(i).name,'lny18.mat')
    cur.roi_params(1).padding=[.7 .5];
  end

  for j=1:length(cur.roi_data)
    cur.roi_data{j}=cur.roi_data{j}(:,:,cur.roi_motifs{j}==tmp_motif_select);
    lag_idx(j)=round(cur.roi_dates{j}(1)-cur.roi_dates{1}(1));
  end

  % easy to assign lag indices, round off day difference between two datenumbers

  [stats(i).peak_stable stats(i).peak_ispeak]=stan_cadata_drift_analyze_peaktime(...
    cur.roi_data,lag_idx,'padding',cur.roi_params(1).padding,...
    'movie_fs',cur.roi_params(1).fs,'smoothing',0,'realign',1,'dist_thresh',.1,'maxlag',maxlag);

end

save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_peaktime_new.mat']),'stats');
