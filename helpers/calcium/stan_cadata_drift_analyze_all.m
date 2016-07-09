function stan_cadata_drift_analyze_all(ext)
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

if nargin<1
    ext='lib';
end

listing=dir(fullfile(dirs.agg_dir,dirs.ca_dir,ext,'*.mat'));
maxlag=.1;
move_thresh=20;

parfor i=1:length(listing)

  % use only the selected motif
  disp([listing(i).name]);
  
  cur=load(fullfile(dirs.agg_dir,dirs.ca_dir,ext,listing(i).name),...
    'roi_data','roi_motifs','roi_params','roi_dates');
  birdid=regexprep(listing(i).name,'\-\w+\.mat','');
  
  cur_roi=load(fullfile(dirs.agg_dir,dirs.ca_dir,['roi_' ext],[ birdid '_roi_stats-' ext '.mat']));
  
  % strip out ROIs that move too much, save indices of the ROIs we keep
  
  lag_idx=zeros(1,length(cur.roi_data));
  len=cellfun(@length,cur.roi_data);
  to_del=len==0;

  if length(strfind(listing(i).name,'lw76'))>0
    tmp_motif_select=1;
    lag_corr=0;
    pxtomicron=.6250;
  else
    tmp_motif_select=motif_select;
    lag_corr=1;
    pxtomicron=1.1;
  end

  % pretty much all of the pads are incorrect, fun...

  if length(strfind(listing(i).name,'lny13'))>0
    cur.roi_params(1).padding=[.7 .677];
    %to_del(5)=true; % something strange happened on day 5
  end

  if length(strfind(listing(i).name,'lny18'))>0
    cur.roi_params(1).padding=[.5 .95];
  end

  if length(strfind(listing(i).name,'lny54rb'))>0
    cur.roi_params(1).padding=[.3 .85];
  end

  for j=1:length(cur.roi_data)
    cur.roi_data{j}=cur.roi_data{j}(:,:,cur.roi_motifs{j}==tmp_motif_select);
    cur.roi_dates{j}=cur.roi_dates{j}(cur.roi_motifs{j}==tmp_motif_select);
    lag_idx(j)=round(min(cur.roi_dates{j})-min(cur.roi_dates{1}));
  end

  cur.roi_data(to_del)=[];
  cur.roi_motifs(to_del)=[];
  cur.roi_params(to_del)=[];
  cur.roi_dates(to_del)=[];
  lag_idx(to_del)=[];

  for j=1:length(cur.roi_data)

    % re-sort by day/night

    [~,idx]=sort(cur.roi_dates{j},'ascend');
    cur.roi_data{j}=cur.roi_data{j}(:,:,idx);
    cur.roi_dates{j}=cur.roi_dates{j}(idx);

  end
  
  nrois=size(cur.roi_data{j},2);
  ndays=length(cur_roi.roi_stats);
  roi_movement=zeros(ndays-1,nrois);
  tmp=cellfun(@range,cur_roi.roi_stats(1).coords,'uniformoutput',false);
  thresholds=(min(cat(1,tmp{:}),[],2))'*pxtomicron;
  
  for j=1:ndays-1
    roi_movement(j,:)=max(abs(cur_roi.roi_stats(j).weighted_com-cur_roi.roi_stats(j+1).weighted_com)')*pxtomicron;
  end
  
  cur.roi_include=~any(roi_movement>repmat(thresholds,[ndays-1 1]));
  
  for j=1:length(cur.roi_data)
      cur.roi_data{j}=cur.roi_data{j}(:,cur.roi_include,:);
  end
  
  % easy to assign lag indices, round off day difference between two datenumbers

  [stats(i).rmat_mu stats(i).pmat stats(i).vmat]=stan_cadata_drift_analyze(...
    cur.roi_data,lag_idx,'padding',cur.roi_params(1).padding,...
    'movie_fs',cur.roi_params(1).fs,'lag_corr',lag_corr,...
    'realign',0,'smoothing',0,'smooth_kernel','b','maxlag',maxlag,'nboots',1e4);

  stats(i).trial_dates=cur.roi_dates;
  stats(i).use_data=cur;
  stats(i).use_roi=cur_roi;
  stats(i).use_roi.pxtomicron=pxtomicron;

end

save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-' ext '.mat']),'stats');
