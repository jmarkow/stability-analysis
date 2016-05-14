%function stan_cadata_drift_analyze_all()
%
%
%
%
% batch analysis of all calcium imaging birds

[options,dirs]=stan_preflight;
ext='lib';
listing=dir(fullfile(dirs.agg_dir,dirs.ca_dir,['roi_' ext],'*.mat'));

for i=1:length(listing)

  % load in the centroids, check movement in centroids across days,
  % anything else???
  tmp=load(fullfile(dirs.agg_dir,dirs.ca_dir,['roi_' ext],listing(i).name))

  nrois=size(tmp.roi_stats(1).weighted_com,1);
  ndays=length(tmp.roi_stats);

  raw_stats{i}=tmp.roi_stats;
  roi_stats(i).com_drift=zeros(ndays-1,nrois);

  for j=1:ndays-1
    roi_stats(i).com_drift(j,:)=sqrt(sum((raw_stats{i}(j).weighted_com-raw_stats{i}(j+1).weighted_com).^2,2));
  end

  roi_stats(i).com=raw_stats{i}(1).weighted_com';
  % also perhaps collect the average com to check if

end

save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_roi_new-' ext '.mat']),'roi_stats');
