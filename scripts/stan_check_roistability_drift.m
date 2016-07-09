%%%

% assumes roi_stats, stats, and fig stats are loaded in
ext='con';
filewrite=true;


[opts,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.datastore_dir,['mu_ca_timecourse-' ext '.mat']))
%load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_roi_new-' ext '.mat']))
load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-' ext '.mat']),'stats');

stable_drift=[];
unstable_drift=[];

for i=[1:4]

    unstable_mat=~fig_stats.drift.unstable{i};
    %use_data=roi_stats(i).com_drift;

    roi_include=stats(i).use_data.roi_include;

    %load(roi_files(i).name);
    [~,nrois]=size(unstable_mat);
    ndays=length(stats(i).use_roi.roi_stats);
    use_data=zeros(ndays-1,nrois);

    for j=1:ndays-1
        use_data(j,:)=sqrt(sum((stats(i).use_roi.roi_stats(j).weighted_com(roi_include,:)-...
          stats(i).use_roi.roi_stats(j+1).weighted_com(roi_include,:)).^2,2))';
    end

    for j=1:size(use_data,1)
        roi_idx=find(unstable_mat(j+1,:)==0&unstable_mat(1,:)==1);
        unstable_drift=[unstable_drift;use_data(j,roi_idx)'];
        roi_idx=find(unstable_mat(j+1,:)==1&unstable_mat(1,:)==1);
        stable_drift=[stable_drift;use_data(j,roi_idx)'];
    end
end

% plot CDF w/ bootstrap?

[p,~,stats]=ranksum(unstable_drift,stable_drift,'tail','right');

if filewrite
  fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,['castability_v_roistability-' ext '.txt']),'w+');
  fprintf(fid,'CA stability v ROI stability: p=%e z=%g\n',p,stats.zval);
  fclose(fid);
end
