%%%

% assumes roi_stats, stats, and fig stats are loaded in
ext='lib';

[opts,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_peaktime_new-' ext '.mat']))
load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_roi_new-' ext '.mat']))
stable_drift=[];
unstable_drift=[];

for i=[1:4]
    unstable_mat=cat(1,stats(i).peak_stable{:});
    ispeak_mat=cat(1,stats(i).peak_ispeak{:});
    use_data=roi_stats(i).com_drift;
    for j=1:size(use_data,1)
        roi_idx=find(unstable_mat(j+1,:)==0&unstable_mat(1,:)==1&ispeak_mat(j+1,:)==1);
        unstable_drift=[unstable_drift;use_data(j,roi_idx)'];
        roi_idx=find(unstable_mat(j+1,:)==1&unstable_mat(1,:)==1&ispeak_mat(j+1,:)==1);
        stable_drift=[stable_drift;use_data(j,roi_idx)'];
    end
end


% plot CDF w/ bootstrap?
