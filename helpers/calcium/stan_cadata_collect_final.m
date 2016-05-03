% simply collects the calcium data from Will's folder, formats and saves to
% a more convenient format for further analysis...

[opts,dirs]=stan_preflight;

listing=robofinch_dir_recurse(pwd,'ROI_data_cleansed.mat');

% all except lw76

for i=1:length(listing)
  load(listing(i).name);
  tokens=regexp(listing(i).name,filesep,'split');
  [roi_data,roi_dates,roi_times,roi_motifs,roi_filenames,roi_params]=...
              stan_cadata_collect_freedomscope_v2_cell(ROI_data_cleansed);
  save(fullfile(dirs.agg_dir,dirs.ca_dir,[lower(tokens{end-1}) '.mat']),'roi_data',...
    'roi_dates','roi_times','roi_motifs','roi_filenames','roi_params');
end
