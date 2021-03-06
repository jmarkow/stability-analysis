function roi_stats=stan_collect_roistats(DIR)
%
%
%
%
%

filter='ROI_DATA.mat';
listing=dir(fullfile(DIR,'*.mat'));
counter=1;
for i=1:length(listing)
  listing(i).name
  if any(strfind(listing(i).name,filter))
    load(fullfile(DIR,listing(i).name),'ROI');
    roi_stats(counter).filename=listing(i).name;
    roi_stats(counter).coords=ROI.coordinates;
    roi_stats(counter).bw_stats=ROI.stats;
    roi_stats(counter).weighted_com=stan_cadata_roi_com(ROI.coordinates,ROI.reference_image);
    roi_stats(counter).max_proj=ROI.reference_image;
    counter=counter+1;
  end
end

%save('roi_spatial_stats.mat','roi_stats');
