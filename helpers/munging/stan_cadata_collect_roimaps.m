function stan_cadata_collect_roimaps()
%% aggregate the ROI spatial data

[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.ca_dir,'ROI_MAPS','LW76_image_roi','roi_data_image_combined.mat'),'rois');
roi_map(1)=rois;
roi_map(1).type='';
roi_map(1).reference_image=[];

clear rois;

load(fullfile(dirs.agg_dir,dirs.ca_dir,'ROI_MAPS','LNY13_image_roi','roi_data_image.mat'),'ROI');
roi_map(2)=ROI;

clear ROI;

load(fullfile(dirs.agg_dir,dirs.ca_dir,'ROI_MAPS','LNY18_image_roi','roi_data_image.mat'),'ROI');
roi_map(3)=ROI;

clear ROI;

save(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_maps.mat'),'roi_map');
