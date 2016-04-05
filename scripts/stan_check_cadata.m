% script for analyzing new calcium data

filename='testing';

%%

% only for lw76
roi_data=cell(1,length(ROI_dat));

for i=1:length(ROI_dat)
    roi_data{i}=cat(3,ROI_dat{i}.align_detrended{:});
end

% assumes data already loaded in

%%

[roi_data,roi_dates,roi_times]=...
 stan_cadata_collect_freedomscope_v2_cell(ROI_data_cleansed(1:3));

%%

figure();
stan_cadata_sortmat(roi_data,'scaling','l','sort_day',1,'smoothing',0,'smooth_kernel','g',...
    'padding',[.2 0.5],'movie_fs',22,'chk_day',1,'fig_row',1,'fig_nrows',2,'realign',0);
stan_cadata_sortmat(roi_data,'scaling','l','sort_day',5,'smoothing',0,'smooth_kernel','g',...
  'padding',[.2 0.5],'movie_fs',22,'chk_day',1,'fig_row',2,'fig_nrows',2,'realign',0);
load custom_colormaps;
colormap(calcium_contrast);