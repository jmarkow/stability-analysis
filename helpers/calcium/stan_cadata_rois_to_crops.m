function [CROP_MOVIE,ROIS]=stan_cadata_crop_movie(FRAMES,ROIS,MOTIF)
% given ROI x and y points, crop out movie
%
%
%

if nargin<3 | isempty(MOTIF)
  MOTIF=2;
end

nrois=length(ROIS.coordinates);
CROP_MOVIE=cell(nrois,length(FRAMES));
CROP_ROIS=cell(nrois,1);
% given a set of ROIS, make a bunch of montage movies motif-selective

if strcmpi(ROIS.type,'imagej')
  for i=1:nrois

    % convert bounding box to ellipse

    [newx,newy]=stan_cadata_ellipse_from_boundingbox(...
      ROIS.coordinates{i}(:,1),ROIS.coordinates{i}(:,2),200);

    ROIS.coordinates{i}=[newx(:) newy(:)];

  end
end

for i=1:nrois
  for j=1:length(FRAMES)
    CROP_MOVIE{i,j}=stan_cadata_crop_movie(FRAMES{j}{MOTIF},ROIS.coordinates{i}(:,1),ROIS.coordinates{i}(:,2),50);
  end
end
