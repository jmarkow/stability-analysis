function [ROIS]=stan_cadata_crop_movie(CDATA)
% given ROI x and y points, crop out movie
%
%
%


% draw bboxes to make "ROIs" for cropped movies

clipping=[5 95];

% make figure, display projection, choose ROI, export as X Y coords and binary mask
clips=prctile(CDATA(:),clipping);
CDATA(CDATA<clips(1))=clips(1);
CDATA=CDATA-clips(1);
CDATA(CDATA>clips(2)-clips(1))=clips(2)-clips(1);

response='c';

select_fig=figure();
imagesc(CDATA);colormap(gray);colorbar();
caxis([0 clips(2)-clips(1)])

ROIS={};

while ishandle(select_fig)
  h=imrect;
  position=wait(h);
  delete(h);

  if ~isempty(position)
    ROIS{end+1}=[position(1) position(2);position(1)+position(3) position(2)+position(4)];
  end
end
