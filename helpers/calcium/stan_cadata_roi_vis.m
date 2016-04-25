function stan_cadata_roi_vis(IMAGE,ROIS,IDX)
%
%
%
%
%
%

IMAGE=double(IMAGE);

imagesc(IMAGE);

caxis([0 256]);
colormap(gray);

hold on;

for i=1:length(ROIS)

  % use the sorting idx

  center=ROIS(IDX(i)).Centroid;

  % label with i, which corresponds to the appropriate idx in other plots

  text(center(1)+10,center(2)+10,num2str(i),'Color','Green','FontSize',15);
  plot(ROIS(IDX(i)).ConvexHull(:,1),ROIS(IDX(i)).ConvexHull(:,2),'r-','linewidth',1.5);

end

axis off;
