%%% assumes the calcium data and the peaktime starts are all loaded in

idx=4;
unstable_mat=cat(1,stats(4).peak_stable{:});
stable_rois=all(unstable_mat==1);
unstable_rois=unstable_mat(1,:)==1&any(unstable_mat(2:end,:)==0);

stable_roi_data=cell(size(roi_data));
unstable_roi_data=cell(size(roi_data));

for i=1:length(roi_data)
  stable_roi_data{i}=roi_data{i}(:,stable_rois,:);
  unstable_roi_data{i}=roi_data{i}(:,unstable_rois,:);
end

figure();

stan_cadata_sortmat(unstable_roi_data,'scaling','l','sort_day',1,'smoothing',0,'smooth_kernel','g',...
    'padding',[.25 0.75],'movie_fs',22,'fig_row',2,'fig_nrows',2,'realign',0);
stan_cadata_sortmat(stable_roi_data,'scaling','l','sort_day',1,'smoothing',0,'smooth_kernel','g',...
    'padding',[.25 0.75],'movie_fs',22,'fig_row',1,'fig_nrows',2,'realign',0);

colormap(jet);
ax=findall(gcf,'type','axes')
linkaxes(ax,'x')
%ylim([3.5 55.5]);
%xlim([.3 .72]);
ylimits=get(ax(5),'ylim');
set(ax(5),'YTick',ylimits,'YTickLabel',[ylimits-min(ylimits)]+[1 0],'FontSize',8);
yh=ylabel(ax(end),'Cell');
ylimits=get(ax(end),'ylim');
set(ax(end),'YTick',ylimits,'YTickLabel',[ylimits-min(ylimits)]+[1 0],'FontSize',8);
yh=ylabel(ax(end),'Cell');
%set(yh,'position',get(yh,'position')+[.3 0 0]);
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.1],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(end))
set(h,'clipping','off');


colormap(jet)