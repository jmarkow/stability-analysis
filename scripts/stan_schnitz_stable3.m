%%% assumes the calcium data and the peaktime starts are all loaded in

[options,dirs]=stan_preflight;
%load(fullfile(dirs.agg_dir,dirs.ca_dir,'lw76.mat'));
%load(fullfile(dirs.agg_dir,dirs.datastore_dir,'mu_ca_timecourse.mat'));
idx=4;
clims=[.25 3];
clims2=[0 6];
motif_select=1;
%unstable_mat=cat(1,stats(idx).peak_stable{:});
% we can define stability using xcorr or peak timing...
scaling='s';

smoothing=0.1;
padding=[.25 .75];


for i=1:length(roi_data)
    roi_data2{i}=roi_data{i}(:,:,roi_motifs{i}==motif_select);
end

figs.schnitz_stable=figure();

% ave_mat=stan_cadata_sortmat(unstable_roi_data,'scaling','l','sort_day',1,'smoothing',0,'smooth_kernel','g',...
%     'padding',padding,'movie_fs',roi_params(1).fs,'fig_row',2,'fig_nrows',2,'realign',0,'maxlag',.1);
% stan_cadata_sortmat(stable_roi_data,'scaling','l','sort_day',1,'smoothing',0,'smooth_kernel','g',...
%     'padding',padding,'movie_fs',roi_params(1).fs,'fig_row',1,'fig_nrows',2,'realign',0,'maxlag',.1);


ave_mat=stan_cadata_sortmat(roi_data2,'scaling',scaling,'sort_day',1,'smoothing',smoothing,'smooth_kernel','g',...
    'padding',padding,'movie_fs',roi_params(1).fs,'fig_row',1,'fig_nrows',2,'realign',0,'maxlag',.1,'threshold',.25);

ax=findall(gcf,'type','axes')
for i=1:length(ax)
    caxis(ax(i),[clims]);
end

stan_cadata_sortmat(roi_data2,'scaling','d','sort_day',1,'smoothing',0,'smooth_kernel','g',...
    'padding',padding,'movie_fs',roi_params(1).fs,'fig_row',2,'fig_nrows',2,'realign',0,'maxlag',.1,'threshold',.25);

ax=findall(gcf,'type','axes')
for i=1:5
  caxis(ax(i),[clims2]);
end


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
h=line([xlimits(1) xlimits(1)+.2],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(end))
set(h,'clipping','off');
colormap(hot)
set(figs.schnitz_stable,'PaperPositionMode','auto','position',[200 200 280 160]);
%
% for i=1:length(ax)
%     caxis(ax(i),[clims]);
% end

figs.schnitz_colorbar=figure();
imagesc([0:1/1e3:1],[],[0:1/1e3:1]);
colormap(hot);
caxis([clims]);
set(gca,'XTick',[0 .5 1]);
set(gca,'YTick',[]);
set(figs.schnitz_colorbar,'PaperPositionMode','auto','position',[200 200 80 30]);

markolab_multi_fig_save(figs.schnitz_stable,pwd,['lw76_combined_stable' sprintf('%g-%g_scale-%s_smoothing-%g',clims(1),clims(2),scaling,smoothing) ],'eps,png,fig','renderer','painters');
%markolab_multi_fig_save(figs.schnitz_colorbar,pwd,['lw76_' sprintf('%g-%g',clims(1),clims(2)) '_colorbar' ],'eps,png,fig','renderer','painters');
