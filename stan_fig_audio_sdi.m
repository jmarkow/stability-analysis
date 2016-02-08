function stan_fig_audio_sdi()

%
%
%

consensus=1;
threshold=5;

[options,dirs]=stan_preflight;

if consensus
  load(fullfile(dirs.agg_dir,dirs.sdi_dir,'sdi_data_consensus.mat'),'data');
else
  load(fullfile(dirs.agg_dir,dirs.sdi_dir,'sdi_data.mat'),'data');
end
load custom_colormaps;

%%%%% plot both, use nice colormap, etc. etc.
ax=[];

fig.sdi=figure();
ax(1)=subplot(2,1,1);

if consensus
  imagesc(data(1).t,data(1).f/1e3,mean(data(1).consensus>threshold,3));
else
  imagesc(data(1).t,data(1).f/1e3,data(1).sdi.im);
end

axis xy;
%colormap(tim_map2);
%caxis([0 .4]);
colorbar();

ax(2)=subplot(2,1,2);
if consensus
  imagesc(data(2).t,data(2).f/1e3,mean(data(2).consensus>threshold,3));
else
  imagesc(data(2).t,data(2).f/1e3,data(2).sdi.im);
end
axis xy;
%colormap(tim_map2);
%caxis([0 .4]);
colorbar();

linkaxes(ax,'xy');
