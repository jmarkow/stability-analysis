%
%
%

[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.sdi_dir,'sdi_data.mat'),'data');
load custom_colormaps;

%%%%% plot both, use nice colormap, etc. etc.
ax=[];

fig.sdi=figure();
ax(1)=subplot(2,1,1);
imagesc(data(1).t,data(1).f/1e3,data(1).sdi.im);
axis xy;
%colormap(tim_map2);
%caxis([0 .4]);
colorbar();

ax(2)=subplot(2,1,2);
imagesc(data(2).t,data(2).f/1e3,data(2).sdi.im);
axis xy;
%colormap(tim_map2);
%caxis([0 .4]);
colorbar();

linkaxes(ax,'xy');
