function stan_fig1a(DATA,CONSENSUS,THRESHOLD)
% Generates Fig. 1a, time-frequency contours
%
% note: time-frequency contours are computed using the acontour library found at https://github.com/jmarkow/acontour

if nargin<3 | isempty(THRESHOLD)
  THRESHOLD=5; % threshold on power to binarize image
end

if nargin<2 | isempty(CONSENSUS)
  CONSENSUS=true;
end

ax=[];

fig=figure();
ax(1)=subplot(2,1,1);

if CONSENSUS
  imagesc(DATA(1).t,DATA(1).f/1e3,mean(DATA(1).consensus>THRESHOLD,3));
else
  imagesc(DATA(1).t,DATA(1).f/1e3,DATA(1).sdi.im);
end

axis xy;
colorbar();
caxis([.1 .8])

ax(2)=subplot(2,1,2);
if CONSENSUS
  imagesc(DATA(2).t,DATA(2).f/1e3,mean(DATA(2).consensus>THRESHOLD,3));
else
  imagesc(DATA(2).t,DATA(2).f/1e3,DATA(2).sdi.im);
end

axis xy;
colorbar();
caxis([.1 .8]);
colormap(hot);

linkaxes(ax,'xy');
xlim([.37 1]);
ylim([0 8]);
set(fig,'units','centimeters','position',[4 4 8.5725 8.2550],'paperpositionmode','auto');
