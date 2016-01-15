function [figs,stats]=stan_audio_overnight_sdi()
%%%
%
%
%

[options,dirs]=stan_preflight;

% load in each mat file, make spectral density image, etc.

% assume we have contours1

% get listing

low=.15;
high=.75;
saturation=1.4;

listing(1).name=fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis','lg318rblk_2012-08-23','mic_data.mat');
listing(2).name=fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis','lg318rblk_2012-08-24','mic_data.mat');

padding=[.2 .2];
freq_band=[1e3 7e3];
thinning=10;

[pathname,filename,ext]=fileparts(listing(1).name);
load(fullfile(pathname,'overnight_sdi_consensus.mat'),'sdi');
thresh=.5;

% chop in half, take day/night night/day

% day/night first day

ntrials1=size(sdi(1).consensus,3);

pool1=1:floor(ntrials1/2);
pool2=ntrials1-(floor(ntrials1/2)-1):ntrials1;

day1=mean(log(sdi(1).consensus(:,:,pool1))>=thresh,3);
night1=mean(log(sdi(1).consensus(:,:,pool2))>=thresh,3);

ntrials2=size(sdi(2).consensus,3);
pool21=1:floor(ntrials2/2);
pool22=ntrials2-(floor(ntrials2/2)-1):ntrials2;

day2=mean(log(sdi(2).consensus(:,:,pool21))>=thresh,3);

im1=markolab_im_overlay({night1,day1},'low',low,'high',high,'saturation',saturation,'mapping',{1,3});
im2=markolab_im_overlay({night1,day2},'low',low,'high',high,'saturation',saturation,'mapping',{1,3});

freq_idx(1)=min(find(sdi(1).f>=freq_band(1)));
freq_idx(2)=max(find(sdi(1).f<=freq_band(2)));

t_idx(1)=min(find(sdi(1).t>=padding(1)));
t_idx(2)=max(find(sdi(1).t<=sdi(1).t(end)-padding(2)));

day1score=day1(freq_idx(1):freq_idx(2),t_idx(1):t_idx(2));
day2score=day2(freq_idx(1):freq_idx(2),t_idx(1):t_idx(2));
night1score=night1(freq_idx(1):freq_idx(2),t_idx(1):t_idx(2));

stats.day1night1=corr(night1score(:),day1score(:),'type','pearson');
stats.day2night1=corr(night1score(:),day2score(:),'type','pearson');

% timecourse

[rows,columns]=size(day1score);
day1timecourse=zeros(1,columns);
day2timecourse=zeros(1,columns);

for i=1:columns
  tmp1=corrcoef(day1score(:,i),night1score(:,i));
  tmp2=corrcoef(day2score(:,i),night1score(:,i));
  day1timecourse(i)=tmp1(2,1);
  day2timecourse(i)=tmp2(2,1);
end

wins_x=[...
  .4916 .5195;...
  .5711 .5850;...
  .72 .7817];
wins_y=[...
  1.9 3.6;...
  1.9 3;
  4.25 6];

% zoom into each window save, add boxes manually?

m=7; % needs to be odd
n=4;
im1_idx=repmat([1:n-1],[floor(m/2) 1])+4.*repmat([0:m/2-1]',[1 n-1])
im2_idx=im1_idx+n*(floor(m/2))

ax=[];
figs.sdi=figure();
ax(1)=subplot(m,n,im1_idx(:));image(sdi(1).t,sdi(1).f/1e3,im1);
axis xy;box off;
ylim([0 9]);
hold on;

for i=1:size(wins_x,1)
  patch([wins_x(i,:) fliplr(wins_x(i,:))],...
    [wins_y(i,1).*ones(1,2) wins_y(i,2).*ones(1,2)],1,'facecolor','none','edgecolor','w')
end

set(gca,'xtick',[],'ytick',[])
ax(2)=subplot(m,n,im2_idx(:));image(sdi(2).t,sdi(2).f/1e3,im2);
axis xy;box off;
set(gca,'XTick',[],'YTick',[])
ylim([0 9]);
hold on;

for i=1:size(wins_x,1)
  patch([wins_x(i,:) fliplr(wins_x(i,:))],...
  [wins_y(i,1).*ones(1,2) wins_y(i,2).*ones(1,2)],1,'facecolor','none','edgecolor','w')
end

set(gca,'xtick',[])
ax(3)=subplot(m,n,(n*(m-1)+1):(n*m-1));
plot(sdi(1).t(t_idx(1):t_idx(2)),day1timecourse);
hold on;
plot(sdi(2).t(t_idx(1):t_idx(2)),day2timecourse,'r-');
ylim([.5 1.02]);
set(gca,'Ticklength',[0 0],'YTick',[.5 1],'FontSize',7)
box off;
linkaxes(ax,'x');
xlim([.2 .82])
subplot(m,n,[1:floor(m/2)]*n);
plot(night1score(1:thinning:end),day1score(1:thinning:end),'k.','markersize',5);
hold on;
plot([0 1],[0 1],'k-')
ylabel('Day 1 (px)');
set(gca,'TickLength',[0 0],'xtick',[0 1],'ytick',[0 1],'FontSize',7)

subplot(m,n,[1:floor(m/2)]*n+n*floor(m/2));
plot(night1score(1:thinning:end),day2score(1:thinning:end),'k.','markersize',5);
hold on;
plot([0 1],[0 1],'k-')
ylabel('Day 2 (px)');
xlabel('Night 1 (px)');
set(gca,'TickLength',[0 0],'xtick',[0 1],'ytick',[0 1],'FontSize',7)

for i=1:size(wins_x,1)

  ax=[];
  figs.([sprintf('zoomfig_%i',i)])=figure();

  ax(1)=subplot(2,1,1);
  image(sdi(1).t,sdi(1).f/1e3,im1);
  axis xy;box off;
  ylim([0 9]);

  ax(2)=subplot(2,1,2);
  image(sdi(2).t,sdi(2).f/1e3,im2);
  axis xy;box off;
  ylim([0 9]);
  linkaxes(ax,'xy');
  axis([wins_x(i,:) wins_y(i,:)])

end
