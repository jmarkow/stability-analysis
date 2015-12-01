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

% sdi=[];
% for i=1:length(listing)
%
% 	disp([listing(i).name]);
% 	load(fullfile(listing(i).name),'agg_audio');
%
% 	[consensus f t]=acontour(agg_audio.data(:,1),agg_audio.fs);
% 	% compute sdi, relevant quantities, store
%
% 	ntrials=size(agg_audio.data,2);
% 	[rows,columns]=size(consensus);
%
% 	consensus=zeros(rows,columns,ntrials,'single');
%
% 	parfor j=1:ntrials
% 		disp([num2str(j) ' of ' num2str(ntrials)]);
% 		[consensus(:,:,j)]=acontour(agg_audio.data(:,j),agg_audio.fs);
% 	end
%
% 	sdi(i).consensus=consensus;
% 	sdi(i).f=f;
% 	sdi(i).t=t;
%
% 	clearvars consensus f t agg_audio;
%
% end
%
% save(fullfile(pathname,'overnight_sdi_consensus.mat'),'sdi','-v7.3');

[pathname,filename,ext]=fileparts(listing(1).name);
load(fullfile(pathname,'overnight_sdi_consensus.mat'),'sdi');
thresh=.9;

% chop in half, take day/night night/day

% day/night first day

ntrials1=size(sdi(1).consensus,3);

pool1=1:floor(ntrials1/2);
pool2=ntrials1-(floor(ntrials1/2)-1):ntrials1;

day1=mean(log(sdi(1).consensus(:,:,pool1))>=thresh,3);
day1=(day1-min(day1(:)))./(max(day1(:))-min(day1(:)));

day1_sc=min(day1,high);
day1_sc=max(day1_sc,low);
day1_sc=day1_sc-low;
day1_sc=day1_sc./(high-low);

night1=mean(log(sdi(1).consensus(:,:,pool2))>=thresh,3);
night1=(night1-min(night1(:)))./(max(night1(:))-min(night1(:)));

night1_sc=min(night1,high);
night1_sc=max(night1_sc,low);
night1_sc=night1_sc-low;
night1_sc=night1_sc./(high-low);

ntrials2=size(sdi(2).consensus,3);
pool21=1:floor(ntrials2/2);
pool22=ntrials2-(floor(ntrials2/2)-1):ntrials2;

day2=mean(log(sdi(2).consensus(:,:,pool21))>=thresh,3);
day2=(day2-min(day2(:)))./(max(day2(:))-min(day2(:)));

day2_sc=min(day2,high);
day2_sc=max(day2_sc,low);
day2_sc=day2_sc-low;
day2_sc=day2_sc./(high-low);

[rows,columns]=size(day1);
im1=zeros(rows,columns,3);
im1(:,:,1)=night1_sc;
im1(:,:,3)=day1_sc;
im1=im1.*saturation;
im1=imfilter(im1,fspecial('disk',2));

im2=zeros(size(im1));
im2(:,:,1)=night1_sc;
im2(:,:,3)=day2_sc;
im2=im2.*saturation;
im2=imfilter(im2,fspecial('disk',2));

freq_idx(1)=min(find(sdi(1).f>=freq_band(1)));
freq_idx(2)=max(find(sdi(1).f<=freq_band(2)));

t_idx(1)=min(find(sdi(1).t>=padding(1)));
t_idx(2)=max(find(sdi(1).t<=sdi(1).t(end)-padding(2)));

day1score=day1(freq_idx(1):freq_idx(2),t_idx(1):t_idx(2));
day2score=day2(freq_idx(1):freq_idx(2),t_idx(1):t_idx(2));
night1score=night1(freq_idx(1):freq_idx(2),t_idx(1):t_idx(2));

ax=[];
figs.scatter=figure();
ax(1)=subplot(2,1,1);
plot(night1score(1:thinning:end),day1score(1:thinning:end),'co');
hold on;
plot([0 1],[0 1],'k-')
ylabel('Day 1 (pixel value)');
set(gca,'TickLength',[0 0],'xtick',[0 1],'ytick',[0 1])

ax(2)=subplot(2,1,2);
plot(night1score(1:thinning:end),day2score(1:thinning:end),'co');
hold on;
plot([0 1],[0 1],'k-')
xlabel('Night 1')
ylabel('Day 2');
set(gca,'TickLength',[0 0],'xtick',[0 1],'ytick',[0 1])

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
  .6441 .7817];
wins_y=[...
  1.9 3.6;...
  1.9 3;
  4.25 6];

% zoom into each window save, add boxes manually?

ax=[];
figs.sdi=figure();
ax(1)=subplot(7,1,1:3);image(sdi(1).t,sdi(1).f/1e3,im1);
axis xy;box off;
ylim([0 9]);
hold on;

for i=1:size(wins_x,1)
  patch([wins_x(i,:) fliplr(wins_x(i,:))],...
    [wins_y(i,1).*ones(1,2) wins_y(i,2).*ones(1,2)],1,'facecolor','none','edgecolor','w')
end

set(gca,'xtick',[])
ax(2)=subplot(7,1,4:6);image(sdi(2).t,sdi(2).f/1e3,im2);
axis xy;box off;
ylim([0 9]);
hold on;

for i=1:size(wins_x,1)
  patch([wins_x(i,:) fliplr(wins_x(i,:))],...
  [wins_y(i,1).*ones(1,2) wins_y(i,2).*ones(1,2)],1,'facecolor','none','edgecolor','w')
end

set(gca,'xtick',[])
ax(3)=subplot(7,1,7);
plot(sdi(1).t(t_idx(1):t_idx(2)),day1timecourse);
hold on;
plot(sdi(2).t(t_idx(1):t_idx(2)),day2timecourse);
ylim([.5 1]);
set(gca,'Ticklength',[0 0],'YTick',[.5 1])
box off;
linkaxes(ax,'x');

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
