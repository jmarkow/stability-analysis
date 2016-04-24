function MOVIE=stan_cadata_crop_movie(FRAMES,X,Y,PAD)
% given ROI x and y points, crop out movie
%
%
%

if nargin<4
  PAD=50;
end

bbox=[min(X) max(X);min(Y) max(Y)];
bbox_pad=round(bbox+[-PAD PAD;-PAD PAD]);

[height,width,nframes]=size(FRAMES);

if any(bbox_pad(1,:)<0|bbox_pad(1,:)>width)
  disp('Out of movie bounds...');
  MOVIE=[];
  return;
end

if any(bbox_pad(2,:)<0|bbox_pad(2,:)>height)
  disp('Out of movie bounds...');
  MOVIE=[];
  return;
end

% if we clear the checks extract here...

MOVIE=FRAMES(bbox_pad(2,1):bbox_pad(2,2),bbox_pad(1,1):bbox_pad(1,2),:);
