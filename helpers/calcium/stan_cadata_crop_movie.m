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

bbox_pad(bbox_pad<0)=1;

if bbox_pad(1,2)>width
    bbox_pad(1,2)=width;
end

if bbox_pad(2,2)>height
    bbox_pad(2,2)=height;
end

% if we clear the checks extract here...

MOVIE=FRAMES(bbox_pad(2,1):bbox_pad(2,2),bbox_pad(1,1):bbox_pad(1,2),:);
