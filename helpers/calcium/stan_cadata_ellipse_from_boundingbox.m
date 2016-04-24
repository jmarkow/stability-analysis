function [X,Y]=stan_cadata_ellipse_from_boundingbox(X,Y,NPOINTS)
% given x,y coordinates from a bounding box (say from ImageJ),
% get the corresponding ellipse

if nargin<3
  NPOINTS=50;
end

phi = linspace(0,2*pi,NPOINTS);

cosphi = cos(phi);
sinphi = sin(phi);

xbar=mean(X);
ybar=mean(Y);

xlength=range(X);
ylength=range(Y);

a=max(xlength,ylength)/2;
b=min(xlength,ylength)/2;

% what's the major axis, that's our orientation

if xlength>ylength
  theta=0;
else
  theta=pi/2;
end

R = [ cos(theta)   sin(theta)
-sin(theta)   cos(theta)];

xy = [a*cosphi; b*sinphi];
xy = R*xy;

X = xy(1,:) + xbar;
Y = xy(2,:) + ybar;
