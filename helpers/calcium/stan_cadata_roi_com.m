function COM=stan_cadata_analyze(COORDS,IMAGE)
%
%
%
%
%
%

if ~isa(IMAGE,'double')
  IMAGE=double(IMAGE);
end

% image is a file, read in and convert to intensity...

nrois=length(COORDS);

% x,y pair (ROI format)

COM=zeros(nrois,2);
[height,width]=size(IMAGE);

for i=1:nrois

  mask=zeros(height,width);
  idx=sub2ind([height width],COORDS{i}(:,2),COORDS{i}(:,1));

  mask(idx)=1;
  newdata=IMAGE.*mask;

  % COM


  ind=1:height;
  tmp=newdata.*repmat(ind(:),[1 width]);
  COM(i,2)=sum(tmp(:))./sum(newdata(:));

  ind=1:width;
  tmp=newdata.*repmat(ind,[height 1]);
  COM(i,1)=sum(tmp(:))./sum(newdata(:));

end
