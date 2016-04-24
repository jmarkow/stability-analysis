function stan_cadata_crop_avis(CROP,FS,ROIS,PAD)
%
%
%
%

% given crop movies (nrois,ndays cell array), make montage movies

% 1: scale the same across days (within ROIs)
% 2:

if nargin<3
  ROIS=[];
end

climits=[.5 99.5]; %middle 95%

[nrois,ndays]=size(CROP);



%v.Quality=100;

for i=1:nrois

    % make the montage movie


    [height,width,nframes]=size(CROP{i,1});
    
    
    montage_movie=cat(2,CROP{i,:});

    % get climits

    tmp=prctile(montage_movie(:),climits);

    % that's our caxis

    if isempty(montage_movie)
        continue;
    end

    montage_movie=montage_movie-tmp(1);
    montage_movie(montage_movie>(tmp(2)-tmp(1)))=tmp(2)-tmp(1);
    montage_movie=montage_movie./(tmp(2)-tmp(1));
    montage_movie=uint8(montage_movie.*256);

    fprintf('Writing movie for roi %i\n',i);

    v=VideoWriter(sprintf('roi_%i',i),'Uncompressed AVI');
    v.FrameRate=FS;
    open(v);

    for j=1:nframes
        % insert borders???

        frame=ind2rgb(medfilt2(montage_movie(:,:,j),[11 11]),colormap('gray(256)'));

        if ~isempty(ROIS)

          newx=round(ROIS.coordinates{i}(:,1)-min(ROIS.coordinates{i}(:,1))+PAD);
          newy=round(ROIS.coordinates{i}(:,2)-min(ROIS.coordinates{i}(:,2))+PAD);

          for k=1:ndays
              % set new left and bottom edges
              
              cur_edge=(k-1)*(width);
              
              for l=1:length(newx)
                frame(newy(l),newx(l)+cur_edge,:)=[1 0 1];
              end
          end
        end

        %writeVideo(v,medfilt2(montage_movie(:,:,j),[11 11]));
        writeVideo(v,frame);

    end

    close(v);

end
