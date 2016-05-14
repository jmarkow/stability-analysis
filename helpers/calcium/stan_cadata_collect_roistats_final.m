function roi_stats=stan_collect_roistats()
%
%
%
%
%
% assumes we're in the roi map directory

listing=dir(pwd);
[pathname,filename,~]=fileparts(pwd);

if strfind(filename,'NEW')
  ext='lib';
else
  ext='con';
end

for i=1:length(listing)

  if listing(i).isdir & listing(i).name(1)~='.'
    birdid=regexprep(lower(listing(i).name),'_[a-z]+','');
    roi_stats=stan_cadata_collect_roistats(fullfile(listing(i).name,'mat'));


    save([birdid '_roi_stats-' ext '.mat'],'roi_stats');
  end

end
