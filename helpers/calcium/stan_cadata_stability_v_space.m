function dist=stan_cadata_stability_v_space()
%% load in the stability analysis

[options,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats.mat'),'stats');
peakstats=load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_stats_peaktime.mat'),'stats');
peakstats=peakstats.stats;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,'cadata_maps'),'roi_map');

%% load in the ROI spatial data

%% euclid distance between ROI centroids...

nboots=1000;

% scramble, repeat nboot times!

%dist.bootstrap=[];

vars={'stable','unstable','between'};

% set i to 1:3 once we get the rest of the data

for i=1:3

  % get centroid pairwise Euclidean distance

  for j=1:nboots
    dist(i).bootstrap{j}=[];
  end


  centroids=cat(1,roi_map(i).stats.Centroid);

  % compare within stable, within unstable and between two groups...

  unstable=peakstats(i).peak_stable(1,:)==1&any(peakstats(i).peak_stable(1:5,:)==0);
  stable=~unstable;

  unstable=find(unstable);
  stable=find(stable);

  centroids_dist.unstable=squareform(pdist(centroids(unstable,:),'Euclidean'));
  centroids_dist.stable=squareform(pdist(centroids(stable,:),'Euclidean'));
  centroids_dist.between=pdist2(centroids(stable,:),centroids(unstable,:),'Euclidean');

  for j=1:length(vars)
    dist(i).(vars{j})=[];
    mat=triu(ones(size(centroids_dist.(vars{j}))),1);
    dist(i).(vars{j})=[centroids_dist.(vars{j})(mat==1)];
  end

  [~,rndidx_all]=sort(rand(size(centroids,1),nboots));

  for j=1:nboots

    pool1=rndidx_all(1:length(unstable),j);
    pool2=rndidx_all(length(unstable)+1:end,j);

    % pairwise

    tmp=pdist2(centroids(pool1,:),centroids(pool2,:));
    mat=triu(ones(size(tmp)),1);
    dist(i).bootstrap{j}=tmp(mat==1);

  end

  % change this to bootstrap analysis

  %pairs.between=nchoosek([unstable(:);stable(:)],2);

end


save(fullfile(dirs.agg_dir,dirs.datastore_dir,'stability_v_space_peaktime.mat'),'dist');
