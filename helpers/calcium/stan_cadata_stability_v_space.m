function dist=stan_cadata_stability_v_space()
%% load in the stability analysis

ext='lib';
[options,dirs]=stan_preflight;

load(fullfile(dirs.agg_dir,dirs.datastore_dir,['mu_ca_timecourse-' ext '.mat']))
%roi_dir=dir(fullfile(dirs.agg_dir,dirs.ca_dir,['roi_' ext],'*.mat'));
roi_files=robofinch_dir_recurse(fullfile(dirs.agg_dir,dirs.ca_dir,['roi_' ext]),'*.mat');

%% load in the ROI spatial data

%% euclid distance between ROI centroids...

nboots=1000;

% scramble, repeat nboot times!

%dist.bootstrap=[];

vars={'stable','unstable','between'};

% set i to 1:3 once we get the rest of the data

for i=1:4

  % get centroid pairwise Euclidean distance

  for j=1:nboots
    dist(i).bootstrap_stable{j}=[];
    dist(i).bootstrap_unstable{j}=[];
    dist(i).bootstrap_between{j}=[];
  end

  load(roi_files(i).name);
  tmp=cellfun(@mean,roi_stats(1).coords,'uniformoutput',false);
  centroids=cat(1,tmp{:});
  % compare within stable, within unstable and between two groups...

  unstable=fig_stats.drift.unstable{i}(1,:)==0&any(fig_stats.drift.unstable{i}(1:end,:)==1);
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

    tmp=pdist2(centroids(pool1,:),centroids(pool2,:),'Euclidean');
    mat=triu(ones(size(tmp)),1);
    dist(i).bootstrap_between{j}=tmp(mat==1);

    tmp=squareform(pdist(centroids(pool1,:),'Euclidean'));
    mat=triu(ones(size(tmp)),1);

    dist(i).bootstrap_unstable{j}=tmp(mat==1);

    tmp=squareform(pdist(centroids(pool2,:),'Euclidean'));
    mat=triu(ones(size(tmp)),1);

    dist(i).bootstrap_stable{j}=tmp(mat==1);

  end

  % change this to bootstrap analysis

  %pairs.between=nchoosek([unstable(:);stable(:)],2);

  dist(i).pval.between=1-mean(median(dist(i).between)<=cellfun(@median,dist(i).bootstrap_between));
  dist(i).pval.stable=1-mean(median(dist(i).stable)<=cellfun(@median,dist(i).bootstrap_stable));
  dist(i).pval.unstable=1-mean(median(dist(i).unstable)<=cellfun(@median,dist(i).bootstrap_unstable));

end


save(fullfile(dirs.agg_dir,dirs.datastore_dir,['stability_v_space_peaktime-' ext '.mat' ]),'dist');
