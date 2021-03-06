%
%
%
%
% batch analysis of all calcium imaging birds

[options,dirs]=stan_preflight;
motif_select=2;
ext='con';
listing=dir(fullfile(dirs.agg_dir,dirs.ca_dir,ext,'*.mat'));
maxlag=.1;
compare_day=1;

for i=1:3

  % use only the selected motif
  disp([listing(i).name]);
  cur=[];
  cur=load(fullfile(dirs.agg_dir,dirs.ca_dir,ext,listing(i).name),...
    'roi_data','roi_motifs','roi_params','roi_dates');

  lag_idx=zeros(1,length(cur.roi_data));
  len=cellfun(@length,cur.roi_data);
  to_del=len==0;

  if length(strfind(listing(i).name,'lw76'))>0
    tmp_motif_select=1;
    lag_corr=0;
  else
    tmp_motif_select=motif_select;
    lag_corr=1;
  end

  % pretty much all of the pads are incorrect, fun...

  if length(strfind(listing(i).name,'lny13'))>0
    cur.roi_params(1).padding=[.7 .677];
    %to_del(5)=true; % something strange happened on day 5
  end

  if length(strfind(listing(i).name,'lny18'))>0
    cur.roi_params(1).padding=[.5 .95];
  end

  if length(strfind(listing(i).name,'lny54rb'))>0
    cur.roi_params(1).padding=[.3 .85];
  end

  for j=1:length(cur.roi_data)
      cur.trial_idx{j}=1:size(cur.roi_data{j},3);
  end
  
  for j=1:length(cur.roi_data)
    cur.roi_data{j}=cur.roi_data{j}(:,:,cur.roi_motifs{j}==tmp_motif_select);
    cur.roi_dates{j}=cur.roi_dates{j}(cur.roi_motifs{j}==tmp_motif_select);
    
    lag_idx(j)=round(min(cur.roi_dates{j})-min(cur.roi_dates{1}));
  end

  cur.roi_data(to_del)=[];
  cur.roi_motifs(to_del)=[];
  cur.roi_params(to_del)=[];
  cur.roi_dates(to_del)=[];
  cur.trial_idx(to_del)=[];
  lag_idx(to_del)=[];

  for j=1:length(cur.roi_data)

    % re-sort by day/night

    [~,idx]=sort(cur.roi_dates{j},'ascend');
    cur.roi_data{j}=cur.roi_data{j}(:,:,idx);
    cur.roi_dates{j}=cur.roi_dates{j}(idx);
    cur.trial_idx{j}=cur.trial_idx{j}(idx);

  end

  % easy to assign lag indices, round off day difference between two datenumbers

  % insert global shift if we need to


  pad_smps=round(cur.roi_params(1).padding*cur.roi_params(1).fs);
  template=zscore(mean(cur.roi_data{compare_day}(pad_smps(1):end-pad_smps(2),:,:),3));

  ndays=length(cur.roi_data);
  nrois=size(cur.roi_data{1},2);
  global_shift=nan(1,ndays);

  for j=1:ndays

    roi_shifts=nan(1,nrois);
    mu=zscore(mean(cur.roi_data{j}(pad_smps(1):end-pad_smps(2),:,:),3));

    for k=1:nrois
      [r,lags]=xcorr(template(:,k),mu(:,k));
      [~,idx]=max(r);
      roi_shifts(k)=lags(idx);
    end

    global_shift(j)=round(median(roi_shifts));

  end

  global_shift

  % take max shift, we'll need to crop out that data...

  crop=max(global_shift)

  % crop cuts in from left and right, adjust pads if necessary...

  if crop>0
    for j=1:ndays
      cur.roi_data{j}=circshift(cur.roi_data{j},global_shift(j),1);
      cur.roi_data{j}=cur.roi_data{j}(crop:end-crop,:,:);
    end
  end

  pad_smps=pad_smps-crop;

  all_ca=cat(3,cur.roi_data{:});
  all_dates=cat(2,cur.roi_dates{:});
  all_trial_idx=cat(2,cur.trial_idx{:});
  ntrials=size(all_ca,3);
  nrois=size(all_ca,2);

  max_smps=round(maxlag*cur.roi_params(1).fs);

  corrmat{i}=zeros(ntrials,ntrials,nrois);
  timemat{i}=zeros(ntrials,ntrials);

  for k=1:nrois

    ca_data=zscore(squeeze(all_ca(pad_smps(1):end-pad_smps(2),k,:)));
    k
    %corrmat{i}(:,:,k)=1-squareform(pdist(ca_data','correlation'));
    for l=1:ntrials
      for m=1:ntrials
        corrmat{i}(l,m,k)=max(xcorr(ca_data(:,l),ca_data(:,m),max_smps,'coeff'));
        if k==1
          timemat{i}(l,m)=all_dates(l)-all_dates(m);
        end
      end
    end

  end

  clear cur;

end

% now for each bird take the upper triangle of the correlation matrix

all_ca=[];
all_time=[];
all_id=[];

for i=1:length(corrmat)
  idx=find(tril(true(size(corrmat{i}(:,:,1))),-1));
  for j=1
    tmp=median(corrmat{i},3);
    all_ca=[all_ca;tmp(idx)];
    tmp=timemat{i};
    all_time=[all_time;tmp(idx)];
    all_id=[all_id;ones(size(tmp(idx)))*i];
  end
end

% plot all points with smoothing line???


%save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-' ext '.mat']),'stats');
