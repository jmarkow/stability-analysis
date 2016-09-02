function [fig,fig_stats]=stan_cadata_drift_analyze_pertime_ave(ext)
%
%
%
%
% batch analysis of all calcium imaging birds

[options,dirs]=stan_preflight;
motif_select=2;
if nargin<1
  ext='con';
end

maxlag=.1;
compare_day=1;
corrmat=[];
timemat=[];
trialmat=[];
filewrite=true;

% TODO make sure we account for actual lag!!!
load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-' ext '.mat']),'stats');

for i=1:3

  % use only the selected motif

  cur=stats(i).use_data;
  lag_idx=zeros(1,length(cur.roi_data));

  if i==4
    tmp_motif_select=1;
    lag_corr=0;
  else
    tmp_motif_select=motif_select;
    lag_corr=1;
  end

  for j=1:length(cur.roi_data)
      cur.trial_idx{j}=1:size(cur.roi_data{j},3);
  end

  for j=1:length(cur.roi_data)
    lag_idx(j)=round(min(cur.roi_dates{j})-min(cur.roi_dates{1}));
  end

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
  all_hrs=hour(all_dates);

  nrois=size(all_ca,2);

  max_smps=round(maxlag*cur.roi_params(1).fs);

  % normalize by averaging peak value in traces nearby in time (bin first)

  for j=1:length(cur.roi_data)

    hrs=hour(cur.roi_dates{j});
    bins=[min(hrs):2:max(hrs)];
    [~,binidx]=histc(hrs,bins);
    ntrials=size(cur.roi_data{j},3);

    for k=1:nrois
      for l=1:ntrials
        cur_hrs=hrs(l);
        hr_diffs=abs(cur_hrs-hrs);
        [~,idx]=sort(hr_diffs,'ascend');

        % take the first 5 (or something)

        peak_norm=median(max(cur.roi_data{j}(pad_smps(1):end-pad_smps(2),k,idx(1:min(10,ntrials)))));
        cur.roi_data{j}(:,k,l)=cur.roi_data{j}(:,k,l)./peak_norm;

      end
    end
  end

  for j=2:length(cur.roi_data)
    for k=1:nrois
      template=squeeze(mean((cur.roi_data{j-1}(pad_smps(1):end-pad_smps(2),k,:)),3));
      template=template-mean(template);

      for l=1:size(cur.roi_data{j},3)

        ca_data=(cur.roi_data{j}(pad_smps(1):end-pad_smps(2),k,l));
        ca_data=ca_data-mean(ca_data);

        lag_idx=round(min(cur.roi_dates{j-1})-min(cur.roi_dates{j}));

        if abs(lag_idx)>1
          continue;
        end

        corrmat{i,j-1}(k,l)=max(xcorr(template(:),ca_data(:),max_smps,'coeff'));
        %corrmat{i,j-1}(k,l)=max(xcorr(template(:),ca_data(:),max_smps,'coeff'));

        %[r,lags]=xcorr(template(:),ca_data(:),max_smps,'coeff');
        %[~,loc]=max(r);
        % scan over max_smps in both directions

        % shifts=-max_smps:max_smps;
        % dist=nan(1,length(shifts));
        %
        % for m=1:length(shifts)
        %   %dist(m)=sqrt(mean((template(:)-circshift(ca_data(:)-mean(ca_data(:)),[shifts(m) 1])).^2));
        %
        %   if shifts(m)>0
        %     shifted_copy=ca_data(1:end-(shifts(m)-1));
        %     template_copy=template(shifts(m):end);
        %   elseif shifts(m)<0
        %     shifted_copy=ca_data(-shifts(m):end);
        %     template_copy=template(1:end-(-shifts(m)-1));
        %   else
        %     shifted_copy=ca_data(:);
        %     template_copy=template(:);
        %   end
        %
        %
        %   dist(m)=mean(abs(template_copy(:)-shifted_copy(:)));
        %
        % end
        %
        % corrmat{i,j-1}(k,l)=min(dist);
        timemat{i,j-1}(k,l)=cur.roi_dates{j}(l)-floor(cur.roi_dates{j}(l));
        trialmat{i,j-1}(k,l)=cur.trial_idx{j}(l); % trial idx for the day
      end
    end
  end


  clear cur;

end

% now for each bird take the upper triangle of the correlation matrix
%%

all_ca=[];
all_time=[];
all_id=[];
all_trial=[];

for i=1:size(corrmat,1)
  for j=1:size(corrmat,2)
    if ~isempty(corrmat{i,j})
      all_ca=[all_ca;mean(corrmat{i,j})'];
      all_time=[all_time;timemat{i,j}(1,:)'];
      all_trial=[all_trial;trialmat{i,j}(1,:)'];
      all_id=[all_id;ones(size(mean(corrmat{i,j})'))*i];
    end
  end
end

% plot all points with smoothing line???

[~,idx]=sort(all_time);
fig.nminus1_regress=figure();
stan_plot_regress(all_time*24,all_ca,ones(size(all_id)));
colormap(bone);
xlim([7 17])
%ylabel('Correlation to n-1 ave (r)')
%xlabel('Time (hr)')
set(gca,'TickLength',[0 0],'YTick',[.2 .9],'FontSize',7);
title('')
ylim([.2 .9])
set(fig.nminus1_regress,'position',[300 700 230 210]);
[fig_stats.timevca.r,fig_stats.timevca.p]=corr(all_time,all_ca,'type','spearman');

if filewrite
  fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,['fig6_unittimevca-' ext '.txt']),'w+');
  fprintf(fid,'Unit time v n-1 ave corr: r=%g p=%e, n=%i',fig_stats.timevca.r,fig_stats.timevca.p,length(all_time));
  fclose(fid);
end

%save(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-' ext '.mat']),'stats');
