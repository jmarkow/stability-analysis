function [fig,fig_stats]=stan_cadata_drift_analyze_max(ext)
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
listing=dir(fullfile(dirs.agg_dir,dirs.ca_dir,ext,'*.mat'));
maxlag=.1;
compare_day=1;
filewrite=true;

all_ca=[];
all_ca_mu=[];
all_hrs=[];
all_id=[];
all_trials=[];
smoothing_tau=.1;
load(fullfile(dirs.agg_dir,dirs.datastore_dir,['cadata_stats_new-' ext '.mat']),'stats');

for i=1:3

  % use only the selected motif

  cur=stats(i).use_data;

  lag_idx=zeros(1,length(cur.roi_data));
  len=cellfun(@length,cur.roi_data);
  to_del=len==0;

  tmp_motif_select=motif_select;
  lag_corr=1;


  for j=1:length(cur.roi_data)
    lag_idx(j)=round(min(cur.roi_dates{j})-min(cur.roi_dates{1}));
  end

  for j=1:length(cur.roi_data)

    % re-sort by day/night

    [~,idx]=sort(cur.roi_dates{j},'ascend');
    cur.roi_data{j}=cur.roi_data{j}(:,:,idx);
    cur.roi_dates{j}=cur.roi_dates{j}(idx);
    %cur.trial_idx{j}=cur.trial_idx{j}(idx);

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
  nrois=size(cur.roi_data{1},2);

  % get max variability when binning by TIME

  % in hrs

  min_n=3;
  bin_size=1;
  bin_c=0:bin_size:24;

  for j=1:ndays

    hrs=(cur.roi_dates{j}-floor(cur.roi_dates{j}))*24;
    [~,bin_idx]=histc(hrs,bin_c);
    bins=unique(bin_idx);

    for k=1:length(bins)
      idx=bin_idx==bins(k);
      if sum(idx)>=min_n
        tmp=squeeze(max(cur.roi_data{j}));
        all_ca=[all_ca;std(tmp(:,idx),[],2)];
        all_ca_mu=[all_ca_mu;mean(tmp(:,idx),2)];
        all_hrs=[all_hrs;ones(size(tmp,1),1).*bins(k)];
        all_id=[all_id;ones(size(tmp,1),1)*i];
        all_trials=[all_trials;ones(size(tmp,1),1)*sum(hrs<=bins(k))];
        %all_trials=[all_trials;cur.trial_idx{j}(idx)'];
      end
    end
  end


end

all_ca=1./all_ca;

% uniq_hrs=unique(all_hrs)
% mu=[];
% for i=1:length(uniq_hrs)
%     mu(i)=mean(all_ca(all_hrs==uniq_hrs(i)));
%     ci(:,i)=bootci(10e3,{@mean,all_ca(all_hrs==uniq_hrs(i))},'type','cper');
% end
%
% fig.shaded_plot=figure();
% x=uniq_hrs;
% markolab_shadeplot(x,ci,[1 0 0],'none');

fig.camax_v_time=figure();
clip=1;
stan_plot_regress(all_hrs,all_ca,ones(size(all_id)),'clip',clip);
colormap(bone);
xlim([7 17])
%ylabel('Peak dF/F_0 SD')
%xlabel('Time (hr)')
set(gca,'TickLength',[0 0],'YTick',[0 clip],'FontSize',7);
title('')
ylim([-.15 clip+.11])
set(fig.camax_v_time,'position',[300 700 230 210]);
[fig_stats.timevcamax.r,fig_stats.timevcamax.p]=corr(all_hrs,all_ca,'type','spearman');
[fig_stats.timevcamax_partial.r,fig_stats.timevcamax_partial.p]=partialcorr(all_hrs,all_ca,all_ca_mu,'type','spearman');

if filewrite
  fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,['fig6_camaxvtime-' ext '.txt']),'w+');
  fprintf(fid,'Unit time v dF/F peak SD: r=%g p=%e\n',fig_stats.timevcamax.r,fig_stats.timevcamax.p);
  fprintf(fid,'Unit time v dF/F peak SD partial: r=%g p=%e\n',fig_stats.timevcamax_partial.r,fig_stats.timevcamax_partial.p);
  fclose(fid);
end
