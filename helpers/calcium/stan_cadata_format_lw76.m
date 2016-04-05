function ROI_DAT=stan_cadata_format_lw76(ROI_DAT)
%
%
%
% convert to Will's format (rois byt trials)

for i=1:length(ROI_DAT)

  ntrials=length(ROI_DAT{i}.align_detrended);
  [nsamples,nrois]=size(ROI_DAT{i}.align_detrended{1});

  %

  ROI_DAT{i}.raw_dat=cell(nrois,ntrials);
  ROI_DAT{i}.raw_time=cell(nrois,ntrials);


  for j=1:ntrials

    timevec=ROI_DAT{i}.frame_idx{j};
    timevec=timevec-timevec(1);
    timevec=timevec./ROI_DAT{i}.fs;

    for k=1:nrois

      ROI_DAT{i}.raw_dat{k,j}=ROI_DAT{i}.align_detrended{j}(:,k)';
      ROI_DAT{i}.raw_time{k,j}=timevec;

    end
  end

end
