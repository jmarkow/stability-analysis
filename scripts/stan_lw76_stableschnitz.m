%% get filenames from new .mat file and data from old "stable"

ndays=length(ROI_dat);

roi_data=cell(1,ndays);
roi_motifs=cell(1,ndays);
roi_dates=cell(1,ndays);
roi_filenames=cell(1,ndays);
roi_times=cell(1,ndays);
roi_params=[];

for i=1:length(ROI_dat)

  ntrials=length(ROI_dat{i}.align_detrended);

  roi_data{i}=cat(3,ROI_dat{i}.align_detrended{:});
  roi_motifs{i}=zeros(1,ntrials);
  roi_dates{i}=zeros(1,ntrials);

  roi_filenames{i}=A{i}.Filenames;
  roi_times{i}=A{i}.Align_frame_idx{1}./24.414e3;

  roi_params(i).fs=22;
  roi_params(i).padding=[.25 .75];

  for j=1:ntrials
      tmp=A{i}.Filenames{j};
      tokens=regexp(tmp,'\_(\d+)\_','tokens');
      roi_dates{i}(j)=datenum([ tokens{1}{1} tokens{2}{1} ],'yyyymmddHHMMSS');
      tokens2=regexp(tmp,'\_(\d+)$','tokens');
      roi_motifs{i}(j)=str2num(tokens2{1}{1});
  end
end

%