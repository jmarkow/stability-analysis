% simply collects the calcium data from Will's folder, formats and saves to
% a more convenient format for further analysis...

[opts,dirs]=stan_preflight;
ext='con';

listing=robofinch_dir_recurse(pwd,['ROI_data_cleansed-' ext '.mat']);

% all except lw76

for i=1:length(listing)-1
  load(listing(i).name);
  tokens=regexp(listing(i).name,filesep,'split');
  [roi_data,roi_dates,roi_times,roi_motifs,roi_filenames,roi_params]=...
              stan_cadata_collect_freedomscope_v2_cell(ROI_data_cleansed);
  save(fullfile(dirs.agg_dir,dirs.ca_dir,[lower(tokens{end-1}) '-' ext '.mat']),'roi_data',...
    'roi_dates','roi_times','roi_motifs','roi_filenames','roi_params');
end

% now collect lw76

%%
clearvars roi_*;
load(listing(end).name);

% parse motif numbers and dates

ndays=length(ROI_data_cleansed);

roi_data=cell(1,ndays);
roi_motifs=cell(1,ndays);
roi_dates=cell(1,ndays);
roi_filenames=cell(1,ndays);
roi_times=cell(1,ndays);
roi_params=[];

for i=1:length(ROI_data_cleansed)

  ntrials=length(ROI_data_cleansed{i}.align_detrended);

  roi_data{i}=cat(3,ROI_data_cleansed{i}.align_detrended{:});
  roi_motifs{i}=zeros(1,ntrials);
  roi_dates{i}=zeros(1,ntrials);

  roi_filenames{i}=ROI_data_cleansed{i}.filename;
  roi_times{i}=ROI_data_cleansed{i}.frame_idx{1}./24.414e3;

  roi_params(i).fs=22.0144;
  roi_params(i).padding=[.25 .75];

  for j=1:ntrials
      tmp=ROI_data_cleansed{i}.filename{j};
      tokens=regexp(tmp,'\_(\d+)\_','tokens');
      roi_dates{i}(j)=datenum([ tokens{1}{1} tokens{2}{1} ],'yyyymmddHHMMSS');
      tokens2=regexp(tmp,'\_(\d+)$','tokens');
      roi_motifs{i}(j)=str2num(tokens2{1}{1});
  end

end

save(fullfile(dirs.agg_dir,dirs.ca_dir,['lw76-' ext '.mat']),'roi_data',...
  'roi_dates','roi_times','roi_motifs','roi_filenames','roi_params');
