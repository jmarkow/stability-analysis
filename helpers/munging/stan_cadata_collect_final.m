% simply collects the calcium data from Will's folder, formats and saves to
% a more convenient format for further analysis...

[opts,dirs]=stan_preflight;
ext='lib';

listing=robofinch_dir_recurse(pwd,['ROI_data_cleansed-' ext '.mat']);

% all except lw76


for i=1:length(listing)-1
  load(listing(i).name);

  tokens=regexp(listing(i).name,filesep,'split');
  [datapath,datafile,~]=fileparts(listing(i).name)
  new_listing=robofinch_dir_recurse(datapath,'*-bad.gif');
  roi_exclude=cell(1,length(ROI_data_cleansed));

  for j=1:length(ROI_data_cleansed)

      to_del=[];

      for k=1:length(ROI_data_cleansed{j}.filename(2,:))

          for l=1:length(new_listing)

            [~,match_file,~]=fileparts(new_listing(l).name);
            match_file=regexprep(match_file,'-bad','');

            if strcmp(ROI_data_cleansed{j}.filename{2,k},match_file)

                fprintf('Removing trial related to file %s\n',match_file);
                fprintf('Day %g trial %g\n',j,k);
                to_del(:,end+1)=k;

            end
          end
      end

      ROI_data_cleansed{j}.raw_time(:,to_del)=[];
      ROI_data_cleansed{j}.raw_dat(:,to_del)=[];
      ROI_data_cleansed{j}.filename(:,to_del)=[];
      ROI_data_cleansed{j}.interp_dff(:,to_del)=[];
      ROI_data_cleansed{j}.interp_raw(:,to_del)=[];
      ROI_data_cleansed{j}.interp_time(to_del)=[];
      ROI_data_cleansed{j}.analogIO_dat(to_del)=[];
      ROI_data_cleansed{j}.analogIO_time(to_del)=[];

  end


  [roi_data,roi_dates,roi_times,roi_motifs,roi_trials,roi_filenames,roi_params]=...
              stan_cadata_collect_freedomscope_v2_cell(ROI_data_cleansed);
  roi_audio=cell(1,length(roi_data));

  for j=1:length(roi_data)

    % get mic data

    roi_audio{j}=cat(2,ROI_data_cleansed{j}.analogIO_dat{:});
    roi_audio{j}=roi_audio{j}(:,roi_trials{j});
    roi_params(j).audio_fs=48e3;

  end

  save(fullfile(datapath,[datafile '-badremoved.mat']),'ROI_data_cleansed','-v7.3');
  save(fullfile(dirs.agg_dir,dirs.ca_dir,[lower(tokens{end-1}) '-' ext '.mat']),'roi_data',...
    'roi_dates','roi_times','roi_motifs','roi_filenames','roi_params','roi_audio','roi_trials');

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
roi_audio=cell(1,ndays);
roi_trials=cell(1,ndays);
roi_params=[];

[datapath,datafile,~]=fileparts(listing(end).name);
new_listing=robofinch_dir_recurse(datapath,'*-bad.gif');
roi_exclude=cell(1,length(ROI_data_cleansed));

 for j=1:length(ROI_data_cleansed)

     to_del=[];

     for k=1:length(ROI_data_cleansed{j}.filename)
         for l=1:length(new_listing)
             [~,match_file,~]=fileparts(new_listing(l).name);
             match_file=regexprep(match_file,'-bad','');
             if strcmp(ROI_data_cleansed{j}.filename{k},match_file)

                fprintf('Removing trial related to file %s\n',match_file);
                fprintf('Day %g trial %g\n',j,k);
                to_del(:,end+1)=k;

             end
         end
     end

    ROI_data_cleansed{j}.align_detrended(to_del)=[];
    ROI_data_cleansed{j}.frame_idx(to_del)=[];
    ROI_data_cleansed{j}.mic_data(to_del)=[];
    ROI_data_cleansed{j}.filename(to_del)=[];

 end

for i=1:length(ROI_data_cleansed)

  ntrials=length(ROI_data_cleansed{i}.align_detrended);

  roi_data{i}=cat(3,ROI_data_cleansed{i}.align_detrended{:});
  roi_motifs{i}=zeros(1,ntrials);
  roi_trials{i}=1:ntrials;
  roi_dates{i}=zeros(1,ntrials);
  roi_audio{i}=cat(1,ROI_data_cleansed{i}.mic_data{:})';

  roi_filenames{i}=ROI_data_cleansed{i}.filename;
  roi_times{i}=ROI_data_cleansed{i}.frame_idx{1}./24.414e3;

  roi_params(i).fs=ROI_data_cleansed{i}.movie_fs;
  roi_params(i).audio_fs=24.414e3;
  roi_params(i).padding=[.25 .75];

  for j=1:ntrials
      tmp=ROI_data_cleansed{i}.filename{j};
      tokens=regexp(tmp,'\_(\d+)\_','tokens');
      roi_dates{i}(j)=datenum([ tokens{1}{1} tokens{2}{1} ],'yyyymmddHHMMSS');
      tokens2=regexp(tmp,'\_(\d+)$','tokens');
      roi_motifs{i}(j)=str2num(tokens2{1}{1});
  end

end

save(fullfile(datapath,[datafile '-badremoved.mat']),'ROI_data_cleansed','-v7.3');
save(fullfile(dirs.agg_dir,dirs.ca_dir,['lw76-' ext '.mat']),'roi_data',...
  'roi_dates','roi_times','roi_motifs','roi_filenames','roi_params','roi_audio','roi_trials');
