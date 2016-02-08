function stan_audio_compute_sdi_all_consensus()

[~,dirs]=stan_preflight();
listing=robofinch_dir_recurse(fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis'),'mic_data.mat');

for i=1:length(listing)

  disp([listing(i).name]);
  [pathname,filename,ext]=fileparts(listing(i).name);

  if exist(fullfile(pathname,'sdi_consensus.mat'),'file')
    continue;
  end

  load(listing(i).name,'agg_audio');

  [consensus,f,t]=acontour(agg_audio.data(:,1),agg_audio.fs);

  % compute sdi, relevant quantities, store

  ntrials=size(agg_audio.data,2);
  [rows,columns]=size(consensus);
  consensus=zeros(rows,columns,ntrials,'single');

  % parfor is choking on the structure, remapping to variables

  use_data=agg_audio.data;
  fs=agg_audio.fs;

  sdi.date=agg_audio.date;
  sdi.fs=agg_audio.fs;

  clearvars agg_audio;

  parfor j=1:ntrials
    disp([num2str(j) ' of ' num2str(ntrials)]);
    tmp=use_data(:,j);
    tmp2=single(acontour(tmp,fs));
    consensus(:,:,j)=tmp2;
  end

  sdi.consensus=consensus;
  sdi.f=f;
  sdi.t=t;

  clearvars consensus f t use_data fs;
  save(fullfile(pathname,'sdi_consensus.mat'),'sdi','-v7.3');

  clearvars sdi;
end
