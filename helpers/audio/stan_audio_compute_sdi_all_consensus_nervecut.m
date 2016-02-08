function stan_audio_compute_sdi_all_consensus_nervecut()

[~,dirs]=stan_preflight();
listing=robofinch_dir_recurse(fullfile(dirs.agg_dir,dirs.nervecut_audio_dir),'*.mat');
storedir='sdi';

savedir=fullfile(dirs.agg_dir,dirs.nervecut_audio_dir,storedir);

if ~exist(savedir,'dir')
  mkdir(savedir);
end

for i=1:length(listing)

  disp([listing(i).name]);
  [pathname,filename,ext]=fileparts(listing(i).name);

  if exist(fullfile(savedir,[filename '_sdi.mat']),'file')
    continue;
  end

  load(listing(i).name,'audio');
  [consensus,f,t]=acontour(audio.data(:,1),audio.fs);

  % compute sdi, relevant quantities, store

  ntrials=size(audio.data,2);
  [rows,columns]=size(consensus);
  consensus=zeros(rows,columns,ntrials,'single');

  % parfor is choking on the structure, remapping to variables

  use_data=audio.data;
  fs=audio.fs;
  sdi.fs=audio.fs;

  clearvars audio;

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
  save(fullfile(savedir,[filename '_sdi.mat']),'sdi','-v7.3');
  clearvars sdi;

end
