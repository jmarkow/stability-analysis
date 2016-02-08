function stan_audio_compute_sdi_all();

[options,dirs]=stan_preflight();

listing=robofinch_dir_recurse(fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis'),'mic_data.mat');

for i=1:length(listing)
  disp([listing(i).name]);
  [pathname,filename,ext]=fileparts(listing(i).name);
  %if exist(fullfile(pathname,'sdi.mat'),'file')
  %  continue;
  %end
  load(listing(i).name,'agg_audio');
  [sdi.histogram sdi.f sdi.t sdi.contours]=zftftb_sdi(agg_audio.data,agg_audio.fs,'spect_thresh',.8);
  sdi.date=agg_audio.date;
  sdi.fs=agg_audio.fs;
  save(fullfile(pathname,'sdi.mat'),'sdi','-v7.3');
end
