%
%
%

[options,dirs]=stan_preflight;

% load in each mat file, make spectral density image, etc.

cur_dir=fullfile(dirs.agg_dir,dirs.sdi_dir,'use_data');
listing=dir(fullfile(cur_dir,'*.mat'));

for i=1:length(listing)

	disp([listing(i).name]);
	load(fullfile(cur_dir,listing(i).name),'agg_audio');

	% compute sdi, relevant quantities, store

	size(agg_audio.data)
	[data(i).sdi data(i).f data(i).t data(i).contours]=zftftb_sdi(agg_audio.data,agg_audio.fs,...
		'spect_thresh',.8,'len',70,'overlap',69,'tscale',1.5);
	clear agg_audio;

end

save(fullfile(cur_dir,'..','sdi_data.mat'),...
	'data','-v7.3');
