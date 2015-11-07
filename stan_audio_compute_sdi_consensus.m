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

	for j=1:size(agg_audio.data,2)
		disp([num2str(j) ' of ' num2str(size(agg_audio.data,2))]);
		[data(i).consensus{j} data(i).f{j} data(i).t{j} data(i).auditory_contours{j} data(i).sonogram{j}]=...
			acontour(agg_audio.data(:,j),agg_audio.fs);
	end

	clear agg_audio;

end

save(fullfile(cur_dir,'..','sdi_data_consensus.mat'),...
	'data','-v7.3');
