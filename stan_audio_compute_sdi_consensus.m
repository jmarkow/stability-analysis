%
%
%

[options,dirs]=stan_preflight;

% load in each mat file, make spectral density image, etc.

cur_dir=fullfile(dirs.agg_dir,dirs.sdi_dir,'use_data');
listing=dir(fullfile(cur_dir,'*.mat'));

% get dims

for i=1:length(listing)

	disp([listing(i).name]);
	load(fullfile(cur_dir,listing(i).name),'agg_audio');

	[consensus f t]=acontour(agg_audio.data(:,1),agg_audio.fs);
	% compute sdi, relevant quantities, store

	ntrials=size(agg_audio.data,2);
	[rows,columns]=size(consensus);

	consensus=zeros(rows,columns,ntrials,'single');

	parfor j=1:ntrials
		disp([num2str(j) ' of ' num2str(ntrials)]);
		[consensus(:,:,j)]=acontour(agg_audio.data(:,j),agg_audio.fs);
	end

	data(i).consensus=consensus;
	data(i).f=f;
	data(i).t=t;

	clearvars consensus f t agg_audio;

end

save(fullfile(cur_dir,'..','sdi_data_consensus.mat'),...
	'data','-v7.3');
