function stan_agg_nervecut_audio()
% copies audio data from nervecut birds, finds last suitable pre and earliest suitable post
%

dirs_name='dirs.txt';
options_name='options.txt';
dry_run=0;
save_dir='mu_agg';
agg_file='roboaggregate.mat';
max_depth=4;

% get options

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);
dirs=stan_read_options(fullfile(cur_path,dirs_name));
options=stan_read_options(fullfile(cur_path,options_name));

% all proc data goes into the same directory

save_dir=fullfile(dirs.agg_dir,dirs.nervecut_audio_dir);
if ~exist(save_dir,'dir'), mkdir(save_dir); end

tmp=dir(dirs.data_dir);
birds={};

for i=1:length(tmp)
	if tmp(i).isdir & tmp(i).name(1)~='.'
		birds{end+1}=fullfile(dirs.data_dir,tmp(i).name);
	end
end


for i=1:length(birds)

	% get directories in this directory
	
	listing=dir(fullfile(birds{i},'barecarbon_nervecut*'));

	if isempty(listing)
		continue;
	end

	dirs={listing(:).name};

	% first directory is pre, second is post
	
	listing=dir(fullfile(birds{i},dirs{1}));
	dates={};

	for j=1:length(listing)
		if listing(j).isdir & listing(j).name~='.' ...
				& ~strcmp(listing(j).name,'stan') ...
				& ~strcmp(listing(j).name,'templates')
			dates{end+1}=fullfile(birds{i},dirs{1},listing(j).name);
		end
	end

	% travel backwards through listing until we get a hit

	motif_list={};

	for j=length(dates):-1:1
		agg_files=robofinch_dir_recurse(dates{j},agg_file);

		for	k=1:length(agg_files)

			tmp=regexp(agg_files(k).name,'\/(\w+)\_roboextract\/','tokens')

			if ~isempty(tmp)
				motif_name=tmp{1}{1};
				
				if ~isempty(motif_list) & any(strcmp(motif_list,motif_name))
					continue;
				end

				load(agg_files(k).name,'audio','file_datenum');

				[~,ntrials]=size(audio.data);
				trial_flag=ntrials>=options.audio_trial_limit;

				% filename is bird name, condition, motif
				tokens=regexp(birds{i},filesep,'split');
				
				savefile= fullfile(save_dir,[ tokens{end} '_' motif_name '.mat' ]);

				if trial_flag
					motif_list{end+1}=motif_name;
					disp(['Saving ' savefile '...']);
					save(savefile,'audio','file_datenum');
				end

				
				% save data

			end
		end
	end


	% now step forward for post data
	
	listing=dir(fullfile(birds{i},dirs{2}));
	dates={};

	for j=1:length(listing)
		if listing(j).isdir & listing(j).name~='.' ...
				& ~strcmp(listing(j).name,'stan') ...
				& ~strcmp(listing(j).name,'templates')
			dates{end+1}=fullfile(birds{i},dirs{2},listing(j).name);
		end
	end


	for j=1:length(dates)

		agg_files=robofinch_dir_recurse(dates{j},agg_file);

		for	k=1:length(agg_files)

			tmp=regexp(agg_files(k).name,'\/(\w+postcut\w+)\_roboextract\/','tokens')

			if ~isempty(tmp)
				motif_name=tmp{1}{1};
				
				if ~isempty(motif_list) & any(strcmp(motif_list,motif_name))
					continue;
				end

				load(agg_files(k).name,'audio','file_datenum');

				[~,ntrials]=size(audio.data);
				trial_flag=ntrials>=options.audio_trial_limit;

				tokens=regexp(birds{i},filesep,'split');

				savefile= fullfile(save_dir,[ tokens{end} '_' motif_name '.mat' ]);

				if trial_flag
					motif_list{end+1}=motif_name;
					disp(['Saving ' savefile '...']);
					save(savefile,'audio','file_datenum');
				end

				% save data

			end
		end
	end


end

