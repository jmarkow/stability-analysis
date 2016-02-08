function KEY=stan_template_key()
%
%
%
%
[options,dirs]=stan_preflight;
temp_dir=fullfile(dirs.agg_dir,dirs.template_dir);

% spits out a key with filenames, birdnames and motif names for template files in the template directory

tmp=dir(fullfile(temp_dir,'*.mat'));
listing={tmp(:).name};


birdname=regexp(listing,'^([a-zA-Z0-9]+)_.*','tokens');
motif_name=regexp(listing,'^[a-zA-Z0-9]+_(.*)','tokens');

KEY=struct();

for i=1:length(listing)
	[bird_name,motif_name]=stan_read_filename(listing{i});
	KEY(i).filename=fullfile(temp_dir,listing{i});
	KEY(i).bird_name=bird_name;
	KEY(i).motif_name=motif_name;
end
