function stan_mat_batch()
% searches for mat files and associated log files, processes and copies to target directory
%

alias_name='aliases.txt';

log_file(1).field='logfile';
log_file(1).filename='log.txt';
log_file(1).multi=0;

% returns filenames to process and their associated log file

filenames=robofinch_dir_recurse(pwd,'*.mat',[],[],log_file);
[log_names,~,log_id]=unique({filenames(:).logfile});

% get aliases

cur_file=mfilename('fullpath')
[cur_path,~,~]=fileparts(cur_file)

[aliases.targets,aliases.sources,aliases.date_targets,aliases.date_sources]=stan_read_aliases(fullfile(cur_path,alias_name));
aliases
pause();
for i=1:length(log_names)	
	
	log_map=stan_read_config(log_names{i});

	pause(.05);

	proc_idx=find(log_id==i);
	files_to_proc=filenames(proc_idx);


	for j=1:length(log_map)

		% first check for datenumbers
	
		if ~isempty(log_map(j).date_num) & ~isempty(aliases.date_sources)
			
			alias_date_idx=(log_map(j).date_num==aliases.date_sources);

			if any(alias_date_idx)
				log_map(j).name=aliases.date_targets{alias_date_idx};
				continue;
			end
		end	
		
		alias_idx=strcmpi(log_map(j).name,aliases.sources);

		if any(alias_idx)
			log_map(j).name=aliases.targets{alias_idx};
		end

	end

	for j=1:length(log_map)
		disp([log_map(j).name]);
	end

	% scan for aliases, convert dates to datenums

	for j=1:length(files_to_proc)

		% process files with the same log id

	end

end

