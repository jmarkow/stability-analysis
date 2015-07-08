function stan_mat_batch()
% searches for mat files and associated log files, processes and copies to target directory
%

alias_name='aliases.txt';
options_name='options.txt';
datefmt='yymmdd_HHMMSS';
recid='barecarbon';
gain_factor=1e3;
newfs=20e3;
anti_alias=9e3;
delim='_';
bird_delim='&';
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
options=stan_read_options(fullfile(cur_path,options_name));

for i=1:length(log_names)	

	[log_map,map]=stan_read_config(log_names{i});

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

	% construct new filename using new format

	base_filename='';	

	for j=1:length(log_map)

		if isempty(log_map(j).ch)
			continue;
		end

		log_map(j).ch.idx=log_map(j).ch.idx-1; % subtract 1 to match NIDAQ indexing

		if j>1
			base_filename=[ base_filename bird_delim ];
		end

		base_filename=[ base_filename log_map(j).name delim recid ];

		if any(log_map(j).ch.ismic)
			base_filename=[ base_filename sprintf('%smic%iadc',delim,log_map(j).ch.idx(log_map(j).ch.ismic)) ];
		end

		if any(~log_map(j).ch.ismic)

			idxs=log_map(j).ch.idx(find(~log_map(j).ch.ismic));
			idxs=sort(idxs);

			if length(idxs)>1
				base_filename=[ base_filename sprintf('%sdata%i-%iadc',delim,idxs(1),idxs(end)) ];
			else
				base_filename=[ base_filename sprintf('%sdata%iadc',delim,idxs) ];
			end
		end

	end

	% scan for aliases, convert dates to datenums

	for j=1:length(files_to_proc)

		% process files with the same log id

		load(files_to_proc(j).name,'data');

		[nsamples,nchannels]=size(data.voltage);

		data2.voltage=data.voltage/gain_factor;	% correct amplifier gain, resample at sensible frequency (~20 kHz)
	
		data2.time=data.time-min(data.time);
		data2.start_time=data.start_time;
		data2.start_time(end)=data2.start_time(end)+data.time(1);
		data2.fs=data.sampling_rate;
		data2.labels=[0:nchannels-1];
		data2.names=map.names;

		decimate_f=round(data2.fs/newfs);
		%cutoff=.8*(data2.fs/2)/decimate_f
		cutoff=anti_alias/(data2.fs/2);
		[b,a]=ellip(4,.2,40,cutoff,'low');

		data2.voltage=downsample(filtfilt(b,a,data2.voltage),decimate_f);
		data2.time=downsample(data2.time,decimate_f);
		data2.fs=newfs;

		new_filename=[ base_filename delim datestr(datenum(data2.start_time),datefmt) '.mat' ];

		% get names for each channel
		
		data=data2;
		clear data2;
		save(fullfile(options.convert_destination,new_filename),'data');
		
	end
end

