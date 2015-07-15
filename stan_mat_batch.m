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
dry_run=0;

par_save = @(FILE,data) save(FILE,'data');
% returns filenames to process and their associated log file

filenames=robofinch_dir_recurse(pwd,'data_*.mat',[],[],log_file);
[log_names,~,log_id]=unique({filenames(:).logfile});

% get aliases

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);

[aliases.targets,aliases.sources,aliases.date_targets,aliases.date_sources]=stan_read_aliases(fullfile(cur_path,alias_name));
options=stan_read_options(fullfile(cur_path,options_name));

for i=1:length(log_names)	

	disp([log_names{i}]);

	[log_map,map]=stan_read_config(log_names{i});

	pause(.05);

	proc_idx=find(log_id==i);
	files_to_proc=filenames(proc_idx);

	[log_path,~,~]=fileparts(log_names{i});

	if exist(fullfile(log_path,'.convert_complete'),'file');
		continue;
	end

	for j=1:length(log_map)

		% first check for datenumbers

		if ~isempty(log_map(j).date_num) & ~isempty(aliases.date_sources)

			alias_date_idx=(log_map(j).date_num==aliases.date_sources);

			if any(alias_date_idx) & ~strcmpi(log_map(j).name,'lhp54')
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

	if dry_run
		continue;
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

	parfor j=1:length(files_to_proc)

		tmp=[];
		data=[];
		data2=[];

		fprintf('Processing %s\n',files_to_proc(j).name);

		% process files with the same log id

		try
			tmp=load(files_to_proc(j).name,'data');
		catch err
			warning('Could not read file %s',files_to_proc(j).name);
			continue;
		end

		if ~isfield(tmp,'data')
			warning('Could not read data from file %s',files_to_proc(j).name);
			continue;
		end

		data=tmp.data;
		tmp=[];

		[nsamples,nchannels]=size(data.voltage);

		data2.voltage=data.voltage;	% resample at sensible frequency 
	
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
		
		% store all relevant info 

		data2.parameters.units=repmat({'Volts'},[1 nchannels]);
		data2.parameters.sensor_range=[-10 10];
		data2.parameters.input_range=[-10 10];
		data2.parameters.units_range=[-10 10];
		data2.parameters.amp_gain=gain_factor;
		data2.parameters.gain_correct=false; % we did not yet adjust the ephys data by amp gain

		new_filename=[ base_filename delim datestr(datenum(data2.start_time),datefmt) '.mat' ];

		% get names for each channel
		
		data=data2;
		data2=[];

		stan_par_save(fullfile(options.convert_destination,new_filename),data);

	end

	fid=fopen(fullfile(log_path,'.convert_complete'),'w');
	fclose(fid);

end

