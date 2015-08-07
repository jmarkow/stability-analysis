function stan_process_birds(DIR,varargin)
%
%
%
%

proc_name='robomulab.mat';
stan_dir='stan';
save_filename='stan_withinbird';
newfs=500;

% load in MU calculation, compare across days

% this script should be run in the base ephys directory (e.g. ephys_barecarbon)
% we're assuming next dir is bird, then rec id, etc.
%

% base_dir
% -->bird
% ----->recid
% ------->date


% recurse and find mu data

if nargin<1 | isempty(DIR), DIR=pwd; end


tmp=dir(DIR);
bird_listing={};

for i=1:length(tmp)
	if tmp(i).isdir & tmp(i).name(1)~='.'
		bird_listing{end+1}=tmp(i).name;
	end
end

%  need a key for motif name and channel

for i=1:length(bird_listing)

	disp([bird_listing{i}]);
	
	% get rec_id
	
	tmp=dir(fullfile(DIR,bird_listing{i}));
	rec_listing={};

	for j=1:length(tmp)
		if tmp(j).isdir & tmp(j).name(1)~='.'
			rec_listing{end+1}=tmp(j).name;
		end
	end


	for j=1:length(rec_listing)

		disp([rec_listing{j}]);	

		cur_bird=bird_listing{i};
		cur_recid=rec_listing{j};
		cur_dir=fullfile(DIR,cur_bird,cur_recid);

		storage_dir=fullfile(cur_dir,stan_dir);
		cur_savefile=[ cur_bird '_' cur_recid '_' save_filename '.mat' ];
		
		%if exist(fullfile(storage_dir,cur_savefile),'file'), continue; end

		proc_files=robofinch_dir_recurse(cur_dir,proc_name);
	
		store=struct();
		motif_list={};
		counter=1;

		for k=1:length(proc_files)

			disp([proc_files(k).name]);

			% extract date/template name from the directory name

			date_token=regexp(proc_files(k).name,[ filesep '(\d+-\d+-\d+)' filesep ],'tokens');
			date_number=datenum(date_token{1}{1});

			motif_token=regexp(proc_files(k).name,[ filesep '(\w+)_roboextract' filesep ],'tokens');
			motif_name=motif_token{1}{1};

			% for each channel, will have motif, channel, date triplet
			%
			% initialize structure if we have a new motif
		
			motif_idx=strcmp(motif_list,motif_name);
			new_motif=isempty(motif_idx)|~any(motif_idx);

			if new_motif
				motif_list{end+1}=motif_name;
				motif_idx=length(motif_list);
				store(motif_idx).bird_id=cur_bird;
				store(motif_idx).rec_id=cur_recid;	
				store(motif_idx).motif_name=motif_name;
				store(motif_idx).ch_list={};
				store(motif_idx).datenums=[];
				store(motif_idx).rms.data={};
				store(motif_idx).rms.parameters={};
				store(motif_idx).spikes.smooth_rate={};	
				store(motif_idx).spikes.times={};
				store(motif_idx).spikes.trial={};
				store(motif_idx).spikes.threshold={};
				store(motif_idx).spikes.parameters={};

			end

			% map to appropriate motif element, now grab channel list from this motif

			ch_list=store(motif_idx).ch_list;
			
			load(proc_files(k).name,'rms','spikes','parameters');

			tmp=parameters.ch_names;
			tmp=strtrim(regexprep(tmp,'\[.*\]',''));
			ch_names=tmp;

			% map to channel

			for l=1:length(ch_names)
				
				ch_idx=strcmp(ch_list,ch_names{l});
				new_ch=isempty(ch_idx)|~any(ch_idx);
				
				if new_ch
					store(motif_idx).ch_list{end+1}=ch_names{l};
					ch_idx=length(store(motif_idx).ch_list);
					col_idx=1;
				else
					col_idx=length(store(motif_idx).spikes.smooth_rate{ch_idx})+1;
				end

				% downsample, assume we've already done the necessary smoothing

				down_fact_spikes=round(spikes{l}.fs/newfs);
				down_fact_rms=round(rms{l}.fs/newfs);

				store(motif_idx).datenums(ch_idx,col_idx)=date_number;	
				store(motif_idx).parameters(ch_idx,col_idx)=parameters;
				store(motif_idx).spikes.smooth_rate{ch_idx}{col_idx}=downsample(spikes{l}.smooth_rate,down_fact_spikes);
				store(motif_idx).spikes.times{ch_idx}{col_idx}=spikes{l}.times;
				store(motif_idx).spikes.trial{ch_idx}{col_idx}=spikes{l}.trial;
				store(motif_idx).spikes.threshold{ch_idx}{col_idx}=spikes{l}.threshold;
				store(motif_idx).spikes.parameters{ch_idx}{col_idx}=spikes{l}.smooth_params;
				store(motif_idx).spikes.parameters{ch_idx}{col_idx}.smooth_fs=newfs;
				store(motif_idx).spikes.parameters{ch_idx}{col_idx}.spike_fs=spikes{l}.fs;
				store(motif_idx).rms.data{ch_idx}{col_idx}=downsample(rms{l}.data,down_fact_rms);
				store(motif_idx).rms.parameters{ch_idx}{col_idx}=rms{l}.params;
				store(motif_idx).rms.parameters{ch_idx}{col_idx}.newfs=newfs;

			end
		end

		if length(proc_files)>0
			if ~exist(storage_dir,'dir'), mkdir(storage_dir), end;
			save(fullfile(storage_dir,cur_savefile),'store','-v7.3');
		end

	end
end

