function store_spikes=stan_ephys_xproj_scan()
% drift in any peaks > 100 Hz

[status,result]=unix(['find ' pwd ' -type f -name "spikedata_ch*.mat"']);
pkfile=regexp(result,'\n','split');
condition='nocar'; % use car or nocar fields

bound=[.1];
binsize=.001;
proc_fs=1e3;
sig=.005;
xval=[-3*sig:1/proc_fs:3*sig];
kernel=normpdf(xval,0,sig);
kernel=kernel./sum(kernel);
store_spikes=[];

for i=1:length(pkfile)-1

	load(pkfile{i},'cluster','proc_data');
	[path,file,ext]=fileparts(pkfile{i});

	tokens=regexp(path,filesep,'split');

	birdid=tokens{end-2};
	birdarea=tokens{end};

	fid=fopen(fullfile(path,'cellinfo.txt'),'r');

	readdata=textscan(fid,'%s%[^\n]','commentstyle','#',...
		'delimiter','\t','MultipleDelimsAsOne',1);

	% close the file

	fclose(fid);

	% read in cluster number from cellinfo.txt

	clusternum=str2num(readdata{2}{find(strcmpi(readdata{1},'cluster:'))});
	channel=str2num(readdata{2}{find(strcmpi(readdata{1},'channel:'))});

	% how many peaks in smooth rate? take multi-peak neurons

	[nsamples]=size(proc_data);
  ntrials=max(cluster.trials{clusternum});

	start_pt=ceil(bound*proc_fs);
	stop_pt=floor(((nsamples/cluster.parameters.fs)*proc_fs)-(bound*proc_fs));

	bin_edges=[start_pt:stop_pt];
	binspikes=zeros(ntrials,length(bin_edges)-1);
	smoothspikes=zeros(size(binspikes));

	% could have option to clean spikes first (for projection neurons especially)

	for j=1:ntrials

		tmp=cluster.times{clusternum}(cluster.trials{clusternum}==j);

		% convert spike time to samples in proc_fs

		tmp=round((tmp./cluster.parameters.fs)*proc_fs);
		bins=histc(tmp,bin_edges);

		binspikes(j,:)=binspikes(j,:)+bins(1:end-1);
		smoothspikes(j,:)=conv(binspikes(j,:),kernel,'same');

	end

	store_spikes(end+1).binned_spikes=binspikes';
	store_spikes(end).smoothed_spikes=smoothspikes';
	store_spikes(end).bin_edges=bin_edges(1:end-1)'/proc_fs;
	store_spikes(end).fs=proc_fs;
	store_spikes(end).windows=cluster.windows{clusternum};
	store_spikes(end).times=cluster.times{clusternum};
	store_spikes(end).trials=cluster.trials{clusternum};
	store_spikes(end).bird=birdid;

end

%randstore=[];
