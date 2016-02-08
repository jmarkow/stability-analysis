function SU_DATA=stan_ephys_getdata_su_fr
%
%
%
%

[options,dirs]=stan_preflight;

if nargin<1 | isempty(COLORS)
	COLORS='jet';
end

kernedges=[-3*options.smooth_sig:1/options.smooth_fs:3*options.smooth_sig];
kernel=(1/(options.smooth_sig*sqrt(2*pi)))*exp((-(kernedges-0).^2)./(2*options.smooth_sig^2));

% get listing of single units
%

tmp=dir(fullfile(dirs.agg_dir,dirs.su_analyze_dir));
unit_dirs={};
unit_names={};

for	i=1:length(tmp)
	if tmp(i).isdir & tmp(i).name(1)~='.'
		unit_dirs{end+1}=fullfile(dirs.agg_dir,dirs.su_analyze_dir,tmp(i).name);
		unit_names{end+1}=tmp(i).name;
	end
end


SU_DATA=struct();

for i=1:length(unit_dirs)

	listing=dir(fullfile(unit_dirs{i},'*.mat'));
	nplots=length(listing);

	plot_spikes=struct();

	tmp=regexp(listing(1).name,'^(\w+)\_','tokens');
	bird_id=tmp{1}{1};
	date_num=[];	

	SU_DATA(i).smooth_fs=options.smooth_fs;
	SU_DATA(i).smooth_sig=options.smooth_sig;

	for j=1:length(listing)
		
		% get cluster number

		load(fullfile(unit_dirs{i},listing(j).name),'cluster','su_info','proc_data');

		tmp=regexp(listing(j).name,'\_(\d+-\d+-\d+)\_','tokens');
		date_num=datenum(tmp{1}{1});

		[nsamples,ntrials]=size(proc_data);
		clust_num=su_info.cluster;

		spikes=cluster.times{clust_num};
		spikes=spikes/cluster.parameters.fs;
		trials=cluster.trials{clust_num};

		uniq_trials=unique(trials);

		smooth_samples=ceil((nsamples/cluster.parameters.fs)*options.smooth_fs);
		smooth_rate=zeros(smooth_samples,ntrials);

		for k=1:length(uniq_trials)
			spikes_smps=round(spikes(trials==uniq_trials(k))*options.smooth_fs);
			spikes_smps(spikes_smps==0)=1;
			smooth_rate(spikes_smps,uniq_trials(k))=1;
			smooth_rate(:,uniq_trials(k))=conv(smooth_rate(:,k),kernel,'same');
		end

		SU_DATA(i).datenum(j)=date_num;
		SU_DATA(i).smooth_rate{j}=smooth_rate;
		
	end

end

