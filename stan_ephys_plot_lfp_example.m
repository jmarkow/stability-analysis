function fig=stan_ephys_plot_lfp(COLORS)
%
%
%
%
%
%

if nargin<1 | isempty(COLORS)
	COLORS='jet';
end


[options,dirs]=stan_preflight;

option_names=fieldnames(options);

% which options specify control raster

idx=regexp(option_names,'lfp_example(\d+)_(\w+)','tokens'); 

% gather all rasters specified in options

lfp=struct();

for i=1:length(idx)
	if ~isempty(idx{i})
		example_number=str2num(idx{i}{1}{1});
		lfp(example_number).(idx{i}{1}{2})=options.(option_names{i});
	end
end

template_key=stan_read_templates();
new_fs=1e3;

for i=1:length(lfp)

	load(fullfile(dirs.agg_dir,lfp(i).path),'store');

	% grab the template
	
	[pathname,filename,ext]=fileparts(lfp(i).path);
	tokens=regexp(filename,'\_','split');
	bird_id=tokens{1};

	load(fullfile(dirs.agg_dir,dirs.template_dir,[ bird_id '_' lfp(i).motif_name '.mat']),'template')

	lfp(i)

	motif_list={store(:).motif_name};
	motif_idx=strcmp(motif_list,lfp(i).motif_name);

	ch_idx=find(store(motif_idx).ch_list==lfp(i).channel);

	lfp_data=store(motif_idx).lfp.data{ch_idx};
	
	[b,a]=ellip(4,.2,40,[25 35]/(options.lfp_fs/2),'bandpass');

	filt_lfp=cellfun(@(x) filtfilt(b,a,x),lfp_data,'uniformoutput',0);
	mu_lfp=cellfun(@(x) mean(zscore(x),2),filt_lfp,'uniformoutput',0);

	% interpolate?

	mu_lfp=zscore(cat(2,mu_lfp{:}));
	
	spacing=6;
	spacing=repmat(spacing,[1 length(lfp(i).days)-1]);
	spacing=[0 cumsum(spacing)];
	
	[nsamples,ntrials]=size(mu_lfp);

	%upsample_factor=new_fs/options.lfp_fs;
	%mu_t2=[1:nsamples*upsample_factor]/new_fs;

	spacing=repmat(spacing,[nsamples 1]);
	mu_t=[1:nsamples]/options.lfp_fs;

	fig.(bird_id)=figure();

	[s,f,t]=zftftb_pretty_sonogram(template.data,template.fs,'filtering',300,'len',70,'overlap',69.5,'clipping',[-3 2],'zeropad',0);
	ax(1)=subplot(3,1,1);
	imagesc(t+options.padding_lfp(1),f/1e3,s);
	axis xy;
	ylim([0 10]);
	colormap(COLORS);
	set(gca,'XTick',[],'YTick',[]);

	ax(2)=subplot(3,1,2:3);
	plot(mu_t,mu_lfp(:,lfp(i).days)+spacing,'b-');

	linkaxes(ax,'x');

	if isfield(lfp(i),'xlim')
		xlim([lfp(i).xlim]);
	end


end
