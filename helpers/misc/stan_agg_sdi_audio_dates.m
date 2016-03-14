function stan_agg_nervecut_audio()
% copies audio data from nervecut birds, finds last suitable pre and earliest suitable post
%

% get options

[options,dirs]=stan_preflight;

% gather dates, copy only mic data (stitch too?)
%

listing=robofinch_dir_recurse(pwd,'aggregated_data.mat');

if isempty(listing)
	return;
end

tmp=regexp(listing(1).name,'((\d+-)+\d+)','match');
date_number=datenum(tmp);
tmp=regexp(listing(1).name,filesep,'split');
birdid=tmp{end-7};

storedir=fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis',[ birdid '_' datestr(date_number,'yyyy-mm-dd') ]);

if ~exist(storedir,'dir')
	mkdir(storedir);
end

timestamps=[];

for i=1:length(listing)
	disp([listing(i).name]);
	vars=whos('-file',listing(i).name);
	varnames={vars(:).name};
	if ~strcmp(varnames,'MIC_DATA')
  	tmp=load(listing(i).name,'agg_file_datenum');
		timestamps=[timestamps agg_file_datenum];
	else
		disp('Legacy...')
		load(listing(i).name,'START_DATENUM')
		timestamps=[timestamps START_DATENUM];
	end
end

save(fullfile(storedir,['mic_data_timestamps.mat']),'timestamps');
