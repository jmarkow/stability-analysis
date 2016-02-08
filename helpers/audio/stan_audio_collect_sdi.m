function stan_audio_collect_sdi()
%
%
%

[options,dirs]=stan_preflight;

% load in each mat file, make spectral density image, etc.

% assume we have contours1

% get listing

workdir=fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis');
files=robofinch_dir_recurse(workdir,'sdi.mat');

% check all files for bird id
birdid=cell(1,length(files));
dates=zeros(1,length(files));
filenames={files(:).name};

for i=1:length(filenames)
	splits=regexp(filenames{i},filesep,'split');
	token_string=splits{end-1};
	tmp=regexp(token_string,'(\w+)\_','tokens');
	birdid{i}=tmp{1}{1};
	dates(i)=datenum(regexp(token_string,'\d+\-\d+\-\d+','match'),'yyyy-mm-dd');
end

load(filenames{1},'sdi');

frac=2;
freq_use=[500 7e3];
t_use=[.2 .2]

tfdensity=[];
simscores=[];

birds=unique(birdid);
tfdensity=cell(1,length(birds));

for i=1:length(birds)

	idx=strcmp(birdid,birds{i});
	contour_files=filenames(idx);

	load(contour_files{1},'sdi');

	crop_f=[];
	crop_t=[];

	if freq_use(1)>0
		crop_f(1)=max(find(sdi.f<=freq_use(1)));
	else
		crop_f(1)=1;
	end

	crop_f(2)=min(find(sdi.f>=freq_use(2)));
	crop_f=crop_f(1):crop_f(2);

	crop_t(1)=max(find(sdi.t<=t_use(1)));
	crop_t(2)=min(find(sdi.t>=sdi.t(end)-t_use(2)));
	crop_t=crop_t(1):crop_t(2);

	% tfdensity{i}.day=zeros(length(crop_f),length(crop_t),length(contour_files));
	% tfdensity{i}.night=zeros(size(tfdensity{i}.day));
	% tfdensity{i}.all=zeros(size(tfdensity{i}.day));
	tfdensity{i}.day=cell(1,length(contour_files));
	tfdensity{i}.night=cell(size(tfdensity{i}.day));

	tfdensity{i}.filenames=contour_files;
	tfdensity{i}.dates=dates(idx);

	for j=1:length(contour_files)

		disp([contour_files{j}]);
		load(contour_files{j},'sdi');
		[ntrials]=size(sdi.contours.im,3);

		group1=1:floor(ntrials/frac);
		group2=ntrials-(floor(ntrials/frac)-1):ntrials;

		% tfdensity{i}.day(:,:,j)=mean(sdi.contours.im(crop_f,crop_t,group1),3);
		% tfdensity{i}.night(:,:,j)=mean(sdi.contours.im(crop_f,crop_t,group2),3);
		% tfdensity{i}.all(:,:,j)=mean(sdi.contours.im(crop_f,crop_t,:),3);

		tfdensity{i}.day{j}=sdi.contours.im(crop_f,crop_t,group1);
		tfdensity{i}.night{j}=sdi.contours.im(crop_f,crop_t,group2);

	end
end

save(fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis','sdi_analysis_data.mat'),'tfdensity','-v7.3');
