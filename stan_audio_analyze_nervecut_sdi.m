function simscores=stan_analyze_nervecut_audio()
%
%
%

% get SAP features, compare pre/post for best match motif
%


save_dir='features';

% get options

[options,dirs]=stan_preflight;

key=stan_read_nervecut_audio;

sdi_dir=fullfile(dirs.agg_dir,dirs.nervecut_audio_dir,'sdi');
listing=dir(fullfile(sdi_dir,'*.mat'));

[birds,~,idx]=unique({key(:).bird_name});
simscores=cell(length(birds),2,2);
spect_thresh=5;
song_thresh=0;

for i=1:2:length(listing)

	listing(i).name
	listing(i+1).name

	precut=load(fullfile(sdi_dir,listing(i).name),'sdi');
	postcut=load(fullfile(sdi_dir,listing(i+1).name),'sdi');

	sdi{1}=precut.sdi;
	sdi{2}=postcut.sdi;

	clearvars precut postcut;

	% get simscores

	ntrials=size(sdi{1}.consensus,3);

	template1=mean(sdi{1}.consensus>spect_thresh,3);
	template2=mean(sdi{2}.consensus>spect_thresh,3);

	template1_env=mean(template1);
	template2_env=mean(template2);

	% align templates using envelope

	[r,lags]=xcorr(template1_env,template2_env);

	[~,loc]=max(abs(r));
	shift=lags(loc);

	if shift<0
		template2=template2(:,-shift:end);
		sdi{2}.consensus=sdi{2}.consensus(:,-shift:end,:);
	elseif shift>0
		template1=template1(:,shift:end);
		sdi{1}.consensus=sdi{1}.consensus(:,shift:end,:);
	end

	% trim to proper size

	npoints=min(size(template1,2),size(template2,2));

	template1=template1(:,1:npoints);
	template2=template2(:,1:npoints);

	sdi{1}.consensus=sdi{1}.consensus(:,1:npoints,:);
	sdi{2}.consensus=sdi{2}.consensus(:,1:npoints,:);

	% edit out silence

	env=zscore(sum(template1));
	idx=env>song_thresh;

	template1=template1(:,idx);
	template2=template2(:,idx);

	sdi{1}.consensus=sdi{1}.consensus(:,idx,:);
	sdi{2}.consensus=sdi{2}.consensus(:,idx,:);

	template1_norm=sum(template1(:).^2);
	template2_norm=sum(template2(:).^2);

	contours=sdi{1}.consensus>spect_thresh;

	for j=1:ntrials
		cur_contour=contours(:,:,j);
		cur_norm=sum(cur_contour(:).^2);
		simscores{ceil(i/2),1,1}(j)=sum(sum(template1.*cur_contour))./sqrt(cur_norm*template1_norm);
		simscores{ceil(i/2),1,2}(j)=sum(sum(template2.*cur_contour))./sqrt(cur_norm*template2_norm);
	end

	ntrials=size(sdi{2}.consensus,3);

	contours=sdi{2}.consensus>spect_thresh;
	for j=1:ntrials
		cur_contour=contours(:,:,j);
		cur_norm=sum(cur_contour(:).^2);
		simscores{ceil(i/2),2,1}(j)=sum(sum(template1.*cur_contour))./sqrt(cur_norm*template1_norm);
		simscores{ceil(i/2),2,2}(j)=sum(sum(template2.*cur_contour))./sqrt(cur_norm*template2_norm);
	end
end
