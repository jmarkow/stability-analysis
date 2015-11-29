%
%
%

[options,dirs]=stan_preflight;

% load in each mat file, make spectral density image, etc.

% assume we have contours1

% get listing

load(fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis','sdi_analysis_data.mat'),'tfdensity');

% analyze circadian variation within-day/overnight shift


score_day2night={};
score_night2night={};
score_day2day={};

for i=1:length(tfdensity)

	ndays=length(tfdensity{i}.day);
	i

	for j=1:ndays

		template=mean(tfdensity{i}.night{j},3);
		norm1=sum(template(:).^2);
		template2=mean(tfdensity{i}.day{j},3);
		norm12=sum(template2(:).^2);

		for k=j:ndays

			lag=(tfdensity{i}.dates(k)-tfdensity{i}.dates(j))+1;

			% get p(score) from each segment

			ntrials=size(tfdensity{i}.day{k},3);

			for l=1:ntrials

				% replace with ssim over average images?

				cur_contour_day=double(tfdensity{i}.day{k}(:,:,l));
				norm2=sum(cur_contour_day(:).^2);

				cur_contour_night=double(tfdensity{i}.night{k}(:,:,l));
				norm3=sum(cur_contour_night(:).^2);

				score_day2night{i}{j,k}(l)=sum(sum(template.*cur_contour_day))/sqrt(norm1*norm2);
				score_day2day{i}{j,k}(l)=sum(sum(template2.*cur_contour_day))/sqrt(norm12*norm2);
				score_night2night{i}{j,k}(l)=sum(sum(template.*cur_contour_night))/sqrt(norm1*norm3);

			end
		end
	end
end
