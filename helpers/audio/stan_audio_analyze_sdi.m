%
%
%

[options,dirs]=stan_preflight;
filewrite=false;

% load in each mat file, make spectral density image, etc.

% assume we have contours1

% get listing

load(fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis','sdi_analysis_data_consensus.mat'),'tfdensity');

% analyze circadian variation within-day/overnight shift
thresh=-.5;

score_day2night={};
score_night2night={};
score_day2day={};
score_all2all={};
frac=2;

for i=1:length(tfdensity)

	ndays=length(tfdensity{i}.all);
	i

	for j=1:ndays

		ntrials=size(tfdensity{i}.all{j},3);

		group1=1:floor(ntrials/frac);
		group2=ntrials-(floor(ntrials/frac)-1):ntrials;

		template1=mean(tfdensity{i}.all{j}(:,:,group1),3);
		norm11=sum(template1(:).^2);

		template2=mean(tfdensity{i}.all{j}(:,:,group2),3);
		norm12=sum(template2(:).^2);

		template3=mean(tfdensity{i}.all{j},3);
		norm13=sum(template3(:).^2);

		idx=zscore(sum(template3))>thresh;

		template1=template1(:,idx);
		template2=template2(:,idx);
		template3=template3(:,idx);

		for k=j:ndays

			lag=(tfdensity{i}.dates(k)-tfdensity{i}.dates(j))+1;

			% get p(score) from each segment

			ntrials=size(tfdensity{i}.all{k},3);

			group1=1:floor(ntrials/frac);
			group2=ntrials-(floor(ntrials/frac)-1):ntrials;

			group_trials=length(group1);

			% get average time difference?

			for l=1:group_trials

				% replace with ssim over average images?

				cur_contour_day=double(tfdensity{i}.all{k}(:,idx,group1(l)));
				norm2=sum(cur_contour_day(:).^2);

				cur_contour_night=double(tfdensity{i}.all{k}(:,idx,group2(l)));
				norm3=sum(cur_contour_night(:).^2);
				%

				score_day2night{i}{j,k}(l)=sum(sum(template2.*cur_contour_day))/sqrt(norm12*norm2);
				score_day2day{i}{j,k}(l)=sum(sum(template1.*cur_contour_day))/sqrt(norm11*norm2);
				score_night2night{i}{j,k}(l)=sum(sum(template2.*cur_contour_night))/sqrt(norm12*norm3);

				%score_day2night{i}{j,k}(l)=corr(template2(:),cur_contour_day(:));
				%score_day2day{i}{j,k}(l)=corr(template1(:),cur_contour_day(:));
				%score_night2night{i}{j,k}(l)=corr(template2(:),cur_contour_night(:));

			end

			for l=1:ntrials
				cur_contour=double(tfdensity{i}.all{k}(:,idx,l));
				norm_both=sum(cur_contour(:).^2);
				score_all2all{i}{j,k}(l)=sum(sum(template3.*cur_contour))/sqrt(norm_both*norm13);
			end

		end
	end
end

% within_day=[];
% within_night=[];
% 
% for i=1:length(score_day2night)
% 	for j=1:ndays-1
% 		within_day=[within_day mean(score_day2night{i}{j,j})];
% 		within_night=[within_night median(score_night2night{i}{j,j+1})];
% 	end
% end
% 
% % lag analysis
% 
% between_day=[];
% between_night=[];
% 
% for i=1:length(score_day2night)
% 	for j=1:ndays-1
% 		between_day=[between_day median(score_day2night{i}{j,j+1})];
% 		between_night=[between_night median(score_night2night{i}{j,j+1})];
% 	end
% end

within_day=[];
between_day=[];

for i=1:length(score_day2night)
	within_day(i)=mean(cat(2,score_night2night{i}{diag(ones(ndays,1),0)==1}));
	between_day(i)=mean(cat(2,score_day2night{i}{diag(ones(ndays-1,1),1)==1}));
end

morning=[];
evening=[];

for i=1:length(score_day2night)
	for j=1:size(score_day2night{i},1)-1
		% add to morning window if hr less than 8
		tmp=datevec(tfdensity{i}.timestamps{j});
		idx=tmp(:,4)>=19;
		norm_mu=mean(score_all2all{i}{j,j}(idx));
		norm_std=std(score_all2all{i}{j,j}(idx));
		for k=j+1
			tmp=datevec(tfdensity{i}.timestamps{k});
			idx=tmp(:,4)<=10;
			morning=[morning mean((score_all2all{i}{j,k}(idx)))];
			idx=tmp(:,4)>=19;
			evening=[evening mean((score_all2all{i}{j,k}(idx)))];
		end
	end
end


[p,h,stats]=signrank(within_day,between_day,'tail','right');

% analyze using evenly-spaced time windows



if filewrite
	save(fullfile(dirs.agg_dir,dirs.datastore_dir,'sdi_scores.mat'),...
		'score_day2night','score_night2night','score_day2day','score_all2all','-v7.3');
	fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,'audio_overnight.txt'),'w+');
	fprintf(fid,'Between day v within day:  p=%e w=%g n=%i',p,stats.signedrank,length(within_day));
	fclose(fid);
end

%%

lags=[1:5];

for i=lags
	lag_day{i}=[];
	lag_night{i}=[];
	lag_all{i}=[];
end

for i=1:length(score_day2night)
	for j=lags
        for k=j:max(lags)
            
            lag=(k-j)+1;
            
            norm_mu=mean(score_night2night{i}{j,j});
            norm_std=std(score_night2night{i}{j,j});
            
            tmp=datevec(tfdensity{i}.timestamps{k});
			
            idx=tmp(:,4)<=9;
            idx2=tmp(:,4)>=19;
            
            lag_day{lag}=[lag_day{lag} mean(score_all2all{i}{j,k}(idx))];
            lag_night{lag}=[lag_night{lag} mean(score_all2all{i}{j,k}(idx2))];
            lag_all{lag}=[lag_all{lag} mean(score_all2all{i}{j,k})];
        
            % normalize??
            
        end
	end
end

lag_day_night=cell(1,length(lag_day)*2);

counter=1;
for i=1:2:length(lags)*2
	lag_day_night{i}=lag_day{counter};
	lag_day_night{i+1}=lag_night{counter};
	counter=counter+1;
end

% all day shifts and night shifts?

plot_mu=cellfun(@mean,lag_day_night);
plot_ci=zeros(2,length(plot_mu));

for i=1:length(plot_mu)
	%plot_ci(:,i)=bootci(1e3,{@mean,lag_day_night{i}},'type','cper');
	sem=std(lag_day_night{i})./sqrt(length(lag_day_night{i}));
	plot_ci(:,i)=[plot_mu(i)+sem;plot_mu(i)-sem];
end
