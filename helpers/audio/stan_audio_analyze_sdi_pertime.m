function [fig,fig_stats]=stan_audio_analyze_sdi_pertime
%
%
%

[options,dirs]=stan_preflight;
filewrite=true;

% load in each mat file, make spectral density image, etc.

% assume we have contours1

% get listing

load(fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis','sdi_analysis_data_consensus.mat'),'tfdensity');

%%
% analyze circadian variation within-day/overnight shift
thresh=-inf;
score={};
score_timediff={};


%
%% all to all???
all_score=[];
all_time=[];
all_id=[];


for i=1:length(tfdensity)

	ndays=length(tfdensity{i}.all);
	i
	%
	% 	for j=1:ndays
	%
	% 		template=mean(tfdensity{i}.all{j}(:,:,:),3);
	%
	% 		ntrials=size(tfdensity{i}.all{j},3);
	% 		idx=zscore(sum(template))>thresh;
	%
	% 		template=template(:,idx);
	%         norm1=sum(template(:).^2);
	%
	% 		for k=1:ntrials
	%
	%             cur_contour=tfdensity{i}.all{j}(:,idx,k);
	% 			norm2=sum(cur_contour(:).^2);
	%
	%  			score{i,j}(k)=sum(sum(template.*cur_contour))/sqrt(norm1*norm2);
	%             %score{i,j-1}(k)=corr(template(:),cur_contour(:));
	% 			score_timediff{i,j}(k)=tfdensity{i}.timestamps{j}(k)-floor(tfdensity{i}.timestamps{j}(k));
	%
	% 		end
	%
	% 	end

	for j=1:ndays

		% get bins, average all-to-all pairwise correlation within bin

		hrs=(tfdensity{i}.timestamps{j}-floor(tfdensity{i}.timestamps{j}))*24;
		%hrs=round(hrs);
		[~,hrs]=histc(hrs,[0:1:24]);
		bins=unique(hrs);

		for k=1:length(bins)
			% grab trials, get pdist, average upper triangle

			trial_idx=hrs==bins(k);
			use_trials=tfdensity{i}.all{j}(:,:,trial_idx);

			template=mean(use_trials,3);
			norm1=sum(template(:).^2);

			tmp=nan(size(use_trials,3),1);
			if size(use_trials,3)>2

				for l=1:size(use_trials,3)
					cur_contour=use_trials(:,:,l);
					norm2=sum(cur_contour(:).^2);
					tmp(l)=sum(sum(template.*cur_contour))./sqrt(norm1*norm2);
				end

				all_score=[all_score;mean(tmp)];
				all_time=[all_time;bins(k)];

			end
		end
	end
end

[~,idx]=sort(all_time);
fig.song_variability=figure();
stan_plot_regress(all_time,all_score);
colormap(bone);
xlim([6 23])
%ylabel('Song consistency')
%xlabel('Time (hr)')
set(gca,'TickLength',[0 0],'YTick',[.2 .9],'FontSize',7);
title('')
ylim([.2 .9])
set(fig.song_variability,'position',[300 700 230 210]);
[fig_stats.songvar.r,fig_stats.songvar.p]=corr(all_time,all_score,'type','spearman');

if filewrite
  fid=fopen(fullfile(dirs.agg_dir,dirs.stats_dir,'fig5_songvar.txt'),'w+');
  fprintf(fid,'Song variability v time: r=%g p=%e',fig_stats.songvar.r,fig_stats.songvar.p);
  fclose(fid);
end


%%
%
% all_score=[];
% all_time=[];
% all_id=[];
%
% for i=1:size(score,1)
%
% 	mu=mean(cat(2,score{i,:}));
% 	sig=std(cat(2,score{i,:}));
% 	peak=max(cat(2,score{i,:}));
% 	trough=min(cat(2,score{i,:}));
%
% 	for j=1:size(score,2)
%
% 		if ~isempty(score{i,j})
%
% 			all_score=[all_score;score{i,j}(:)];
% 			all_time=[all_time;score_timediff{i,j}(:)];
% 			all_id=[all_id;i];
%
% 		end
% 	end
% end
