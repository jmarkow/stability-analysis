%
%
%

[options,dirs]=stan_preflight;
filewrite=false;

% load in each mat file, make spectral density image, etc.

% assume we have contours1

% get listing

load(fullfile(dirs.agg_dir,dirs.sdi_dir,'analysis','sdi_analysis_data_consensus.mat'),'tfdensity');

%%
% analyze circadian variation within-day/overnight shift
thresh=-inf;
score={};
score_timediff={};

for i=1:length(tfdensity)

	ndays=length(tfdensity{i}.all);
	i

	for j=2:ndays

		template=mean(tfdensity{i}.all{j-1}(:,:,:),3);
		
		ntrials=size(tfdensity{i}.all{j},3);
		idx=zscore(sum(template))>thresh;

		template=template(:,idx);
        norm1=sum(template(:).^2);

		for k=1:ntrials
			
            cur_contour=tfdensity{i}.all{j}(:,idx,k);
			norm2=sum(cur_contour(:).^2);
            
			score{i,j-1}(k)=sum(sum(template.*cur_contour))/sqrt(norm1*norm2);
			score_timediff{i,j-1}(k)=tfdensity{i}.timestamps{j}(k)-floor(tfdensity{i}.timestamps{j}(k));
            
		end

	end

end


%%

all_score=[];
all_time=[];
all_id=[];

for i=1:size(score,1)
	for j=1:size(score,2)
		if ~isempty(score{i,j})
            
			all_score=[all_score;score{i,j}(:)];
			all_time=[all_time;score_timediff{i,j}(:)];
			all_id=[all_id;ones(size(score{i,j}(:)))*i];
%             
		end
	end
end
