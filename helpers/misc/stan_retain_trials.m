function RETAIN_TRIALS=axcorr_retain_trials(ID)

[uniq_id,~,id_idx]=unique(ID);
nwin=length(uniq_id);
id_idx_count=zeros(1,nwin);

for i=1:nwin
	id_idx_count(i)=sum(id_idx==i);
end

% retain this many trials from each cell

min_trials=min(id_idx_count);
RETAIN_TRIALS=[];

for i=1:nwin
	tmp=find(id_idx==i);
	idx=1:length(tmp);
	RETAIN_TRIALS=[RETAIN_TRIALS;tmp(idx<=min_trials)];
end
