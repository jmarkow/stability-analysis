%%%
% quick check of calcium data

mu1=mean(zscore(roi_data{1}),3);

[~,idx]=max(mu1);
[~,idx2]=sort(idx);

figure();
ax=[];
for i=1:length(roi_data)
  ax(i)=subplot(length(roi_data),1,i);
  mu=mean(zscore(roi_data{i}),3);
  imagesc(mu(:,idx2)');
end
linkaxes(ax,'x');
