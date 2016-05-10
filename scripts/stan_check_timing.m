%%%
% quick check of calcium data

roi_times{1}=roi_times{1}-min(roi_times{1});
left_edge=max(find(roi_times{1}<.2));
right_edge=min(find(roi_times{1}>max(roi_times{1})-.7));

motif_selection=2;

mu1=mean(zscore(roi_data{1}(left_edge:right_edge,:,roi_motifs{1}==motif_selection)),3);
[~,idx]=max(mu1);
[~,idx2]=sort(idx);

figure();
ax=[];
mu={};
for i=1:length(roi_data)
  ax(i)=subplot(length(roi_data),1,i);
  mu{i}=mean(zscore(roi_data{i}(left_edge:right_edge,:,roi_motifs{i}==motif_selection)),3);
  imagesc(filter(ones(10,1)/10,1,zscore(mu{i}(:,idx2)))');
  caxis([0 2]);
end
linkaxes(ax,'x');
