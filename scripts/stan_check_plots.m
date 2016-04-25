%%
legend_string=cell(1,ndays);
for i=1:ndays
    legend_string{i}=sprintf('Day %i',i);
end

legend_string_df=cell(1,ndays-1);
for i=1:ndays-1
    legend_string_df{i}=sprintf('Day %i-%i',i,i+1);
end

% plot com for each day overlaid

figs.com=figure('position',[400 400 300 300],'paperpositionmode','auto');
h=plot(com./movie_fs,[1:nrois]);
set(gca,'ydir','rev','ticklength',[0 0]);
box off;
L=legend(h,legend_string);
xlabel('Time (s)');
ylabel('ROI');

figs.com_residue=figure('position',[400 400 300 300],'paperpositionmode','auto');
h=plot(diff(com,[],2)./movie_fs,[1:nrois]);
set(gca,'ydir','rev','ticklength',[0 0]);
box off;
L=legend(h,legend_string_df);
xlabel('Time (s)');
ylabel('ROI');

all_figs=fieldnames(figs);
for i=1:length(all_figs)
    %markolab_multi_fig_save(figs.(all_figs{i}),pwd,[ bird_name '_' all_figs{i}],'eps,png,fig');
end
