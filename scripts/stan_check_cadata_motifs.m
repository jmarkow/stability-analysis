% script for analyzing new calcium data
% assumes data already loaded in from ROI_data_cleansed.mat

%%

% only for lw76

bird_name='lw76';
movie_fs=100;
motif_selection=[2];
plot_data=true;
%load custom_colormaps;

%% arrange the data into a simpler format

if strcmpi(bird_name,'lw76');
  roi_data=cell(1,length(ROI_dat));
  for i=1:length(ROI_dat)
    roi_data{i}=cat(3,ROI_dat{i}.align_detrended{:});
  end
  movie_fs=22;
end

%%

for ii=1:length(motif_selection)
    if ~strcmpi(bird_name,'lw76')

        [roi_data,roi_dates,roi_times,roi_motifs,roi_filenames]=...
            stan_cadata_collect_freedomscope_v2_cell(ROI_data_cleansed);

        if ~isempty(motif_selection)
            for i=1:length(roi_data)
                roi_data_motifs{ii}{i}=roi_data{i}(:,:,roi_motifs{i}==motif_selection(ii));
            end
        end
    else
        roi_data_motifs{ii}=roi_data;
    end

    figs.schnitzer=figure('position',[400 400 600 300],'paperpositionmode','auto');
    [ave_mat{ii},inc_rois{ii}]=stan_cadata_sortmat(roi_data_motifs{ii},'scaling','l','sort_day',1,'smoothing',0,'smooth_kernel','g',...
        'padding',[.25 .75],'movie_fs',movie_fs,'chk_day',1,'fig_row',1,'fig_nrows',1,'realign',1);

end

%%
%stan_cadata_sortmat(roi_data,'scaling','l','sort_day',3,'smoothing',0,'smooth_kernel','g',...
%  'padding',[.25 0.75],'movie_fs',movie_fs,'chk_day',1,'fig_row',2,'fig_nrows',2,'realign',0);

ax=findall(figs.schnitzer,'type','axes')
linkaxes(ax,'xy')
ylimits=ylim();
set(ax(end),'YTick',ylim(),'YTickLabel',[ylim()-min(ylim())]+1,'FontSize',12);
yh=ylabel(ax(end),'Cell');
xlimits=xlim();
h=line([xlimits(1) xlimits(1)+.1],[ylimits(2)+3 ylimits(2)+3],'color','k','parent',ax(1))
set(h,'clipping','off');

load custom_colormaps;
colormap(calcium_contrast);

%%
% get com for each day and plot to check for systematic shifts (indicative
% of a/v sync issues)

ndays=length(ave_mat{1});

[nsamples,nrois]=size(ave_mat{1}{1});
ind=repmat([1:nsamples]',[1 nrois]);
com=zeros(nrois,ndays);

% get the sort indices

[~,loc]=max(ave_mat{1}{1}); % get max
[~,loc]=sort(loc); % sort max ascending

for i=1:ndays
    %[~,loc]=max(ave_mat{i});
    den=sum(ave_mat{1}{i}(:,loc));
    com(:,i)=sum(ind.*ave_mat{1}{i}(:,loc))./den;
    %com(:,i)=loc/movie_fs;
end

if plot_data
    
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
end
