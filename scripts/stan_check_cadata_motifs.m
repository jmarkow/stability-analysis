% script for analyzing new calcium data
% assumes data already loaded in from ROI_data_cleansed.mat


% only for lw76

tokens=regexp(pwd,filesep,'split');
bird_name=lower(tokens{end});
movie_fs=100;
motif_selection=[2];
plot_data=true;
%load custom_colormaps;

%% arrange the data into a simpler format

if strcmpi(bird_name,'lw76');
  roi_data=cell(1,length(ROI_A));
  for i=1:length(ROI_A)
    roi_data{i}=cat(3,ROI_A{i}.align_detrended{:});
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
    [ave_mat{ii},inc_rois{ii},sorting_idx{ii}]=stan_cadata_sortmat(roi_data_motifs{ii},'scaling','l','sort_day',1,'smoothing',0,'smooth_kernel','g',...
        'padding',[.3 .65],'movie_fs',movie_fs,'chk_day',1,'fig_row',1,'fig_nrows',1,'realign',1);

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
h=line([xlimits(1) xlimits(1)+.1],[ylimits(2)+1 ylimits(2)+1],'color','k','parent',ax(1))
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


for i=1:ndays
    %[~,loc]=max(ave_mat{i});
    den=sum(ave_mat{1}{i}(:,sorting_idx{1}));
    com(:,i)=sum(ind.*ave_mat{1}{i}(:,sorting_idx{1}))./den;
    %com(:,i)=loc/movie_fs;
end

tdiff=diff(com/movie_fs,[],2);

%%

% now get peak times using the method from the paper


[peak.times{1} peak.vals{1}]=fb_compute_peak_simple(ave_mat{1}{1}(:,sorting_idx{1}),...
	'thresh_t',.1,'debug',0,'onset_only',0,'thresh_hi',.5,'thresh_int',5,'thresh_dist',.2,...
	'fs',movie_fs); % thresh_int previously

for i=2:ndays
   [peak.times{i} peak.vals{i}]=fb_compute_peak_simple(ave_mat{1}{i}(:,sorting_idx{1}),...
	'thresh_t',.1,'debug',0,'onset_only',0,'thresh_hi',.5,'thresh_int',5,'thresh_dist',.2,...
	'fs',movie_fs); % thresh_int previously
end

%%
tdiff_peaks=nan(ndays-1,nrois);

for i=2:ndays
    for j=1:nrois

        mindist=inf;

        for k=1:length(peak.times{i}{j})

            tmp2=min(abs(peak.times{i}{j}(k)-peak.times{1}{j}))/movie_fs;

            if tmp2<mindist
                mindist=tmp2;
            end

            tdiff_peaks(i-1,j)=mindist;

        end
    end
end

tdiff_peaks(tdiff_peaks==inf)=nan;
tdiff_peaks=tdiff_peaks';
