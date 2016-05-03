%%
% collect the pns from pbio paper and check for drift in
% putative x projector activity

[options,dirs]=stan_preflight;
load(fullfile(dirs.agg_dir,dirs.datastore_dir,'pbio_pns.mat'));

%% get the per-trial peak values...

win=.0025;
sliding_win=15;
sliding_overlap=12;
peak_changes=[];
colors={'r','b','g','y'};

for i=1:length(store_spikes)
    
    % get the average rate
    proc_fs=store_spikes(i).fs;
    mu=mean(store_spikes(i).smoothed_spikes*proc_fs,2);
    [~,locs]=findpeaks(mu,'minpeakheight',100);
    npeaks=length(locs);
    
    [nsamples,ntrials]=size(store_spikes(i).smoothed_spikes);
    
    if npeaks>1
       
        % now for each peak, grab peak height in a small window
        % does this peak systematically change over the course of the day?
        
        peak_changes(end+1).peak_vals=nan(ntrials,npeaks);
        peak_changes(end).peak_locs=nan(ntrials,npeaks);
        peak_changes(end).peak_locs_mu=locs;
        peak_changes(end).cell_id=i;
        
        for j=1:npeaks
            
            % grab n msec to the left and right of the peak
            
            left_edge=round(((locs(j)/proc_fs)-win)*proc_fs);
            right_edge=round(((locs(j)/proc_fs)+win)*proc_fs);
            
            if left_edge<1 | right_edge>nsamples
                continue;
            end
            
            slice_data=zscore(store_spikes(i).smoothed_spikes);
            mu_slice=slice_data(left_edge:right_edge,:);
            [tmp_vals]=max(mu_slice);
            
            peak_changes(end).peak_vals(:,j)=tmp_vals;
            %peak_changes(end).peak_locs(:,j)=tmp_locs;
            
        end        
    end
end

%%

for i=1:length(peak_changes)
    
    [ntrials,npeaks]=size(peak_changes(i).peak_vals);
    
    % first and second half???
    
    steps=1:sliding_win-sliding_overlap:ntrials-sliding_win;
    nsteps=length(steps);
    smooth_mu=zeros(npeaks,nsteps);
    smooth_ci=cell(npeaks);
    smooth_x=zeros(1,nsteps);
    
    for j=1:nsteps
        
        points=steps(j):steps(j)+sliding_win;
        smooth_x(j)=mean(points);
        use_points=peak_changes(i).peak_vals(points,:);
        smooth_mu(:,j)=mean(use_points);
    
        for k=1:npeaks
           smooth_ci{k}(:,j)=bootci(1000,{@mean,use_points(:,k)},'type','per','alpha',.05);
        end
        
    end
    
   figs.([sprintf('cell%i',peak_changes(i).cell_id)])=figure();
   
   for j=1:npeaks
    markolab_shadeplot(smooth_x,smooth_ci{j},colors{j},'none');
    hold on;
    plot(smooth_x,smooth_mu(j,:),'k-');
   end
   
   title([sprintf('Cell %i',peak_changes(i).cell_id)]);
   ylabel('FR (Z)');
   xlabel('Trial');
   set(gca,'TickDir','out');
   box off;
    
    % plot with 95% CI, slide a window across the data...
    
end


% save the figures

%%

cells=fieldnames(figs);

for i=1:length(cells)
   
   % change width and save??
   
   set(figs.(cells{i}),'paperpositionmode','auto','position',[107   359   408   283]);
   markolab_multi_fig_save(figs.(cells{i}),pwd,[ cells{i} '_peakheights' ] ,'eps,png,fig');
  
end