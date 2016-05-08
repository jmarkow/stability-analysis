%%%% plots the lhp33 ephys example
% assumes roboaggregate.mat already loaded in...

[b,a]=ellip(3,.2,40,[300]/(ephys.fs/2),'high');

% high pass the data

ephys_filtered=filtfilt(b,a,double(ephys.data(:,:,13)));

% get rms

rms_tau=.01;
rms_smps=round(rms_tau.*ephys.fs);
rms_filter=ones(rms_smps,1)./rms_smps;

ephys_use=ephys_filtered(:,motif_number==1);
ephys_rms=(filter(rms_filter,1,ephys_use.^2));

trials_to_use=[80 140];

load custom_colormaps;

figs.msv=figure();
imagesc(ephys.t,[],ephys_rms');
xlabel('Time (s)');
ylabel('Trial');
hold on;
h=scatter(ones(1,length(trials_to_use)).*-.05,trials_to_use,'>','filled');
set(h,'Clipping','off');
set(gca,'Clipping','off','ytick',[1 size(ephys_rms,2)],'xtick',[0 .8]);
colormap(fee_map);
c=colorbar();
c.Label.String='Mean Squared Voltage (uV)';
xlim([0 .8]);
set(figs.msv,'paperpositionmode','auto','position',[500 500 246 180]);


figs.singletrials=figure();
nplots=length(trials_to_use);
ax=[];
for i=1:length(trials_to_use)
  ax(i)=subplot(nplots,1,i);
  plot(ephys.t,ephys_use(:,trials_to_use(i)),'k-');
  box off;
  %title([sprintf('Trial %i',trials_to_use(i))]);
end
linkaxes(ax,'xy');
xlim([0 .8]);
set(ax(1),'XTick',[])
set(ax(2),'xtick',[0 .8]);
xlabel('Time (s)');
ylabel('V (uV)');
set(figs.singletrials,'paperpositionmode','auto','position',[500 500 165 180]);
%%

% figure();
% 
% motifs=unique(motif_number);
% motifs(motifs>4)=[];
% 
% ax=[];
% 
% for i=1:length(motifs)
%   ax(i)=subplot(length(motifs),1,i);
%   imagesc(ephys.t,[],(filter(rms_filter,1,ephys_filtered(:,motif_number==motifs(i)).^2))');
%   title(sprintf('Motif %i',motifs(i)));
% end
% 
% linkaxes(ax,'x');
% set(ax(1:end-1),'xtick',[]);
% xlabel('Time (s)');
% ylabel('Trial');

%%
fignames=fieldnames(figs);
for i=1:length(fignames)
    markolab_multi_fig_save(figs.(fignames{i}),pwd,[ 'lhp33_' fignames{i} ],'eps,png,fig','renderer','painters');
end