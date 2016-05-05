%%%% plots the lhp33 ephys example
% assumes roboaggregate.mat already loaded in...

[b,a]=ellip(3,.2,40,[300]/(ephys.fs/2),'high');

% high pass the data

ephys_filtered=filtfilt(b,a,double(ephys.data(:,:,12)));

% get rms

rms_tau=.01;
rms_smps=round(rms_tau.*ephys.fs);
rms_filter=ones(rms_smps,1)./rms_smps;

ephys_use=ephys_filtered(:,motif_number==1);
ephys_rms=filter(rms_filter,1,ephys_use.^2);


%%

trials_to_use=[57 59];

load custom_colormaps;

figure();imagesc(ephys.t,[],ephys_rms');
xlabel('Time (s)');
ylabel('Trial');
colormap(fee_map);
c=colorbar();
c.Label.String='MS (uV^2)';

figure();
nplots=length(trials_to_use);
ax=[];
for i=1:length(trials_to_use)
  ax(i)=subplot(nplots,1,i);
  plot(ephys.t,ephys_use(:,trials_to_use(i)));
  box off;
  title([sprintf('Trial %i',trials_to_use(i))]);
end
linkaxes(ax,'xy');
set(ax(1),'XTick',[])
xlabel('Time (s)');
ylabel('V (uV)');

%%

figure();

motifs=unique(motif_number);
motifs(motifs>4)=[];
motifs

ax=[];

for i=1:length(motifs)
  ax(i)=subplot(length(motifs),1,i);
  imagesc(ephys.t,[],(filter(rms_filter,1,ephys_filtered(:,motif_number==motifs(i)).^2))');
  title(sprintf('Motif %i',motifs(i)));
end

linkaxes(ax,'x');
set(ax(1:end-1),'xtick',[]);
xlabel('Time (s)');
ylabel('Trial');
