function fig=stan_ephys_plot_phasehist(PHASEHIST)
%
%
%
%
%
%

[options,dirs]=stan_preflight;

cell_types={'pn','int','mu'};
labels={'HVC_X','Interneurons','Multi-unit'};

%counter=2;

%subplot(7,1,1);
%xpoints=-pi:pi/32:3*pi;
%plot(cos(xpoints),'k','linewidth',1.5);
%axis tight;
%axis off;

for i=1:length(cell_types)

	PHASEHIST.(cell_types{i}).mag=repmat(PHASEHIST.(cell_types{i}).mag,[1 2]);
	PHASEHIST.(cell_types{i}).bins=[PHASEHIST.(cell_types{i}).bins PHASEHIST.(cell_types{i}).bins+2*pi];

	%ax(1)=subplot(7,1,counter:counter+1);
	
	ax(i)=subplot(length(cell_types),1,i);
	
	switch lower(cell_types{i}(1))
		case 'p'
			facecolor='r';
			edgecolor='k';
		case 'i'
			facecolor='b';
			edgecolor='k';
		case 'm'
			facecolor=[.8 .8 .8];
			edgecolor='k';
	end

	markolab_stairplot(PHASEHIST.(cell_types{i}).mag./max(PHASEHIST.(cell_types{i}).mag),...
		PHASEHIST.(cell_types{i}).bins,'method','a','facecolor',facecolor,'edgecolor','none');
	set(gca,'XTick',[-pi:pi:3*pi],'YTick',[0 1],'TickDir','out','TickLength',[0 0],'FontSize',7,...
		'XTickLabel',{'-\pi','0','\pi','2\pi','3\pi'});

	xlim([-pi 3*pi]);
	ylim([0 1]);

	if i==1
		yh=ylabel('LFP pow.');
		set(yh,'position',get(yh,'position')+[.2 0 0]);
	end

	if i==length(cell_types)
		xlabel('LFP phase (radians)');
	else
		set(gca,'XTick',[]);
	end

	set(gca,'layer','top')
	%counter=counter+2;

end

linkaxes(ax,'xy');
