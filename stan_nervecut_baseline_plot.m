function stan_nervecut_baseline_plot(BASELINE,NERVECUT)
% takes collected stats, plots and performs hypothesis tests
% 
%
%



dirs_name='dirs.txt';
save_name='nervecut_acoustic_features';
r_cutoff=.4;
% get options

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);
dirs=stan_read_options(fullfile(cur_path,dirs_name));


% bin the baseline data, form 95% confidence interval to compare with nervecut data

win_steps=[0 15 30 60 100];
win_steps(1)=-inf;
win_steps(end)=inf;

first_day=BASELINE(:,1)==0;
artifacts=BASELINE(:,2)<=r_cutoff;
artifacts2=isnan(BASELINE(:,2));

to_del=first_day|artifacts|artifacts2
BASELINE(to_del,:)=[];

[~,bins]=histc(BASELINE(:,1),win_steps);

figure();
for i=1:length(win_steps)-1

	cur_bin=[win_steps(i) win_steps(i+1)];
	cur_bindata=BASELINE(bins==i,2);
	cur_bindata(isnan(cur_bindata))=[];
	
	bin_mu=median(cur_bindata);
	ci=bootci(1e3,{@median,cur_bindata});

	bin_conf(i)=prctile(cur_bindata,.005)
	cur_bin(cur_bin==-inf)=0;
	cur_bin(cur_bin==inf)=150;

	markolab_shadeplot([cur_bin],[1 1;bin_conf(i) bin_conf(i)],[.95 .95 .95],'none');
	hold on;
	plot([cur_bin],[bin_conf(i) bin_conf(i)],'r-','linewidth',1.5)

end

for i=1:length(win_steps)-2
	plot([win_steps(i+1) win_steps(i+1)],[bin_conf(i) bin_conf(i+1)],'r-','linewidth',1.5)
end

h=scatter(NERVECUT(:,1),NERVECUT(:,2),40,NERVECUT(:,5),'markerfacecolor','flat');
colors=get(h,'cdata')
xdata=get(h,'xdata');
cmap=parula(length(unique(colors)));

for i=1:length(xdata)

	% plot confidence intervals
	
	plot([xdata(i) xdata(i)],[NERVECUT(i,3) NERVECUT(i,4)],'k-','color',cmap(colors(i),:),'linewidth',1.5);

end

box on;
set(gca,'TickDir','out','TickLength',[0 0]);
