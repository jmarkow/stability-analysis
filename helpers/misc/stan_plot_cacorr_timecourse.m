function FIG=stan_plot_corr_timecourse(STATS,COLORS,NBOOTS)
%
%
% defaults to 95%

interp_factor=8;
FIG=figure();

for i=1:length(STATS)

  cur_data=STATS(i).rmat_mu.lag.all;

  x=find(cellfun(@length,cur_data)>0);

  mu=ones(1,length(x));
  mu_ci=ones(2,length(x));

  % day one from bootstrap

  for j=1:length(x)
    mu_ci(:,j)=bootci(NBOOTS,{@mean,cur_data{x(j)}(:)},'type','cper','alpha',.01);
    mu(:,j)=mean(cur_data{x(j)}(:));
  end

  x=x-1;
  xx=x(1):1/interp_factor:x(end);

  %xx=0:1/interp_factor:length(x);

  mu_interp=interp1(x(:),mu(:),xx(:),'spline');
  mu_ci_interp=interp1(x(:),mu_ci',xx(:),'spline');
  %
  markolab_shadeplot(xx,mu_ci_interp',COLORS(i,:),'k');
  hold on;
  plot(xx,mu_interp,'k--');

  %markolab_shadeplot(x,mu_ci,COLORS(i,:),'k');
  %hold on;
  %plot(x,mu,'k--');
  plot(x,mu,'ko','markerfacecolor','w');

end

yh=ylabel('Correlation (R)')
set(yh,'position',get(yh,'position')+[.16 0 0]);
ylim([.5 1]);
ylimits=ylim();
set(gca,'YTick',ylimits,'TickLength',[0 0],'FontSize',7);
