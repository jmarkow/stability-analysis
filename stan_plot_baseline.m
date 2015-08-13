function stan_plot_baseline(BASELINE)
%
%
%
%
%

ci_inv=.01; % alpha level (/2 to get two-tailed)
ylimit_rounding=1e-1;
r_cutoff=.4; % below this value typically due to equipment failure
save_name='baseline_regression';
options_name='options.txt';
dirs_name='dirs.txt';

cur_file=mfilename('fullpath');
[cur_path,~,~]=fileparts(cur_file);
options=stan_read_options(fullfile(cur_path,options_name));
dirs=stan_read_options(fullfile(cur_path,dirs_name));



% first column is days since 1, second is regression value, third is bird ID


% remove points where x=0 (by definition == 1, artifacts appear to have r<.4)

first_day=BASELINE(:,1)==0;
artifacts=BASELINE(:,2)<=r_cutoff;
artifacts2=isnan(BASELINE(:,2));

to_del=first_day|artifacts|artifacts2
BASELINE(to_del,:)=[];

% regress, compute 95% CI


b=regress(BASELINE(:,2),[ones(size(BASELINE(:,1))) BASELINE(:,1)]);
npoints=length(BASELINE(:,1));

pred_x=[min(BASELINE(:,1)):.5:max(BASELINE(:,1))];
pred_y=b(1)+pred_x*b(2);

% first get the residuals

estimate=b(1)+BASELINE(:,1)*b(2);
res=BASELINE(:,2)-estimate;
sse=sum(res.^2);

t_stat=tinv((1-ci_inv/2),npoints-2);
pred_x_mu=mean(pred_x);

% get the confidence interval

conf=t_stat*sqrt(((1/(npoints-2))*sse)*...
	((1/npoints)+((pred_x-pred_x_mu).^2)/sum((BASELINE(:,1)-pred_x_mu).^2)));

% interval is predicted values +/- confidence 


fig=figure();
%markolab_shadeplot(pred_x(:)',[pred_y(:)'+conf;pred_y(:)'-conf],[.8 .8 .8],'none');
scatter(BASELINE(:,1),BASELINE(:,2),30,BASELINE(:,3));
hold on;
plot(pred_x,pred_y,'r-','linewidth',1.5);
plot(pred_x,pred_y+conf,'r--','linewidth',1.5);
plot(pred_x,pred_y-conf,'r--','linewidth',1.5);

ylimits=ylim()
ylimits=round(ylimits/ylimit_rounding)*ylimit_rounding;
ylim(ylimits);

set(gca,'TickDir','out','TickLength',[0 0],'ytick',[ylimits(1):.1:ylimits(2)]);
xlim([min(pred_x)-5 max(pred_x)+5]);
set(gca,'xtick',[0:20:max(pred_x)]);
ylabel('Rate correlation (R)');
xlabel('Time (days)');

set(fig,'position',[200 200 220 190],'paperpositionmode','auto');
markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),save_name,'eps,png,fig');

