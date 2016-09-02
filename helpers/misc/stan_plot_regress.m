function b=stan_plot_baseline(X,Y,C,varargin)
%
%
%
%
%

if nargin<3 | isempty(C)
	C=ones(size(X));
end

r_cutoff=.4; % below this value typically due to equipment failure
x_int=.01;
ylimit_rounding=1e-1;
save_name='baseline_regression';
[options,dirs]=stan_preflight;
ci_inv=.01;
shaded_conf=1;
linewidth=.5;
markersize=30;
clip=[];
plot_mode='scatter';
nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'r_cutoff'
			r_cutoff=varargin{i+1};
		case 'x_int'
			x_int=varargin{i+1};
		case 'ylimit_rounding'
			ylimit_rounding=varargin{i+1};
		case 'save_name'
			save_name=varargin{i+1};
		case 'shaded_conf'
			shaded_conf=varargin{i+1};
		case 'markersize'
			markersize=varargin{i+1};
        case 'clip'
            clip=varargin{i+1};
        case 'plot_mode'
            plot_mode=varargin{i+1};
            
	end
end

% first column is days since 1, second is regression value, third is bird ID

% remove points where x=0 (by definition == 1, artifacts appear to have r<.4)

[b]=regress(Y,[ones(size(X)) X]);
npoints=length(X);

pred_x=[min(X):x_int:max(X)];
pred_y=b(1)+pred_x*b(2);

% first get the residuals

b
estimate=b(1)+X*b(2);
res=Y-estimate;
sse=sum(res.^2);
mse=sse/(npoints-2);
t_stat=tinv((1-ci_inv/2),npoints-2)
pred_x_mu=mean(pred_x);

% get the confidence interval

%conf=t_stat*sqrt(((1/(npoints-2))*sse)*...
%	((1/npoints)+((pred_x-pred_x_mu).^2)/sum((X-pred_x_mu).^2)));
conf=t_stat*sqrt((mse*((1/npoints)+((pred_x-pred_x_mu).^2)/(sum((X-pred_x_mu).^2)))));
%conf_low = b(1)+bint(2,1)*pred_x;

% interval is predicted values +/- confidence 

%fig=figure();

ngrps=length(unique(C));

if ~isempty(clip)
    tmp=find(Y>clip);
    Y_clip=Y(tmp);
    Y_clip=ones(size(Y_clip))*(clip+.1*clip);
    X_clip=X(tmp)+.015*range(X)*randn(size(tmp));
    X(tmp)=[];
    Y(tmp)=[];
    C_clip=C(tmp);
    C(tmp)=[];
end

h=scatter(X,Y,markersize,C,'markerfacecolor','flat');
hold on;

if ~isempty(clip)
    h2=scatter(X_clip,Y_clip,markersize*.5,C_clip,'markerfacecolor','none');
end

colormap(parula(ngrps));

if shaded_conf
    markolab_shadeplot(pred_x(:)',[pred_y(:)'+conf;pred_y(:)'-conf],[.85 .85 .85],'none');
    %markolab_shadeplot(pred_x(:)',[conf_hi;conf_low],[.85 .85 .85],'none');
	hold on;
end

if ~shaded_conf
 	plot(pred_x,pred_y+conf,'r--','linewidth',linewidth);
 	plot(pred_x,pred_y-conf,'r--','linewidth',linewidth);
    %plot(pred_x,conf_hi,'r--','linewidth',linewidth);
   	%plot(pred_x,conf_low,'r--','linewidth',linewidth);
end

plot(pred_x,pred_y,'r-','linewidth',linewidth);

ylimits=ylim()
ylimits=round(ylimits/ylimit_rounding)*ylimit_rounding;
ylim(ylimits);
set(gca,'layer','top');
title([num2str(b)])
%set(gca,'TickDir','out','TickLength',[0 0],'ytick',[ylimits(1):.1:ylimits(2)]);
%xlim([min(pred_x)-5 max(pred_x)+5]);
%set(gca,'xtick',[0:20:max(pred_x)]);
%ylabel(y_label);
%xlabel(x_label);

%set(fig,'position',[200 200 220 190],'paperpositionmode','auto');
%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),save_name,'eps,png,fig');

