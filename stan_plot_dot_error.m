function stan_dot_error_plot(X,Y,CI,C,varargin)
% takes collected stats, plots and performs hypothesis tests
% 
%
%


if nargin<3 | isempty(C)
	C=ones(size(MU));
end

colors='';
markersize=20;
symbols={'o','x','s','d','*','+','p','h'};
nparams=length(varargin);
line_alpha=.5;

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'colors'
			colors=varargin{i+1};
		case 'markersize'
			markersize=varargin{i+1};
		case 'smooth_line'
			smooth_line=varargin{i+1};
		case 'line_alpha'
			line_alpha=varargin{i+1};
	end
end

grps=unique(C);
ngrps=length(grps);
cmap=parula(ngrps);

for i=1:length(X)

	% plot confidence intervals
	
	h2=plot([X(i) X(i)],[CI(1,i) CI(2,i)],'k-','color',cmap(C(i),:),'linewidth',1);
	col=h2.Color;
	h2.Color=[ col line_alpha ];
	hold on;

end

h=scatter(X,Y,markersize,C,'markerfacecolor','flat');

box on;
set(gca,'TickDir','out','TickLength',[0 0]);

