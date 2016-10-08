function stan_fig2a(BTLFP,LFP_FS)
% Generates Fig. 2a, burst-triggered LFPs for HVC projection neurons, interneurons, and multi-unit

if nargin<2 | isempty(LFP_FS)
	LFP_FS=1e3;
end

fig=figure();

nsamples=size(BTLFP.int.lfp_win,1);
xwin=floor(nsamples/2);
win_t=[-xwin:xwin]/LFP_FS;
win_t=win_t*1e3;

% find neurons with >1 burst per song (putative X projector)

peakmask=BTLFP.pn.peakcount>1; % pn
win_id=BTLFP.pn.cell_id(peakmask);
inc_trials=stan_retain_trials(win_id);

ax(1)=subplot(3,1,1);
peak_idx=find(peakmask);
mu=mean(BTLFP.pn.lfp_win(:,peak_idx(inc_trials)),2);
lim=max(abs(mu));
lim_rnd=ceil(lim*10)/10;

plot(win_t,mu,'b-');
ylim([-lim_rnd lim_rnd]);
ylimits=ylim();
hold on;

plot([0 0],[ylimits],'k-');
set(gca,'YTick',ylimits,'xtick',[],'FontSize',7);
yh=ylabel('Amp. (Z)');
set(yh,'position',get(yh,'position')+[-.2 0 0]);
%title('HVC_X');

% now get the interneurons

win_id=BTLFP.int.cell_id;
inc_trials=stan_retain_trials(win_id);

% window center

% 2 x 2 grid, mean waveform on top, aligned PLI belo

ax(2)=subplot(3,1,2);

mu=mean(BTLFP.int.lfp_win(:,inc_trials),2);
lim=max(abs(mu));
lim_rnd=ceil(lim*10)/10;

plot(win_t,mu,'b-');
ylim([-lim_rnd lim_rnd]);
ylimits=ylim();
hold on;

plot([0 0],[ylimits],'k-');
set(gca,'YTick',ylimits,'xtick',[],'FontSize',7);
%title('Interneurons');

inc_trials=stan_retain_trials(BTLFP.mu.peak_id);

% finally multi-unit

ax(3)=subplot(3,1,3);
mu=mean(BTLFP.mu.lfp_win(:,inc_trials),2);
lim=max(abs(mu));
lim_rnd=ceil(lim*10)/10;

plot(win_t,mu,'b-');
ylim([-lim_rnd lim_rnd]);
ylimits=ylim();

hold on;
plot([0 0],[ylimits],'k-');
set(gca,'YTick',ylimits,'xtick',[-100 0 100],'FontSize',7);
%title('Multi-unit');

xlabel('Time (ms)');
linkaxes(ax,'x');

set(fig,'units','centimeters','position',[4 4 3 7],'paperpositionmode','auto');
