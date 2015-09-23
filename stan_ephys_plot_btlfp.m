function fig=stan_btlfp_plot(BTLFP)
%
%
%
%
%
%

[options,dirs]=axcorr_preflight;
nbootstraps=100;

% function for converting to Rayleigh Z
% convert to Rayleigh Z or some unbiased estimator

% convert hil windows to PLI, then convert to unbiased PLI estimator w/ CIs

% first get cell ids, grab same number of cell ids

nsamples=size(BTLFP.int.hil_win,1);
xwin=floor(nsamples/2);
win_t=[-xwin:xwin]/options.lfp_fs;
win_t=win_t*1e3;

peakmask=BTLFP.pn.peakcount>1; % pn
win=BTLFP.pn.hil_win(:,peakmask);
win_id=BTLFP.pn.cell_id(peakmask);
inc_trials=axcorr_retain_trials(win_id);

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

win=BTLFP.int.hil_win;
win_id=BTLFP.int.cell_id;
inc_trials=axcorr_retain_trials(win_id);

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

inc_trials=axcorr_retain_trials(BTLFP.mu.peak_id);

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
