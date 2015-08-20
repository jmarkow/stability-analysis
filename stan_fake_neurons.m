function fig=stan_fake_neurons(COLORS)
%
%
%
%
%

if nargin<1
	COLORS='jet';
end

% for figure 1 schematic

[options,dirs]=stan_preflight;

pulse_t=400;
pulse_width=50;
example='lpi36_syllablesBandC_fs20e3_normamp0.mat';
example2='lpi36_syllablesBandC_postcut_fs20e3_normamp0.mat';

timesteps=1e3;

int_start=50;
int_width=900;
int_lambda=1/5; % firing rate
int_n=100; % number of neurons to simulate
noise_amp=.1;
noise_offset=.5;
smoothing=10;
randmod_n=6;
randmod_w=30;
int_w=[30 50];
int_breakpoints=5;
smoothing_kernel=normpdf([-30:1:30],0,8);
smoothing_kernel=smoothing_kernel./sum(smoothing_kernel);

pn_neuron=[zeros(1,pulse_t) ones(1,pulse_width) zeros(1,timesteps-pulse_width-pulse_t)];
int_flat=[zeros(1,int_start) ones(1,int_width).*int_lambda zeros(1,timesteps-int_width-int_start)];

% create the interneuron PSTH
%

s=rng;
rng(1,'twister');

int_randmod=[zeros(1,int_start) ones(1,int_width).*int_lambda zeros(1,timesteps-int_width-int_start)];
%int_breakpos=sort(randi(timesteps-int_w(2),int_breakpoints,1));
int_breakpos=round(linspace(100,timesteps-100,int_breakpoints))+randi(200)-100;

for i=1:int_breakpoints
	int_randmod(int_breakpos(i):int_breakpos(i)+randi(diff(int_w))+int_w(1))=0;
end

int_randmod(pulse_t-50:pulse_t+pulse_width+50)=int_lambda; % ensure it doesn't stop during the projection neuron

rng(s);

% stereotyped rasters

int_smoothspikes_stereotyped=cell(1,3);

for i=1:3
	int_rate=repmat(int_randmod,[int_n 1]);
	int_rate2spikes=poissrnd(int_rate);

	for j=1:size(int_rate2spikes,1)
		int_smoothspikes_stereotyped{i}(j,:)=conv(int_rate2spikes(j,:),smoothing_kernel,'same');
	end
end


int_smoothspikes_drift=cell(1,3);

for i=1:3

	rng(i+1,'twister');

	int_randmod=[zeros(1,int_start) ones(1,int_width).*int_lambda zeros(1,timesteps-int_width-int_start)];
	int_breakpoints=randi(5);
	int_breakpos=round(linspace(100,timesteps-100,int_breakpoints))+randi(200)-100;
	
	for j=1:int_breakpoints
		int_randmod(int_breakpos(j):int_breakpos(j)+randi(diff(int_w))+int_w(1))=0;
	end

	int_rate=repmat(int_randmod,[int_n 1]);
	int_rate2spikes=poissrnd(int_rate);

	for j=1:size(int_rate2spikes,1)
		int_smoothspikes_drift{i}(j,:)=conv(int_rate2spikes(j,:),smoothing_kernel,'same');
	end
end

% use an example sonogram, make 3 x 3 fig

load(fullfile(dirs.agg_dir,dirs.template_dir,example),'template');
pre_template=template;

load(fullfile(dirs.agg_dir,dirs.template_dir,example2),'template');
post_template=template;
% get length of template, align firing rates

[s,f,t]=zftftb_pretty_sonogram(pre_template.data,pre_template.fs,'filtering',300,'clipping',[-3 2],'zeropad',0);
[s2,f2,t2]=zftftb_pretty_sonogram(post_template.data,post_template.fs,'filtering',300,'clipping',[-3 2],'zeropad',0);
npoints=size(int_smoothspikes_drift{1},2);
new_fs=1/(t(end)/npoints);
new_t=[0:npoints-1]/new_fs;

fig=figure();
counter=1;
width_exp=1.15;


for i=1:3
	ax(counter)=subplot(3,3,counter);
	if i==1
		imagesc(t,f/1e3,s);
	else
		imagesc(t,f/1e3,s2);
	end
	colormap(COLORS);
	axis xy
	ylim([0 10]);
	set(gca,'xtick',[],'ytick',[]);


	counter=counter+1;
	ax(counter)=subplot(3,3,counter);
	plot(new_t,mean(int_smoothspikes_stereotyped{i}));
	ylim([0 .25]);
	axis off;

	counter=counter+1;
	ax(counter)=subplot(3,3,counter)
	plot(new_t,mean(int_smoothspikes_drift{i}),'r-');
	ylim([0 .25]);
	axis off;

	counter=counter+1;
end

for i=1:length(ax)
	pos=get(ax(i),'position');
	width_diff=width_exp*pos(3)-pos(3);
	set(ax(i),'position',[ pos(1)-width_diff pos(2) pos(3)*width_exp pos(4)]);
end

linkaxes(ax,'x');
%set(fig,'units','centimeters','position',[3 3 10 6],'paperpositionmode','auto');
%markolab_multi_fig_save(fig,fullfile(dirs.agg_dir,dirs.fig_dir),'fake_neurons','eps,png,fig,pdf','renderer','painters');

