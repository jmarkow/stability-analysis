%
%
%
%
%

% for figure 1 schematic

%[options,dirs]=stan_preflight;

pulse_t=400;
pulse_width=50;
example='lpi36_syllablesBandC_fs20e3_normamp0.mat';
example2='lpi36_syllablesBandC_postcut_fs20e3_normamp0.mat';

timesteps=1e3;

int_start=50;
int_width=900;
int_lambda=1/5; % firing rate
int_n=100; % number of neurons to simulate
noise_amp=.6;
noise_offset=0;
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
int_breakpos=round(linspace(100,timesteps-100,int_breakpoints))+randi(200)-100;

for i=1:int_breakpoints
	int_randmod(int_breakpos(i):int_breakpos(i)+randi(diff(int_w))+int_w(1))=0;
end

int_randmod(pulse_t-50:pulse_t+pulse_width+50)=int_lambda; % ensure it doesn't stop during the projection neuron

rng(s);

% stereotyped rasters

int_smoothspikes_stereotyped=cell(1,3);

for i=1:50
	int_rate=repmat(int_randmod,[int_n 1]);
	int_rate2spikes=poissrnd(int_rate)+randn(size(int_rate))*noise_amp+noise_offset;

	for j=1:size(int_rate2spikes,1)
		int_smoothspikes_stereotyped{i}(j,:)=conv(int_rate2spikes(j,:),smoothing_kernel,'same');
	end
end


int_smoothspikes_drift=cell(1,3);

for i=1:3

	rng(i+1,'twister');

	int_randmod=[zeros(1,int_start) ones(1,int_width).*int_lambda zeros(1,timesteps-int_width-int_start)];
	int_breakpoints(i,:)=randi(5);
	int_breakpos=round(linspace(100,timesteps-100,int_breakpoints(i,:)))+randi(200)-100;

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

mu=cellfun(@mean,int_smoothspikes_stereotyped,'uniformoutput',false);
mu=cat(1,mu{:})';
mu_x=1:size(mu,1);

colors=colormap([ 'winter(' num2str(size(mu,2)) ')' ]);

fig=figure();

% randomly choose who lives/dies

goes=binornd(1:length(int_breakpos),.5./[1:length(int_breakpos)])>0;

circ_x=0:pi/10:2*pi;
target_goes=size(mu,2)/2+10;
target_comes=size(mu,2)/2-10;

for i=1:size(mu,2)
	% minina are still at the breakpoints cuz
	minima=int_breakpos;
	cur_x=mu_x+(i-1)*2;
	cur_x=cur_x+randn*2;
	cur_y=mu(:,i)+i*.1;

	plot(cur_x,mu(:,i)+cur_y,'k-','color',colors(i,:));
	hold on;

	for j=1:length(minima)
		if (goes(j) & i<target_goes) | (~goes(j) & i>target_comes)
			point_x=cur_x(minima(j));
			point_y=cur_y(minima(j));

			h1=patch((sin(circ_x)*4.8)+point_x,(cos(circ_x)*.01)*4.8+point_y,ones(size(circ_x)),'r','edgecolor','none')

			if goes(j)
				perc_complete=(target_goes-i)/target_goes;
				alpha(h1,min(tanh(perc_complete*1.7),1))
			else
				%perc_complete=((size(mu,2)-i)/size(mu,2))*2;
				dist=size(mu,2)-target_comes;
				perc_complete=(i-target_comes)/dist
				alpha(h1,min(tanh(perc_complete)*1.7,1))
			end
		end
	end
end
