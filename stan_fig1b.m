function stan_fig1b
%
%
%

load custom_colormaps;
fake_fig=stan_fake_neurons(fee_map);
tightfig(fake_fig);
set(fake_fig,'units','centimeters','position',[4 4 5.5 3.5],'paperpositionmode','auto');

