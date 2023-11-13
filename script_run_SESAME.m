
clc

%% Load sample variables:
load('vars_for_SESAME.mat');

%% set optional input parameters

cfg.n_samples = 50; %number of samples

cfg.t_start = 160; %initial time pont of the data
cfg.t_stop = 200; %final time pont of the data

cfg.evol_exp = 0; %number of iterations (if > 0), or adaptive number of iterations (if == 0)
cfg.noise_std = 1.2e-13; %estimated noise on the data; if not specified it is self estimated by the algorithm
cfg.verbose = true;

%% run SESAME

posterior = inverse_SESAME(data, LF, sourcespace, cfg);

%% save result

TIME = clock;
save(strcat([date,'_',num2str(TIME(4)),num2str(TIME(5)),num2str(TIME(6)),'_SESAME.mat']),'posterior');


%% visualization

inverse_SESAME_viewer(posterior);
