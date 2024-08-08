%% clear all
clear all, clc
%% Location of all channels
cfg = [];
cfg.layout = 'biosemi128.lay';
layout = ft_prepare_layout(cfg);

%% Time-frequency for individual data (multiplot, singleplot)
% Load Data 
cd(tfr_data_file) % go to folder 
s33 = load('subj33_tfr_miss.mat'); % subject file

% Configuration topoplot alpha/beta
cfg=[];
cfg.zlim = [-0.7 0.7];
cfg.colorbar = 'yes';
cfg.colormap = cm;
cfg.foi = [8 30]; 
cfg.toi = [0 3]; 

ft_topoplotTFR(cfg,s33.tfr_miss_baseline);
ft_singleplotTFR(cfg,s33.tfr_miss_baseline);

%% Time-frequency and Topoplot for grand average
% Load Data 
cd(output_folder_avg) % go to folder 
load('GA_hits.mat');
load('GA_miss.mat');

% Configuration singleplot_TFR
cfg=[];
cfg.zlim=[-0.5 0.5];
cfg.layout = ['biosemi128.lay'];
cfg.channel = 'all'; 
cfg.colorbar = 'yes';
cfg.foi = [8 30]; % doesn't work
cfg.toi = [0 4]; % doesn't work
cfg.colormap = cm;

ft_singleplotTFR(cfg,GA_hits);
ft_singleplotTFR(cfg,GA_miss);

% Configuration multiplot_TFR
cfg=[];
cfg.zlim=[-0.5 0.5];
cfg.xlim= [0 3];
cfg.ylim = [8 20];
cfg.layout = ['biosemi128.lay'];
cfg.channel = 'all'; 
cfg.colorbar = 'yes';
cfg.colormap = cm;

ft_topoplotTFR(cfg,GA_hits);
ft_topoplotTFR(cfg,GA_miss);


%% Time-frequency and Topoplot for differences between conditions

% Compute difference between conditions
GA_sme = GA_hits;
GA_sme.powspctrm = GA_hits.powspctrm - GA_miss.powspctrm;

cfg=[];
cfg.zlim=[-0.15 0.15];
cfg.xlim= [0 4];
cfg.ylim = [8 30];
cfg.layout = ['biosemi128.lay'];
cfg.channel = 'all'; 
cfg.colorbar = 'yes';
cfg.colormap = cm;

ft_singleplotTFR(cfg,GA_sme);
ft_topoplotTFR(cfg, GA_sme);


% cfg.maskparameter = [-0.05 0.05];
% cfg.maskstyle = 'outline';
%% Plots statistics for Hanning

% Prepare cfg for multiplotTFR
cfg = [];
cfg.parameter = 'stat';  % Field to plot
cfg.maskparameter = 'mask';  % Use the mask field to plot significant clusters
cfg.maskstyle = 'outline';  % Outline significant areas
cfg.zlim = 'maxabs';  % Color scale limits
cfg.layout = 'biosemi128.lay';  
cfg.colorbar = 'yes';  % Include colorbar
cfg.colormap = cm;

figure;
ft_multiplotTFR(cfg, stat_freq);
title('Statistical Results (Hits vs Misses)'); % not sure how to interpret this


% Prepare cfg for topoplot
cfg = [];
cfg.parameter = 'stat';
cfg.layout = 'biosemi128.lay';
cfg.colormap = cm;
cfg.colorbar = 'yes';


figure;
ft_topoplotTFR(cfg, stat_freq);
title('Hanning: Topographic Plot of Statistical Differences');
% overlay EEG channels
hold on;  % Keep the current plot
ft_plot_layout(layout, 'label', 'yes', 'box', 'no');
hold off;


% Prepare cfg for clusterplot
cfg = [];
cfg.alpha     = 0.05;
cfg.zlim=[-0.15 0.15];
cfg.alpha = 0.05;
cfg.layout    = 'biosemi128.lay';
cfg.colormap = cm;

figure;
ft_clusterplot(cfg, stat_freq);
title('Cluster Plot for Statistical Differences');

% Time-frequency over individual channels
cfg = [];
cfg.channel = 'A30';  % Replace with your channel of interest
cfg.parameter = 'stat';
figure;
ft_singleplotER(cfg, stat_freq);


%% Plot data at the channel level and overlay averaged

cfg = [];
cfg.parameter = 'powspctrm';
cfg.xlim = [-1 3];
cfg.ylim = [-0.5 0.5];
cfg.channel = 'A30';
cfg.foi = [8 30]; 

figure;
ft_singleplotER(cfg, tfr_hits_all{:});
title('Hits: Power Time-series at channel level');
hold on;
ft_singleplotER(cfg, GA_hits);

cfg = [];
cfg.xlim = [-1 3];
cfg.ylim = [-0.5 0.5];
cfg.channel = 'A30';
cfg.foi = [8 30]; 

figure;
ft_singleplotER(cfg, GA_miss); 
ft_singleplotER(cfg, GA_hits);


figure;
ft_singleplotER(cfg, GA_miss); 
ft_singleplotER(cfg, GA_hits);

