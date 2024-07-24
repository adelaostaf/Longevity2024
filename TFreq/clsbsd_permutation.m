%% Cluster-based permutation statistics
% 
% This script computes a dependent-samples cluster-based permutation test

clear all, close all

%% Setup the environment

% Set up folder
exp_folder = 'Z:/longevity_2024/'; % path to analyse folder
bsl_data_file = 'Z:/longevity_2024/results/Adela/baseline'; 
average_data_file = 'Z:/longevity_2024/results/Adela/average'; % where to find grand averaged data

% Load Data
% go to the data folder
cd(bsl_data_file);

% load pwrspctm
load('tfr_hits_all.mat');
load('tfr_miss_all.mat');

% load neighbours
load(fullfile(exp_folder, 'analyses\plug-ins\biosemi\neighbours.mat')) % load as neighbours

%% Dependent-samples cluster-based permutation test
% hits = tfr_hits_all
% misses = tfr_miss_all

rng(457); % You can use any integer as the seed value

cfg = [];
cfg.channel          = 'all';
cfg.latency          = 'all';
cfg.frequency        = [8 20]; % Frequency of interest 

cfg.statistic        = 'ft_statfun_depsamplesT';
cfg.alpha            = .05;
cfg.tail             = -1;  % one-sided test
cfg.correcttail      = 'alpha'; %
cfg.ivar             = 1;
cfg.uvar             = 2;

cfg.method           = 'montecarlo';
cfg.design           = [ones(1, length(tfr_hits_all)) ones(1, length(tfr_miss_all)) * 2; 1:length(tfr_hits_all) 1:length(tfr_miss_all)];

cfg.neighbours = neighbours; 

cfg.correctm         = 'cluster';
cfg.numrandomization = 2000; % 
cfg.clusterthreshold = 'nonparametric_common';
cfg.clusteralpha     = 0.05;

cfg.clusterstatistic = 'wcm';
cfg.minnbchan        = 3; % justify why this was the case, nr neighbours that exceeds threshold
cfg.clustertail      = -1;  % one-sided test

[stat] = ft_freqstatistics(cfg, tfr_hits_all{:},tfr_miss_all{:});

% cd('H:\data\')
% save cluster stat

pow_poscluster = [];
pow_negcluster = [];

% identify clusters with a pvalue <.05
if isfield(stat, 'negclusterslabelmat')
    neg_cluster_pvals = [stat.negclusters(:).prob]; % clusters pvalues
    neg_clust = find(neg_cluster_pvals < 0.05); % are there any clusters with p<.05?
    pow_negcluster       = ismember(stat.negclusterslabelmat, 1);
end
if isfield(stat, 'posclusterslabelmat')
    pos_cluster_pvals = [stat.posclusters(:).prob]; % clusters pvalues
    pos_clust = find(pos_cluster_pvals < 0.05); % are there any clusters with p<.05?
    pow_poscluster       = ismember(stat.posclusterslabelmat, 1);
end
% 1 neg sig cluster found

%% to know cluster extension
% to do for each significant cluster

clear c* onset_cl end_cl e* on*
clall= find(pow_negcluster==1);
[cx_chan, cy_freq, cz_time] = ind2sub(size(pow_negcluster), clall);

% to know which channels belong to the cluster
chans_idx = intersect(1:length(stat.label),cx_chan);
clust_chans = stat.label(chans_idx);

% to know which frequencies belong to the cluster
chans_idy = intersect(1:length(stat.freq),cy_freq);
clust_freq = stat.freq(chans_idy);

% to know which timepoints belong to the cluster
chans_idz = intersect(1:length(stat.time),cz_time);
clust_time = stat.time(chans_idz);
%% Extract t-values for significant cluster points
t_values_cluster = zeros(length(clall), 4); % Initialize a matrix to hold t-values and their corresponding indices
for i = 1:length(clall)
    t_values_cluster(i, 1) = cx_chan(i); % Channel index
    t_values_cluster(i, 2) = cy_freq(i); % Frequency index
    t_values_cluster(i, 3) = cz_time(i); % Time index
    t_values_cluster(i, 4) = stat.stat(cx_chan(i), cy_freq(i), cz_time(i)); % T-value
end

% Optionally, create a table for better readability
t_values_table = array2table(t_values_cluster, 'VariableNames', {'ChannelIdx', 'FrequencyIdx', 'TimeIdx', 'TValue'});

%% Save t_values_table

tvalinfo_folder = fullfile(exp_folder, '\results\Adela\trialinfo');
csv_filename = fullfile(tvalinfo_folder, 'tvalues_info.csv'); %t-values for critical cluster


% Save tfr_miss_all
tvalal_filename = fullfile(tvalinfo_folder, 'tval_info');
save( tvalal_filename,'t_values_table');

% Save the counts table as a CSV file
writetable(t_values_table, 't_values_table.csv');


%% plotting
% subtraction
% load GA
cd(output_folder_avg);

% file name
load('GA_hits.mat');
load('GA_miss.mat');

cfg =[];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
subtraction    = ft_math(cfg, GA_hits, GA_miss);

% ugly adjustement 
subtraction2 = subtraction;
subtraction2.time = [];
subtraction2.powspctrm  = []; 
ct = 1;
for tp = 1:length(subtraction.time)
    if stat.time(ct) == subtraction.time(tp)
        subtraction2.time(ct) = subtraction.time(tp);

        subtraction2.powspctrm(:,:,ct) = subtraction.powspctrm(:,:,tp);
        ct = ct + 1;
    end
end

% found cluster
data_cluster = double(pow_negcluster);
%clusters = {'poscluster', 'negcluster'};

% find time points included in the cluster
tpcluster = find(squeeze(sum(sum(data_cluster,1),2)))';
it = 1; 
figure
for tp = tpcluster % Create a topography for each time point within the cluster
    temp_cluster = data_cluster(:,:,tp);

    % Find channels and frequencies included in the cluster in each time point
    [chcluster, freqcluster] = find(temp_cluster);
    chcluster = unique(chcluster);
    freqcluster = unique(freqcluster);

    plotcluster = subtraction2;
    plotcluster.powspctrm = plotcluster.powspctrm(:,freqcluster,tp);
    plotcluster.time = plotcluster.time(tp);
    plotcluster.freq = plotcluster.freq(freqcluster);

    subplot(4,4,it)
    cfg                  = [];
    cfg.layout           = 'biosemi128.lay';
    cfg.parameter        = 'powspctrm';
    cfg.marker           = 'off';
    cfg.highlight        = 'on';
    cfg.highlightsymbol  = '.';
    cfg.highlightchannel = chcluster;
    cfg.highlightcolor   = [0 0 0];
    cfg.figure           = 'gca';
    cfg.highlightsize    = 7;
    cfg.colormap         = cm; 
    cfg.zlim             =[-0.1 0.1];
    ft_topoplotTFR(cfg,plotcluster)
    

    title([num2str(round(plotcluster.time,2)) ' s; '])

    if it == 16
        it = 1;
        figure
    else
        it = it+1;
    end

end 

%% Plot probability matrix: Error using image

% Color data must be an m-by-n-by-3 or m-by-n matrix.
figure;
imagesc ( stat.time,0:400, -log10(stat.prob))
colorbar
%% Plot spatial distribution of EEG channels part of the largest cluster: Does not show
cfg = [];
cfg.style     = 'blank';
cfg.layout    = 'biosemi128.lay';
cfg.highlight = 'on';
cfg.highlightchannel = find(any(stat.mask,2));
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_miss)
title('Nonparametric: significant with cluster-based multiple comparison correction')
%% Experiment Mapping effect: Doesnt work

cfg = [];
cfg.channel = clust_chans;
cfg.latency = [min(clust_time) max(clust_time)];
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'yes';
cfg.avgoverfreq = 'yes';
cfg.method = 'analytic';
cfg.statistic = 'cohensd';
cfg.ivar             = 1;
cfg.uvar             = 2;
cfg.design           = [ones(1, length(tfr_hits_all)) ones(1, length(tfr_miss_all)) * 2; 1:length(tfr_hits_all) 1:length(tfr_miss_all)];

effect_rectangle = ft_freqstatistics(cfg, tfr_hits_all{:},tfr_miss_all{:});
disp(effect_rectangle)

%% Save structures

% Specify the folder path where you want to save the file
results_folder = fullfile(exp_folder, '\results\Adela'); 

% Combine the folder path with the file name
filePath = fullfile(results_folder, 'stat.mat');

% Save the structure to the specified folder
save(filePath, 'stat');
save(filePath, 'effect_rectangle');

