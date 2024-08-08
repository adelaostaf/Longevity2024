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
cfg.frequency        = [8 30]; % Frequency of interest 

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
%% Identify significant clusters

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


%% Extract unique channels, frequencies, and timepoints for the current cluster
clall = find(pow_negcluster == 1);
[cx_chan, cy_freq, cz_time] = ind2sub(size(pow_negcluster), clall);

% Initialize cell arrays to store results for each cluster
cluster_channels = {};
cluster_frequencies = {};
cluster_timepoints = {};

% Create a structure to store the detailed information
cluster_info = struct('ClusterID', [], 'Timepoint', [], 'Channels', [], 'Frequencies', []);

% Loop through each significant cluster
unique_clusters = unique(pow_negcluster(clall));
for i = 1:length(unique_clusters)
    cluster_id = unique_clusters(i);
    
    % Find all indices belonging to the current cluster
    current_cluster_indices = find(pow_negcluster == cluster_id);
    [current_cx_chan, current_cy_freq, current_cz_time] = ind2sub(size(pow_negcluster), current_cluster_indices);
    
    % Extract unique channels, frequencies, and timepoints for the current cluster
    chans_idx = unique(current_cx_chan);
    clust_chans = stat.label(chans_idx);
    
    chans_idy = unique(current_cy_freq);
    clust_freq = stat.freq(chans_idy);
    
    chans_idz = unique(current_cz_time);
    clust_time = stat.time(chans_idz);
    
    % Store the results in cell arrays
    cluster_channels{i} = clust_chans;
    cluster_frequencies{i} = clust_freq;
    cluster_timepoints{i} = clust_time;
    
    % Display the overview for the current cluster
    fprintf('Cluster %d:\n', cluster_id);
    
    % For each timepoint in the cluster, display the contributing channels and frequencies
    unique_times_in_cluster = unique(current_cz_time);
    for t = 1:length(unique_times_in_cluster)
        time_idx = unique_times_in_cluster(t);
        relevant_indices = find(current_cz_time == time_idx);
        relevant_chans = unique(current_cx_chan(relevant_indices));
        relevant_freqs = unique(current_cy_freq(relevant_indices));
        
        fprintf('Timepoint: %.3f\n', stat.time(time_idx));
        fprintf('Channels contributing:\n');
        disp(stat.label(relevant_chans));
        fprintf('Frequencies contributing:\n');
        disp(stat.freq(relevant_freqs));
        
        % Add information to the structure
        new_entry = struct('ClusterID', cluster_id, 'Timepoint', stat.time(time_idx), ...
                           'Channels', {stat.label(relevant_chans)}, 'Frequencies', {stat.freq(relevant_freqs)});
        cluster_info = [cluster_info; new_entry];
    end
    
    fprintf('-------------------------\n');
end

% Convert the structure to a table
cluster_table = struct2table(cluster_info);

% Save the structure as a MAT file
clustter_filename = fullfile(exp_folder, 'results', 'Adela', 'trialinfo.mat');
save(clustter_filename, 'cluster_info');

%% To compute mean, sd, 95% CI for negative cluster 

% Initialize arrays to store statistics for significant clusters
mean_hits = [];
mean_misses = [];
std_hits = [];
std_misses = [];
sem_hits = [];
sem_misses = [];
ci_diff = [];

% Extract data points from the significant negative cluster
if any(pow_negcluster(:))
    % Identify all significant clusters
    clall = find(pow_negcluster == 1);
    [cx_chan, cy_freq, cz_time] = ind2sub(size(pow_negcluster), clall);
    
    % Extract the power values from the significant cluster for hits and misses
    sig_cluster_data_hits = [];
    sig_cluster_data_misses = [];
    for i = 1:length(clall)
        sig_cluster_data_hits = [sig_cluster_data_hits; tfr_hits_all{1}.powspctrm(cx_chan(i), cy_freq(i), cz_time(i))];
        sig_cluster_data_misses = [sig_cluster_data_misses; tfr_miss_all{1}.powspctrm(cx_chan(i), cy_freq(i), cz_time(i))];
    end
    
    % Compute mean, standard deviation, and SEM for hits and misses
    mean_hits = mean(sig_cluster_data_hits);
    mean_misses = mean(sig_cluster_data_misses);
    
    std_hits = std(sig_cluster_data_hits);
    std_misses = std(sig_cluster_data_misses);
    
    sem_hits = std_hits / sqrt(length(sig_cluster_data_hits));
    sem_misses = std_misses / sqrt(length(sig_cluster_data_misses));
    
    % Compute the difference in means
    mean_diff = mean_hits - mean_misses;
    
    % Compute the standard error of the difference
    sem_diff = sqrt(sem_hits^2 + sem_misses^2);
    
    % Compute the 95% CI for the difference in means
    alpha = 0.05;
    t_crit = tinv(1 - alpha/2, length(sig_cluster_data_hits) + length(sig_cluster_data_misses) - 2);
    ci_diff = [mean_diff - t_crit * sem_diff, mean_diff + t_crit * sem_diff];
end

% Display the results
disp('Negative Significant Cluster Statistics:');
disp(['Mean (Hits): ', num2str(mean_hits)]);
disp(['Mean (Misses): ', num2str(mean_misses)]);
disp(['Mean Difference: ', num2str(mean_diff)]);
disp(['Standard Deviation (Hits): ', num2str(std_hits)]);
disp(['Standard Deviation (Misses): ', num2str(std_misses)]);
disp(['Standard Error of the Mean (Hits): ', num2str(sem_hits)]);
disp(['Standard Error of the Mean (Misses): ', num2str(sem_misses)]);
disp(['Standard Error of the Difference: ', num2str(sem_diff)]);
disp(['95% CI for Difference: [', num2str(ci_diff(1)), ', ', num2str(ci_diff(2)), ']']);


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

    subplot(1,1,it)
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
    

    title([num2str(round(plotcluster.time,2)) ' sec '])

    if it == 1
        it = 1;
        figure
    else
        it = it+1;
    end

end 

%% Alternative to plotting to save each figure individually
% Initialize the figure counter
it = 1;

% Directory to save the figures
save_dir = 'Z:/longevity_2024/results/Adela/cluster_tpoplots'; % Specify the directory where you want to save the figures

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

    % Create a new figure for each topography
    figure;
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

    title([num2str(round(plotcluster.time, 2)) ' s ']);

    % Define the filename and save the figure as a PNG
    filename = fullfile(save_dir, sprintf('topography_%02d.png', it));
    saveas(gcf, filename);

    % Close the figure after saving to free up memory
    close(gcf);

    it = it + 1;
end

%% Calculate Cohen's d

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

%% Identify clustersize

if isfield(stat, 'negclusterslabelmat')
    % Extract the label matrix for negative clusters
    neg_clusters_labelmat = stat.negclusterslabelmat;

    % Find the unique cluster labels (excluding 0, which indicates no cluster)
    unique_neg_labels = unique(neg_clusters_labelmat);
    unique_neg_labels(unique_neg_labels == 0) = []; % Remove the zero label

    % Initialize a structure to store cluster sizes
    neg_cluster_sizes = struct();

    % Iterate through each unique cluster label
    for i = 1:length(unique_neg_labels)
        label = unique_neg_labels(i);
        % Calculate the size of the cluster by counting occurrences of the label
        cluster_size = sum(neg_clusters_labelmat(:) == label);
        
        % Store the cluster size in the structure
        neg_cluster_sizes.(['Cluster_' num2str(label)]) = cluster_size;

        % Display the cluster size
        fprintf('Negative Cluster %d: Size = %d\n', label, cluster_size);
    end
else
    disp('No negative clusters found.');
end


%% Save structures

% Specify the folder path where you want to save the file
results_folder = fullfile(exp_folder, '\results\Adela'); 

% Combine the folder path with the file name
filePath = fullfile(results_folder, 'stat.mat');

% Save the structure to the specified folder
save(filePath, 'stat');
save(filePath, 'effect_rectangle');

