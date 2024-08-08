
%% clear all
clear all, clc

%% Setup the environment

% Set up folder
exp_folder = 'Z:/longevity_2024/'; % path to analyse folder
baseline_data_file = 'Z:/longevity_2024/results/Adela/baseline'; % where to find baseline corrected data

% Load Data
% go to the data folder
cd(baseline_data_file);

% file name
load('tfr_miss_all.mat');
load('tfr_hits_all.mat');


%% Create grand average with ft_grandaverage 
% Within-Subjects: test difference between average power spectra 

cfg = [];
cfg.foilim = [8 30]; % specify a subset of frequencies like alpha, beta 
cfg.toilim = [-1 4]; % specify subset of latencies 
cfg.channel = 'all'; % specify subset of electrodes
cfg.parameter = 'powspctrm';


GA_hits = ft_freqgrandaverage(cfg, tfr_hits_all{:});
GA_miss = ft_freqgrandaverage(cfg, tfr_miss_all{:});

GA_sme = GA_hits;
GA_sme.powspctrm = GA_hits.powspctrm - GA_miss.powspctrm;


%% Save Data

% Save grand average for hits and miss
output_folder_avg = fullfile(exp_folder, '\results\Adela\average');

% Save grandavg_hits
grandavg_hits_filename = fullfile(output_folder_avg, 'GA_hits');
save( grandavg_hits_filename,'GA_hits');

% Save grandavg_miss
grandavg_miss_filename = fullfile(output_folder_avg, 'GA_miss');
save( grandavg_miss_filename,'GA_miss');
%% Define channel bundles

a_electrodes = {'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14', 'A15', 'A16', 'A17', 'A18', 'A19', 'A20', 'A21', 'A22', 'A23', 'A24', 'A25', 'A26', 'A27', 'A28', 'A29', 'A30', 'A31', 'A32'};
b_electrodes = {'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', ...
                       'B11', 'B12', 'B13', 'B14', 'B15', 'B16', 'B17', 'B18', 'B19', 'B20', ...
                       'B21', 'B22', 'B23', 'B24', 'B25', 'B26', 'B27', 'B28', 'B29', 'B30', ...
                       'B31', 'B32'};
                   
c_electrodes = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', ...
                      'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', ...
                      'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', ...
                      'C31', 'C32'};

d_electrodes = {'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', ...
                       'D11', 'D12', 'D13', 'D14', 'D15', 'D16', 'D17', 'D18', 'D19', 'D20', ...
                       'D21', 'D22', 'D23', 'D24', 'D25', 'D26', 'D27', 'D28', 'D29', 'D30', ...
                       'D31', 'D32'};
negcluster_electrodes  = {};

%% Script to export data for time, frequency and electrode ROIs

nsubjs=numel(tfr_hits_all);
toi=[0 3];
foi=[8 12];
eoi=a_electrodes;

dum=tfr_hits_all{1,1};

export_mat=zeros(nsubjs,2);
for n=1:nsubjs
    toi_idx=find(dum.time == toi(1)):1:find(dum.time == toi(2));
    foi_idx=find(dum.freq == foi(1)):1:find(dum.freq == foi(2));
    eoi_idx=find(ismember(dum.elec.label, eoi));
    
    export_mat(n,1)= mean(mean(mean(tfr_hits_all{1,n}.powspctrm(eoi_idx,foi_idx,toi_idx),1),2),3);
    export_mat(n,2)= mean(mean(mean(tfr_miss_all{1,n}.powspctrm(eoi_idx,foi_idx,toi_idx),1),2),3);
end

% write the matrix directly to a CSV file (without headers)
% writematrix(export_mat, 'GA_d_electrodes.csv');

%% Visualization
cfg=[];
cfg.zlim=[-0.5 0.5];
cfg.layout = ['biosemi128.lay'];
cfg.channel = 'all'; 
cfg.colormap = cm;
ft_multiplotTFR(cfg,GA_hits);
ft_multiplotTFR(cfg,GA_miss);

cfg=[];
cfg.zlim=[-0.15 0.15];
% cfg.maskparameter = [-0.05 0.05];
% cfg.maskstyle = 'outline';
cfg.layout = ['biosemi128.lay'];
cfg.channel = 'all'; 
cfg.colormap = cm;
ft_multiplotTFR(cfg, GA_sme);
cfg.channel = 'all'; 
cfg.colormap = cm;
ft_multiplotTFR(cfg, GA_sme);
