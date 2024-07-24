%% clear all
clear all, clc

%% Setup the environment

% Set up folder
exp_folder = 'Z:/longevity_2024/'; % path to analyse folder
preprocessed_eeg_file = 'Z:/longevity_2024/data/EEG_data'; % where to find preprocesed data

% edit static list
subjects = {
    'subj01', 'subj02', 'subj04', 'subj05', 'subj08', 'subj09', 'subj10', ...
    'subj12', 'subj13', 'subj14', 'subj15', 'subj16', 'subj17', 'subj19', 'subj20', ...
    'subj21', 'subj22', 'subj23', 'subj24', 'subj25', 'subj26', 'subj27', 'subj28', 'subj29', 'subj30', ...
    'subj32', 'subj33', 'subj34', 'subj35', 'subj36', 'subj37', 'subj38', 'subj40', ...
    'subj41', 'subj43', 'subj46', 'subj47', 'subj49',  ...
     'subj52', 'subj53', 'subj54', 'subj56', 'subj57', 'subj58', 'subj59', 'subj60', ...
    'subj61', 'subj62', 'subj63'
};


% initialize Fieldtrip
addpath(fullfile(exp_folder, '\analyses\plug-ins\fieldtrip-20231025'))
addpath(fullfile(exp_folder, '\analyses\functions'))
ft_defaults

%% Loop

for i = 1:length(subjects)

    subj_id = subjects{i};
    task = 'ratings_preprocessed';

% Load Data
% go to the data folder
cd(preprocessed_eeg_file)

% file name
filename = [subj_id '_' task '.mat'];
load(filename)

cfg = [];
cfg.dataset = filename;

%% Preprocess the data for hits/ misses
cfg.trials =  find(data_clean.trialinfo(:, 7) == 1); % Hit column
remembered_trials  = ft_preprocessing(cfg, data_clean);

cfg.trials =  find(data_clean.trialinfo(:, 7) == 0); % Miss column
forgotten_trials  = ft_preprocessing(cfg, data_clean);

%% Time frequency analysis with a Hanning taper and fixed window length
cfg            = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmconvol';
cfg.taper      = 'hanning';
cfg.toi        = -1 : 0.05 : 4;
cfg.foi        = 8:30;
cfg.t_ftimwin  = 6./cfg.foi;
cfg.keeptrials = 'no';  % no need to save trial info
tfr_miss      = ft_freqanalysis(cfg, forgotten_trials);

tfr_hits      = ft_freqanalysis(cfg, remembered_trials);

%% Do ft_baseline 

cfg              = [];
cfg.baseline     = [-0.5 -0.1]; % 
cfg.baselinetype = 'relchange';  % we use a relative baseline
cfg.xlim         = [-1 4]; % time window we want to look at
cfg.ylim         = [8 30];  % plotting alpha/beta bands
cfg.zlim         = 'maxabs';
cfg.marker       = 'on';
cfg.colorbar     = 'yes';
cfg.layout       = 'biosemi128.lay';

tfr_miss_baseline = ft_freqbaseline(cfg, tfr_miss);
tfr_hits_baseline = ft_freqbaseline(cfg, tfr_hits);

% 
%% Save individual data for hits and miss

output_folder = fullfile(exp_folder, '\results\Adela\tfr');

% set filenames and save data
hits_filename = fullfile(output_folder, [subj_id '_tfr_hits.mat']);
save(hits_filename, 'tfr_hits_baseline');

miss_filename = fullfile(output_folder, [subj_id '_tfr_miss.mat']);
save(miss_filename, 'tfr_miss_baseline');

% Store the baseline corrected data in the cell arrays
tfr_miss_all{i} = tfr_miss_baseline;
tfr_hits_all{i} = tfr_hits_baseline;

end
%% Save the baseline corrected data for all participants 
output_folder_baseline = fullfile(exp_folder, '\results\Adela\baseline');

% Save tfr_miss_all
miss_all_filename = fullfile(output_folder_baseline, 'tfr_miss_all');
save( miss_all_filename,'tfr_miss_all');

% Save tfr_hits_all
hits_all_filename = fullfile(output_folder_baseline, 'tfr_hits_all');
save( hits_all_filename,'tfr_hits_all');

%% END OF SCRIPT 
% the other stuff below is failed attempts. only keeping it for inspiration

%% Multi plot TFR
% Define channel bundles

channels_bundle_A = {'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12', 'A13', 'A14', 'A15', 'A16', 'A17', 'A18', 'A19', 'A20', 'A21', 'A22', 'A23', 'A24', 'A25', 'A26', 'A27', 'A28', 'A29', 'A30', 'A31', 'A32'};
b_bundle_electrodes = {'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', ...
                       'B11', 'B12', 'B13', 'B14', 'B15', 'B16', 'B17', 'B18', 'B19', 'B20', ...
                       'B21', 'B22', 'B23', 'B24', 'B25', 'B26', 'B27', 'B28', 'B29', 'B30', ...
                       'B31', 'B32'};
                   
c_bundle_electrodes = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', ...
                      'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', ...
                      'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', ...
                      'C31', 'C32'};

d_bundle_electrodes = {'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', ...
                       'D11', 'D12', 'D13', 'D14', 'D15', 'D16', 'D17', 'D18', 'D19', 'D20', ...
                       'D21', 'D22', 'D23', 'D24', 'D25', 'D26', 'D27', 'D28', 'D29', 'D30', ...
                       'D31', 'D32'};
                   
