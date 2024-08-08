%% clear all
clear all, clc

%% Setup the environment

% Set up folder
exp_folder = 'Z:/longevity_2024/'; % path to analyse folder
preprocessed_eeg_file = fullfile(exp_folder, '/data/EEG_data'); % individual preprocesed data folder
baseline_data_file = fullfile(exp_folder,'/results/Adela/baseline'); % baseline corrected data ( tfr & wave) folder
output_folder_avg = fullfile(exp_folder, '/results/Adela/average'); % grand average for hits and miss folder
tfr_data_file = fullfile(exp_folder,'/results/Adela/tfr') % individual tfr data folder
stat_file = fullfile(exp_folder,'/results/Adela/statistics'); % statistics for tfr and wave

% Colormap
addpath(fullfile(exp_folder,'/analyses/scripts/Adela/brewermap')) % brewermap folder
cm = brewermap([],'-Blues');

%% initialize Fieldtrip
addpath(fullfile(exp_folder, '\analyses\plug-ins\fieldtrip-20231025'))
addpath(fullfile(exp_folder, '\analyses\functions'))
ft_defaults
