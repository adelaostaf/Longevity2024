%% Preprocessing RETRIEVAL block - longevity EEG 
%
clear; clc

%% Set paths
% Experiment folder path
exp_folder = 'Z:\longevity_2024\'; % path to analyse folder
eeg_folder = 'Y:\'; % where to find the raw data
csv_folder = fullfile(exp_folder,'\data\beh_data'); % where to find the behavioural data folder

% set info
subj_id = 'subj36'; % change with the ID of the subj to analyse
task = 'test';

% initialize Fieldtrip
addpath(fullfile(exp_folder, '\analyses\plug-ins\fieldtrip-20231025'))
addpath(fullfile(exp_folder, '\analyses\functions'))
ft_defaults

%% 1.Step: Import data and filter continuous data with the following:
% HPF: 0.5
% LPF: 80
% BSF: [49 51]
%
% go to the data folder
cd(eeg_folder)

% file name
filename = [subj_id '_' task '_raw.bdf'];

% apply filters
cfg           = [];
cfg.dataset   = filename;
cfg.channel = 'EEG';
cfg.bpfilter  = 'yes'; % enables bandpass filter
cfg.bpfreq = [0.5 80];% [high-pass low-pass] in Hz
cfg.bpfiltord = 3; % filter order
cfg.bsfilter = 'yes'; % enables bandstop filter to remove 50Hz noise
cfg.bsfreq = [49 51]; % Hz
data_filtered = ft_preprocessing(cfg);

% check after filtering
cfg = [];
cfg.channel = 'all';
cfg.method = 'mtmfft';
cfg.output = 'pow';
cfg.taper = 'boxcar';
cfg.foi= 0.5:0.5:100;
cfg.tapsmofrq = 2;
freq_filter = ft_freqanalysis(cfg, data_filtered);

figure;
hold on;
plot(freq_filter.freq, freq_filter.powspctrm)

%% 2. Step: Cut data into trials with -2000:1:5000 ms
%
%
% 
cfg = [];
cfg.dataset = filename;
cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype  = 'STATUS';
cfg.trialdef.eventvalue = 63492; % stimulus onset
cfg.trialdef.prestim    = 2; % pre-stimulus window in seconds
cfg.trialdef.poststim   = 5;% post-stimulus window in seconds
cfg = ft_definetrial(cfg);

% extract trials info
trl = cfg.trl;

% retrieve trial info from csv
cfg = [];
cfg.csv_folder = csv_folder;
cfg.task = 'TestData';
cfg.subj_id = subj_id;
cfg.group = 'short_group';
combined_trl_info = trialfun_retrieval_longevity_TIME(cfg);
trl = [trl, combined_trl_info]; % add the condition (0,180) to the trl structure

% create trials
cfg = [];
cfg.trl = trl;
trials = ft_redefinetrial(cfg, data_filtered);

%% 3. Step: Reject bad channels and very bad trials using ft_rejectvisual
%
%
% remove extreme trials/channels
cfg          = [];
cfg.method   = 'summary';
artrej_1       = ft_rejectvisual(cfg, trials);

% trial by trial - remove noisy trials and bad channels
cfg          = [];
cfg.method   = 'trial';
cfg.channel = 'EEG';
artrej_2 = ft_rejectvisual(cfg, artrej_1);

% get the removed channels
bad_channels = setxor(artrej_2.label(:,1),trials.label(:,1));
bad_chan_idx = find(~ismember(trials.label(:,1),artrej_2.label(:,1)));
good_channels = setxor(trials.label(:,1), bad_channels);
good_chan_idx = find(~ismember(trials.label(:,1), bad_channels));

% remove bad muscle
cfg = [];
cfg.artfctdef.zvalue.channel = good_channels;
cfg.artfctdef.zvalue.cutoff = 4;
cfg.artfctdef.zvalue.trlpadding = 0;
cfg.artfctdef.zvalue.fltpadding = 0;
cfg.artfctdef.zvalue.artpadding = 0.1;
cfg.artfctdef.zvalue.bpfilter = 'yes';
cfg.artfctdef.zvalue.bpfreq = [110 140];
cfg.artfctdef.zvalue.bpfiltord = 3;
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.hilbert  = 'yes';
cfg.artfctdef.zvalue.boxcar = 0.2;
cfg.artfctdef.zvalue.interactive = 'yes';
[~, artifact_muscle] = ft_artifact_zvalue(cfg, artrej_2);

% remove the identified trials with muscle artifacts
cfg = [];
cfg.artfctdef.reject = 'complete';
cfg.artfctdef.muscle.artifact = artifact_muscle;
artrej_2_ma = ft_rejectartifact(cfg, artrej_2);

% visualise the data
cfg = [];
cfg.channel = good_channels;
cfg.layout = 'biosemi128.lay';
cfg.viewmode = 'vertical';
cfg.allowoverlap = 'yes';
artif = ft_databrowser(cfg,artrej_2_ma);

%% 4. ICA: runica; reject eye blinks and horizontal eye movements, and heart beats (if you see them) and backproject data;
%
%
% run the ICA on downsampled data
cfg = [];
cfg.resamplefs = 300;
cfg.detrend = 'no';
data_res = ft_resampledata(cfg, artrej_2_ma);

% perform ICA
cfg = [];
cfg.method = 'runica';
cfg.channel = good_channels; 
comp = ft_componentanalysis(cfg, data_res);

% plot ICA components for visual inspection
figure;
cfg = [];
cfg.component = 1:20; % components to plot
cfg.comment = 'no';
cfg.marker = 'no';
cfg.layout = 'biosemi128.lay';
ft_topoplotIC(cfg, comp)

% look at the time course of the components
cfg =[];
%cfg.channel = [1:2 18]; % components to plot
cfg.layout = 'biosemi128.lay';
cfg.viewmode = 'component';
ft_databrowser(cfg, comp);

% remove the components and project back
cfg = [];
cfg.component = [1 12 19]; % components to remove
% apply it on the non-downsampled data:
data_ica = ft_rejectcomponent(cfg, comp, artrej_2_ma); 

% have a look 
cfg = [];
cfg.channel = 'all';
cfg.layout = 'biosemi128.lay';
cfg.allowoverlap = 'yes';
cfg.viewmode = 'vertical';
after_ica = ft_databrowser(cfg,data_ica);

%% 5. Interpolate missing channels and re-reference to average reference
%
%
%
% re-add bad channels to the trial struct and apply a weighted interpolation
data_to_interpolate = prepare_channels_to_interpolate(data_ica, trials.label, bad_channels, good_chan_idx, bad_chan_idx);

% load electrode configuration
load(fullfile(exp_folder, 'analyses/plug-ins/biosemi/electrodes.mat'))

% define neighbours channels
cfg_neighb = [];
cfg_neighb.method = 'triangulation';
cfg_neighb.feedback ='yes';
cfg_neighb.elec = elec; % elec structure extracted from polhemus
neighbours = ft_prepare_neighbours(cfg_neighb, data_to_interpolate);

% interpolate bad channels
cfg = [];
cfg.badchannel = bad_channels;
cfg.elec = elec;
cfg.neighbours = neighbours;
data_fixed = ft_channelrepair(cfg,data_to_interpolate);

% re-reference to average
cfg =[];
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'avg';
data_reref = ft_preprocessing(cfg, data_fixed);

%% 6. reject any remaining artefacts (muscle or other stuff)
%
%
%
cfg = [];
cfg.channel = 'all';
cfg.layout = 'biosemi128.lay';
cfg.allowoverlap = 'yes';
cfg.viewmode = 'vertical';
artrej_3 = ft_databrowser(cfg, data_reref);

visual_artifacts = artrej_3.artfctdef.visual.artifact;
% remove artifacts
cfg = [];
cfg.artfctdef.reject = 'complete';
cfg.artfctdef.visual.artifact = visual_artifacts;
data_clean = ft_rejectartifact(cfg, data_reref);

% save pre-processed data
output_folder = fullfile(exp_folder,['data\EEG_data\' subj_id '_' task '_preprocessed']);
save(output_folder,"data_clean");

%% preprocessing completed 