%%
% [ripple_stage_04a_dimensionality] find ripples in hippocampal channels,
% extract and realign data based on ripple events.
% Do PCA to estimate dimensionality of correct and incorrect.
%                  Bernhard Staresina [bernhard.staresina@psy.ox.ac.uk]
%                  Casper Kerren      [kerren@cbs.mpg.de]

clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

% [~,ftpath]=ft_version;

%% path settings
paths = config_paths();

settings = [];
settings.base_path_castle = paths.base_path;

settings.data_dir           = paths.data_dir;
settings.save_dir           = paths.save_dir;
settings.data_dir_channels  = paths.channels_dir;
settings.anatomy_dir        = paths.anatomy_dir;
settings.AAL_dir            = paths.AAL_dir;
settings.SPM_dir            = paths.SPM_dir;


addpath('/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/scripts_to_publish')

load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;

settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
subjects                    = {'CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK'};
SubjectIDs                  = {'01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK'};

settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};
settings.healthyhemi        = {'R' 'LR' 'R' 'L' 'L' 'R' 'R' 'R' 'R' 'R' 'R' 'LR'};

addpath(genpath([paths.MVPA_Light_master]))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([paths.subfunctions])
addpath(genpath([paths.help_functions]))
addpath(genpath([paths.plotting]))


%% pre-decoding settings
% ripple extraction

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1; % remove co-occuring ripples
settings.time_to_excl_RT            = .25; % exclude last 250 ms of trials, to make sure ripple event was in trial
settings.solo_ripple                = 1; % pick one (maxEnv) ripple per trial if multiple ripple events are found
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur
settings.do_surrogates              = 0; % switch time of ripples between trials, 1 == for all trials, 2 == for correct trials only

%% decoding settings

% data preprocessing
settings.ori_fsample            = 1000; % original sample frequency
settings.do_resample            = 100; % [] or sample frequency

% baseline and zscoring settings
settings.zscore_data4class      = 1;
settings.bs_correct             = 1;
settings.bs_period              = [-.2 0]; % [-.5 -.1]
settings.bs_trim                = 0; % can be 0. amount of % to trim away when calculating the baseling


% smoothing options
settings.do_smoothdata          = 1; % use matlabs smoothdata function for running average
settings.smooth_win             = .200; % .100, [] running average time window in seconds;

% time of interest
settings.TOI_train              = [-.5 6]; % to be able to have -1 to 1 from start of retrieval to end (0-5seconds)
settings.timesteps_train        = ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample); % in s. if you want it to take less sample points multiply [e.g., ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample)*2

settings.TOI_test               = [-1.2 1.2]; % time around ripple (I take this time window to get a proper estimate around the edges too. Only look at -1 to 1 later.
settings.timesteps_test         = ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample); % in s

settings.classifier             = 'lda';
settings.metric                 = 'auc';

%% channel settings

settings.channel_sel            = 1; % 1 exclude hippo, 2 only hippo, 3 all channels

%% settings pca

settings.smooth_before_dim      = 1; % smooth Nans before doing dim reduction
settings.nu_time_points         = 60; % time of sliding window in ms
settings.prc_overlap            = .9; % percentage overlap sliding window
settings.decode_components      = 0; % Decode the PCA components I picked.
%% start for loop

timeaxis = settings.TOI_test(1):settings.timesteps_test:settings.TOI_test(2);
freqaxis = settings.TOI_train(1):settings.timesteps_train:settings.TOI_train(2);

perf    = cell(1,numel(subjects));
channs  = cell(1,numel(subjects));

RT_all_subj_correct     = cell(1,numel(subjects));
RT_all_subj_incorrect   = cell(1,numel(subjects));
ripple_time_correct     = cell(1,numel(subjects));
ripple_time_incorrect   = cell(1,numel(subjects));

numWorkers = 8; %

parpool('local', numWorkers);

tic
parfor isubject = 1:numel(subjects)


    fprintf('processing subject %01d/%02d\n',isubject,numel(subjects));

    %% LOAD data



    tmp             = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',subjects{isubject}]);
    data            = tmp.data;
    onsets_session  = tmp.onsets_session;
    tmp             = [];



    %% load channels and restrict to channels that are in hippocampus

    SubjectID               = SubjectIDs{isubject};

    tmp                     = load([settings.base_path_castle,'ripple_project_publication_for_replication/templates/channels_hipp_ripples']);
    channels_hipp_ripples   = tmp.channels_hipp_ripples;
    tmp                     = [];
    channels                = channels_hipp_ripples(isubject,:);
    channels                = channels(~cellfun(@isempty,channels));

    cfg         = [];
    cfg.channel = channels;
    data        = ft_selectdata(cfg, data);

    % If NaNs in recording, interpolate these.

    cfg             = [];
    cfg.prewindow   = 3;
    cfg.postwindow  = 3;
    data            = ft_interpolatenan(cfg,data);


    %% Load data to realign based on cue onset encoding and based on ripples

    data_in                 = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',subjects{isubject}],'data','onsets_session');
    data                    = data_in.data;
    data_in                 = [];
    onsets                  = onsets_session - data.time{1}(1)*data.fsample;
    data.sampleinfo         = 1+data.sampleinfo - data.sampleinfo(1);
    data.time{1}            = 1/data.fsample+data.time{1}-data.time{1}(1);

    %% channel selection for data

    % 1 exclude hippo, 2 only hippo, 3 all channels
    tmp                             = load([settings.data_dir_channels,'/channels_to_exclude_all_hipp_both_hem.mat']);
    channels_to_exclude_all_hipp    = tmp.channels_to_exclude_all_hipp;
    tmp                             = load([settings.data_dir_channels,'/channels_hipp_ripples.mat']);
    channels_hipp_ripples           = tmp.channels_hipp_ripples;

    cfg         = [];
    cfg.channel = data.label;

    switch settings.channel_sel
        case 1
            cfg.channel = setdiff(setdiff([data.label],char(channels_to_exclude_all_hipp{isubject,:})),settings.scalp_channels);
        case 2
            cfg.channel = intersect(cellstr(setdiff([data.label],settings.scalp_channels)),char(channels_hipp_ripples{isubject,:}));
        case 3
            cfg.channel = setdiff(data.label,settings.scalp_channels);
    end

    %--- load channel info and coordinates
    [~,~,entries]   = xlsread(fullfile(settings.anatomy_dir,'well01_ripples_ROIs_w_labels.xlsx'));

    subject_colums  = entries(1,:);
    column_names    = entries(2,:);
    these_labels    = entries(3:end,strcmp(subject_colums,['s' SubjectIDs{isubject}]) & strcmp(column_names,'label'));
    these_x         = cell2mat(entries(3:end,strcmp(subject_colums,['s' SubjectIDs{isubject}]) & strcmp(column_names,'x')));
    these_y         = cell2mat(entries(3:end,strcmp(subject_colums,['s' SubjectIDs{isubject}]) & strcmp(column_names,'y')));
    these_z         = cell2mat(entries(3:end,strcmp(subject_colums,['s' SubjectIDs{isubject}]) & strcmp(column_names,'z')));

    data = ft_selectdata(cfg, data);

    % write out info of retained channels

    for ichannel = 1:numel(data.label)

        idx = strcmp(data.label{ichannel},these_labels);

        channs{isubject}(ichannel).names  = data.label{ichannel};
        channs{isubject}(ichannel).coords = [these_x(idx) these_y(idx) these_z(idx)];
    end

    %% Load subject file and change RT

    [numbers,strings] = xlsread([settings.base_path_castle,'well01_behavior_all.xls']);

    strings     = strings(2:end,:);

    if isnan(numbers(1,1))
        numbers         = numbers(2:end,:);
    end


    sel         = find(strcmp(strings(:,2),SubjectID));
    sel         = sel(1:numel(onsets));
    trls_enc    = strcmp(strings(sel,4),'encoding');
    trls_ret    = strcmp(strings(sel,4),'retrieval');


    Memory          = cell(size(strings(sel,12)));  % Convert Memory to a cell array of the same size
    Memory          = strings(sel,12);
    Memory          = Memory(trls_ret);
    idx_trial       = find(trls_ret);
    Memory(:, 2)    = num2cell(idx_trial);  % Assign idx_trial to the second column


    RT                          = numbers(sel,11);
    RT(RT==-1 & trls_enc==1)    = 3; % -1 no press in time - set to 3s at encoding

    if settings.full_enc_trial
        RT(trls_enc==1) = 3; % [optional] set all encoding to 3s
    end

    RT(RT==-1 & trls_ret==1) = 5; % -1 no press in time - set to 5s at encoding and 5s at retrieval

    trialinfo = [];

    for itrial = 1:numel(sel)

        trialinfo(itrial).SubjectID       = strings(sel(itrial),2);
        trialinfo(itrial).RunNumber       = numbers(sel(itrial),3);
        trialinfo(itrial).ExpPhase        = strings(sel(itrial),4);
        trialinfo(itrial).TrialNumber     = numbers(sel(itrial),5);
        trialinfo(itrial).EventNumber     = itrial;
        trialinfo(itrial).BlockType       = strings(sel(itrial),6);
        trialinfo(itrial).Subcat          = strings(sel(itrial),7);
        trialinfo(itrial).Word            = strings(sel(itrial),8);
        trialinfo(itrial).OldNew          = strings(sel(itrial),9);
        trialinfo(itrial).Response        = strings(sel(itrial),10);
        trialinfo(itrial).Memory          = strings(sel(itrial),12);
        trialinfo(itrial).RT              = RT(itrial);
    end

    sel     = [];
    strings = [];


    %% Realign trials based on stimulus onsets

    pretrig      = round(abs(settings.TOI_train(1)) * data.fsample);
    posttrig     = round(abs(settings.TOI_train(2)) * data.fsample);

    cfg          = [];
    cfg.trl      = [onsets'-pretrig onsets'+posttrig -pretrig*ones(numel(onsets),1)];

    data_stimuli = ft_redefinetrial(cfg,data);

    data_stimuli.trialinfo = {};
    for itrial = 1:numel(trialinfo)
        data_stimuli.trialinfo{itrial,1} = trialinfo(itrial);
    end

    %% [optional] resample

    if any(settings.do_resample)
        cfg             = [];
        cfg.resamplefs  = settings.do_resample;
        data_stimuli    = ft_resampledata(cfg, data_stimuli);
    end

    %% Time-lock data

    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.removemean  = 'no';


    data_stimuli    = ft_timelockanalysis(cfg, data_stimuli);

    fsample         = round(1/(data_stimuli.time(2)-data_stimuli.time(1)));


    %% [optional] Running mean filter

    if ~isempty(settings.smooth_win)
        data_stimuli.trial = smoothdata(data_stimuli.trial,3,'movmean',settings.smooth_win*fsample);
    end

    %% [optional] BL correct (us pre-stim baseline also for ripple-locked data)

    if settings.bs_correct == 1

        % training data
        bl_idx                  = nearest(data_stimuli.time,settings.bs_period(1)):nearest(data_stimuli.time,settings.bs_period(2));
        bldat                   = trimmean(data_stimuli.trial(:,:,bl_idx),settings.bs_trim,'round',3);
        blmat                   = repmat(bldat,[1 1 size(data_stimuli.trial,3)]);
        data_stimuli.trial      = data_stimuli.trial - blmat;

    end


    %% time-lock data

    cfg             = [];
    cfg.keeptrials  = 'yes';
    data_stimuli    = ft_timelockanalysis(cfg,data_stimuli);


    %%%%%%%%%%%%%%%%%%%%%%
    %%% CLASSIFICATION %%%
    %%%%%%%%%%%%%%%%%%%%%%

    trlinfo = cell2mat(data_stimuli.trialinfo);

    cfg = [];
    sel1    = ismember([trlinfo.ExpPhase],'retrieval');
    cfg.trials = sel1;
    cfg.latency = [0 6];
    data_retrieval = ft_selectdata(cfg,data_stimuli);

    cfg = [];
    sel1    = ismember([trlinfo.ExpPhase],'encoding');
    cfg.latency = [-.5 3];
    cfg.trials = sel1;
    data_encoding = ft_selectdata(cfg,data_stimuli);

    % Now divide retrieval from 1-5 second after cue onset at retrieval and
    % define that as a "ripple".
    % I want to take 1 second before to 1 second after the fictitious
    % ripple.

    % Sampling rate
    fsample = round(1/(data_retrieval.time(2) - data_retrieval.time(1)));

    % We only consider centres between 1 and 5 s after cue onset
    t_min = 1;   % s
    t_max = 5;   % s
    n_win = 60;  % number of pseudo-ripple positions you want

    % Candidate centre times in seconds
    centre_times = linspace(t_min, t_max, n_win);

    % Number of samples corresponding to 1 s on each side
    win_half_s   = 1;                       % 1 second before and after
    win_half_samp = round(win_half_s * fsample);
    n_time       = numel(data_retrieval.time);

    % Get indices in the time axis for each centre
    centre_idx = zeros(1, n_win);
    for iw = 1:n_win
        centre_idx(iw) = nearest(data_retrieval.time, centre_times(iw));
    end


    perf_correct     = [];
    perf_correct_var = [];
    perf_incorrect   = [];
    perf_incorrect_var = [];

    tic
    for iw = 1:n_win
        iw

        c_idx = centre_idx(iw);

        % Define start and end indices for ±1 s
        idx_start = c_idx - win_half_samp;
        idx_end   = c_idx + win_half_samp;


        % Convert back to seconds for ft_selectdata
        t_start = data_retrieval.time(idx_start);
        t_end   = data_retrieval.time(idx_end);

        cfg = [];
        cfg.latency = [t_start t_end];
        % keep all retrieval trials (you already selected retrieval in data_retrieval)
        data_test = ft_selectdata(cfg, data_retrieval);


        %%

        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%

        trlinfo = cell2mat(data_encoding.trialinfo);

        sel1    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'color');
        sel2    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'scene');

        dataToClassifyTraining      = cat(1,data_encoding.trial(sel1,:,:),data_encoding.trial(sel2,:,:));
        clabelTraining              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
        samples_train               = nearest(data_encoding.time,settings.TOI_train(1)):settings.timesteps_train*fsample:nearest(data_encoding.time,settings.TOI_train(2));

        dataToClassifyTraining      = dataToClassifyTraining(:,:,samples_train);

        %% 1. category PCA [colours vs. scenes; coarse]

        trlinfo = cell2mat(data_test.trialinfo);

        %% PCA on the data to get the eigenvalues that explain 90% of the variance.
        % Do it with a sliding window of 60ms with 100Hz)



        for icond = 1:2

            if icond == 1

                sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Memory],'SourceCorrect');
                sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & ismember([trlinfo.Memory],'SourceCorrect') ;



            elseif icond == 2

                sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color')& (ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));
                sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene')& (ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));



            end


            clabelTest              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
            dataToClassifyTest      = cat(1,data_test.trial(sel1,:,:),data_test.trial(sel2,:,:));


            sumsel1                 = sum(sel1);
            sumsel2                 = sum(sel2);

            fs              = 1/(data_test.time(2)-data_test.time(1));
            chunk_size      = round(settings.nu_time_points/((1/fs)*1000)); % Size of each chunk
            overlap_size    = round(chunk_size * settings.prc_overlap);
            num_chunks      = floor((size(dataToClassifyTest, 3) - overlap_size) / (chunk_size - overlap_size));
            TOI_ripple      = linspace(data_test.time(1), data_test.time(end),num_chunks);

            explained_variances     = [];
            how_much_variance       = [];
            accuracy                = [];
            ripple_to_decode        = [];
            for i = 1:num_chunks
                % Calculate the start and end indices of the current chunk
                start_idx   = (i - 1) * (chunk_size - overlap_size) + 1;
                end_idx     = start_idx + chunk_size - 1;

                % Extract data for the current chunk and reshape
                chunk_data = dataToClassifyTest(:, :, start_idx:end_idx);
                chunk_data = reshape(chunk_data, size(dataToClassifyTest, 1), []);

                if settings.smooth_before_dim == 1 % smooth NaNs through linear interpolation
                    for ismooth = 1:size(chunk_data, 1)

                        valid_indices           = ~isnan(chunk_data(ismooth, :));
                        chunk_data(ismooth, :)  = interp1(find(valid_indices), chunk_data(ismooth, valid_indices), 1:size(chunk_data, 2), 'linear', 'extrap');
                    end
                end

                % Perform PCA and compute explained variance
                [coefficients, ~, latent, ~, explained] = pca(chunk_data);

                % Compute explained variance
                explained_variance_pca = latent / sum(latent);

                % use a data-driven approach to get the first elbow point where
                % least variance is explained.
                curvature           = diff(diff(explained_variance_pca));
                [~, elbow_index]    = max(curvature);
                elbow_component     = elbow_index + 1; % Add 1 because of diff operation


                explained_variances(i, :) = elbow_component;

                how_much_variance(i,:) = sum(explained_variance_pca(1:elbow_component));

                % Do PCA inverse transformation for later decoding
                if settings.decode_components == 1
                    selected_components     = coefficients(:,1:elbow_component);
                    transformed_data        = chunk_data * selected_components;
                    reconstructed_data      = transformed_data * selected_components';
                    reconstructed_data      = reshape(reconstructed_data,[size(dataToClassifyTest(:, :, start_idx:end_idx))]);
                    ripple_to_decode(:,:,i) = nanmean(reconstructed_data,3); % take mean of those time points used in sliding window
                end
            end


            accuracy = explained_variances;


            if settings.decode_components == 1

                if settings.zscore_data4class
                    dataToClassifyTraining  = zscore(dataToClassifyTraining);
                    ripple_to_decode      = zscore(ripple_to_decode);
                end

                cfg                         = [];
                cfg.classifier              = settings.classifier;
                cfg.metric                  = settings.metric;
                [accuracy_dec, ~]           = mv_classify_timextime(cfg, dataToClassifyTraining, clabelTraining, ripple_to_decode, clabelTest);

            end

            if icond == 1

                perf_correct(iw,:) = accuracy;
                perf_correct_var(iw,:) = how_much_variance;

                if settings.decode_components == 1
                    perf{isubject}.dec.correct.accuracy = accuracy_dec;
                end
            elseif icond == 2

                perf_incorrect(iw,:) = accuracy;
                perf_incorrect_var(iw,:) = how_much_variance;


                if settings.decode_components == 1
                    perf{isubject}.dec.incorrect.accuracy   = accuracy_dec;
                end
            end

            % add info


        end

    end

    perf{isubject}.correct.accuracy     = perf_correct;
    perf{isubject}.correct.exl_var      = perf_correct_var;
    perf{isubject}.incorrect.accuracy    = perf_incorrect;
    perf{isubject}.incorrect.exl_var     = perf_incorrect_var;

    perf{isubject}.channelcount_test   = size(data_test.trial,2);
    perf{isubject}.channels_test       = data_test.label;
    perf{isubject}.time_train          = TOI_ripple;
    perf{isubject}.time_test           = TOI_ripple;



    %         clearvars -except perf settings isubject subjects SubjectIDs RT_all_subj_correct RT_all_subj_incorrect
end
toc

delete(gcp);

return
%% Stats and plots

load '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/scripts_to_publish/to_plot_for_paper/perf_dimensionality.mat'
explained_var_corr   = [];
explained_var_incorr = [];
figure;
toi_to_take = nearest(perf_dimensionality{1, 1}.time_test,-1):nearest(perf_dimensionality{1, 1}.time_test,1);
for participant = 1:numel(subjects)
    
    dtx = perf{1, participant}.correct.accuracy;

    noise=(rand(size(dtx,1),size(dtx,2))-0.5)./1.5;
    dtx=perf{1, participant}.correct.accuracy-noise;

    mean_shuffle = mean(dtx);
    std_shuffle = std(dtx);

    explained_var_corr(participant,:)      = (perf_dimensionality{1, participant}.correct.accuracy(toi_to_take,:)'-mean_shuffle)./std_shuffle;
    subplot(4,3,participant)
    plot(explained_var_corr(participant,:))

    empirical_correct(participant,:) = perf_dimensionality{1, participant}.correct.accuracy(toi_to_take,:);
  
end

correct_incorrect = {};

for isubject = 1:numel(subjects)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = linspace(-1,1,size(explained_var_corr,2));

    correct_incorrect{1}.individual(isubject,1,:)   = explained_var_corr(isubject,:);
    %     correct_incorrect{1}.individual(isubject,1,:)   = explained_var_corr(isubject,:);
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end

correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);


correct_incorrect{2}                = correct_incorrect{1};


for isubject = 1:numel(subjects)

    correct_incorrect{2}.individual(isubject,1,:)           = zeros(size(explained_var_corr(isubject,:)));
    %     correct_incorrect{2}.individual(isubject,1,:)           = explained_var_incorr(isubject,:);

end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);

%% FT stats (dimensionality)

xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.correcttail         = 'alpha'; % alpha prob no
cfg.neighbours          = [];
cfg.minnbchan           = 0;
cfg.avgovertime         = 'no';    % 'no' 'yes'
cfg.avgoverchan         = 'no';
cfg.computecritval      = 'yes';

cfg.numrandomization    = 'all';%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'individual';

nSub = size(correct_incorrect{1, 1}.individual  ,1);
% set up design matrix
design = zeros(2,2*nSub);
for i = 1:nSub
    design(1,i) = i;
end
for i = 1:nSub
    design(1,nSub+i) = i;
end
design(2,1:nSub)        = 1;
design(2,nSub+1:2*nSub) = 2;

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
[Fieldtripstats] = ft_timelockstatistics(cfg, correct_incorrect{:});
length(find(Fieldtripstats.mask))


%% plot significant vals



d = squeeze(empirical_correct);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

figure;
boundedline(correct_incorrect{1}.time,m,s,'b');
plot(correct_incorrect{1}.time,m,'k','linewidth',2);
hold on


stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 2;

plot(correct_incorrect{1}.time,sigline,'r','linewidth',4);

set(gca,'FontSize',16,'FontName','Arial')
xlabel('ripple time (s)')
ylabel('dimensionality difference')
set(gca,'TickDir','out')

axis tight
vline(0)


xlim([cfg.latency(1), cfg.latency(end)])



%% relate sigline to reaction time on a group level
d = squeeze(correct_incorrect{1}.individual);

RT = RT_correct;

dimensionality_change = mean(d(:,sigline==2),2);

figure;
scatter(RT, dimensionality_change, 'filled', 'MarkerFaceColor', '#0072BD');
xlabel('RT', 'FontSize', 12);
ylabel('Dimensionality Change', 'FontSize', 12);
title('Correlation', 'FontSize', 14);
grid on;
box on;
hold on;

% Fit a linear regression line
p = polyfit(RT, dimensionality_change, 1);
f = polyval(p, RT);
plot(RT, f, 'r-', 'LineWidth', 1.5);


% Add legend
legend('Data', 'Linear Fit', 'Location', 'best');

% Customize the plot appearance
set(gca, 'FontSize', 10);  % Set font size for axis labels
set(gcf, 'Color', 'w');    % Set background color of the figure to white


yfit = polyval(p, RT);
yresid = dimensionality_change - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(dimensionality_change)-1) * var(dimensionality_change);
rsq = 1 - SSresid/SStotal;
disp(['R-squared: ', num2str(rsq)]);

[rho_rt_dim, p_rt_dim] = corr(RT',dimensionality_change, 'type', 'Spearman');

%% correlate RT across time
data_to_correlate   = squeeze(correct_incorrect{1}.individual);
RT = RT_correct;

fs              = 1/(correct_incorrect{1, 1}.time(2)-correct_incorrect{1, 1}.time(1));
chunk_size      = round(100/((1/fs)*1000)); % Size of each chunk
overlap_size    = round(chunk_size * settings.prc_overlap);
num_chunks      = floor((size(correct_incorrect{1}.individual, 3) - overlap_size) / (chunk_size - overlap_size));
TOI_corr        = linspace(correct_incorrect{1, 1}.time(1), correct_incorrect{1, 1}.time(end),num_chunks);

rho_across_time = [];
p_across_time   = [];
for itime = 1:num_chunks
    start_idx   = (itime - 1) * (chunk_size - overlap_size) + 1;
    end_idx     = start_idx + chunk_size - 1;


    [rho_tmp, p_tmp]        = corr(mean(data_to_correlate(:,start_idx:end_idx),2), RT','type', 'spearman');
    rho_across_time(itime)  = rho_tmp;
    p_across_time(itime)    = p_tmp;
end

rho_across_time_perm = [];
p_across_time_perm   = [];
for nu_perm = 1:500
    rand_rt = randperm(12);

    for itime = 1:num_chunks
        start_idx   = (itime - 1) * (chunk_size - overlap_size) + 1;
        end_idx     = start_idx + chunk_size - 1;

        [rho_tmp, p_tmp]                    = corr(mean(data_to_correlate(rand_rt,start_idx:end_idx),2), RT','type', 'spearman');
        rho_across_time_perm(nu_perm,itime) = rho_tmp;
        p_across_time_perm(nu_perm,itime)   = p_tmp;
    end
end

xlimits_ind             = nearest(TOI_corr, -1):nearest(TOI_corr, 1);


alpha = 0.05;
z_threshold = norminv(1 - alpha/2);

zvalue_rho = (rho_across_time(xlimits_ind)-(mean(rho_across_time_perm(:,xlimits_ind))))./std(rho_across_time_perm(:,xlimits_ind));

xlimits = TOI_corr(xlimits_ind);
plot(xlimits, zvalue_rho,'linewidth', 3)

hold on
plot(xlimits, z_threshold * ones(size(xlimits)), 'r--');  % positive threshold
plot(xlimits, -z_threshold * ones(size(xlimits)), 'r--'); % negative threshold
below_threshold = zvalue_rho < -z_threshold;

scatter(xlimits(below_threshold), zvalue_rho(below_threshold), 'r', 'filled', 'MarkerFaceAlpha', 0.5)

ylim([-3 3])
xlim([-1 1])
set(gca,'FontSize',16)
set(gca,'TickDir','out')

hold off

xlabel('Ripple time (sec)')
ylabel('Z-transformed correlation')
title('Z-value of correlation Across Time')

legend('Z-value', 'Positive Threshold', 'Negative Threshold', 'Location', 'NorthEast')


%% plot t line

figure
plot(Fieldtripstats.time,Fieldtripstats.stat,'k','linewidth',2)
hold on
sigline05 = nan(1,numel(Fieldtripstats.prob));
sigline01 = nan(1,numel(Fieldtripstats.prob));
p05 = Fieldtripstats.prob < .05;
p01 = Fieldtripstats.prob < .01;

sigline05(p05) = Fieldtripstats.stat(p05);
sigline01(p01) = Fieldtripstats.stat(p01);

plot(Fieldtripstats.time,sigline05,'r','linewidth',5)
plot(Fieldtripstats.time,sigline01,'y','linewidth',2)

vline(0)
hline(0)

set(gca,'FontSize',8)
xlabel('time (sec)')
set(gca,'TickDir','out')

title('dimensionality');
xlim([cfg.latency(1), cfg.latency(end)])



%% correlate fine-grained decoding with dimensionality

% FT stats (dimensionality)

xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.correcttail         = 'alpha'; % alpha prob no
cfg.neighbours          = [];
cfg.minnbchan           = 0;
cfg.avgovertime         = 'no';    % 'no' 'yes'
cfg.avgoverchan         = 'no';
cfg.computecritval      = 'yes';

cfg.numrandomization    = 500;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'individual';

nSub = size(correct_incorrect{1, 1}.individual  ,1);
% set up design matrix
design = zeros(2,2*nSub);
for i = 1:nSub
    design(1,i) = i;
end
for i = 1:nSub
    design(1,nSub+i) = i;
end
design(2,1:nSub)        = 1;
design(2,nSub+1:2*nSub) = 2;

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
[Fieldtripstats] = ft_timelockstatistics(cfg, correct_incorrect{:});
length(find(Fieldtripstats.mask))


%% plot significant vals

d = squeeze(correct_incorrect{1}.individual-correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

figure;
boundedline(perf{1,1}.time_train,m,s,'k');
plot(perf{1,1}.time_train,m,'k','linewidth',2);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time,sigline,'r','linewidth',4);

set(gca,'FontSize',16,'FontName','Arial')
xlabel('ripple time (s)')
ylabel('dimensionality difference')
set(gca,'TickDir','out')

axis tight
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])

% load original data
load('mask_t_vals_fine_enc_ripple.mat')

mask_t_vals_fine          = squeeze(mask_t_vals_fine_enc_ripple);

% Load data from decoding analysis
load('correct_dec_fine_all.mat')
load('incorrect_dec_fine_all.mat')

correct_dec_all     = correct_dec_fine_all;
incorrect_dec_all   = incorrect_dec_fine_all;

time_dec = linspace(-.5,3,351);
idx_time = nearest(time_dec,-.2):nearest(time_dec,3);

correct_dec_all = correct_dec_all(:,idx_time,:);

incorrect_dec_all = incorrect_dec_all(:,idx_time,:);

for isubject = 1:size(settings.subjects,1)

    tmp_sub = [];
    tmp_sub = squeeze(correct_dec_all(isubject, :,:));
    tmp_sub(~mask_t_vals_fine) = NaN;
    correct_dec_all(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(incorrect_dec_all(isubject, :,:));
    tmp_sub(~mask_t_vals_fine) = NaN;
    incorrect_dec_all(isubject, :,:) = tmp_sub;
end

correct_dec_mean    = nanmean(nanmean(correct_dec_all,3),2);
incorrect_dec_mean  = nanmean(nanmean(incorrect_dec_all,3),2);

dimensionality_change = mean(d(:,sigline==0),2);

[r, p] = corr(dimensionality_change,correct_dec_mean-incorrect_dec_mean,'tail','right');

% Number of data points
n = numel(dimensionality_change);

% Convert r to Fisher's z-score
z = atanh(r);

effect_size = z * sqrt(n - 3);

disp(['Effect size (Cohen''s d): ' num2str(effect_size)]);

%% visualise included channels

% elec_size   = 15;
% transp      = 0.25;
% extracolor  = .5;
%
% epos_all    = [];
%
% for isubject = 1:numel(subjects)
%     epos_all = [epos_all;cell2mat({channs{isubject}.coords}')];
% end
%
% colvec = ones(1,size(epos_all,1));
%
% views = [-90 0;0 0;180 -90];
%
% for iview=1:size(views,1)
%     figure
%     plot_ecog(colvec, ...
%         fullfile(ftpath,'template/anatomy/'),...
%         epos_all,[-max(colvec) max(colvec)+extracolor], transp, views(iview,:), elec_size);
%     colorbar off
% end




%%