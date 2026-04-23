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


settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
subjects                    = {'CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK'};
SubjectIDs                  = {'01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK'};


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
settings.TOI_train              = [-.5 3];
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

    %% Find ripples

    [inData,alldat] = detect_ripples(data, settings);

    %% remove false positives from ripple data

    if settings.remove_falsepositives

        alldat = remove_false_positives(alldat);

    end

    %% delete co-occuring ripples

    if settings.remove_ripple_duplicates

        alldat = remove_ripple_duplicates(alldat);

    end


    %% Load data to realign based on cue onset encoding and based on ripples


    data_in                 = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',subjects{isubject}],'data','onsets_session');
    data                    = data_in.data;
    data_in                 = [];
    onsets_1                = onsets_session - data.time{1}(1)*inData.fsample;

    onsets                  = onsets_session - data.time{1}(1)*inData.fsample;

    data.sampleinfo         = 1+data.sampleinfo - data.sampleinfo(1);
    data.time{1}            = 1/data.fsample+data.time{1}-data.time{1}(1);

     %% Exclude artifacts and interpolate the NaNs with 3 before and after.
    for itrial= 1:numel(data.artifact)
        itrial
        tmp = [];
        tmp = data;
        tmp.trial = [];
        tmp.trial = data.trial{1,1}(itrial,:);
        tmp.artifact = data.artifact{1,itrial};
        for iarti = 1:size(tmp.artifact,1)
            tmp_arti = [];
            tmp_arti = tmp.artifact(iarti,:);
            tmp.trial(1,tmp_arti(1):tmp_arti(2)) = NaN;
            
        end
        replace = find(isnan(tmp.trial));
        for iarti = 2:numel(replace)
            if replace(1,iarti) - replace(1,iarti-1) < 5
                tmp.trial(1,replace(1,iarti-1):replace(1,iarti)) = NaN;
            end
        end
        replace = [];
        replace = find(isnan(tmp.trial));
        for iarti = 1:numel(replace) % add mean of 3 around NaN;
            if replace(1,iarti) < 4 || replace(1,iarti)>numel(tmp.trial)-3
                mean_of = nanmean(tmp.trial(1,replace(1,:))); % if beginning or end of vector, take mean of the whole trial
                tmp.trial(1,replace(1,iarti)) = mean_of;
            else
                mean_of = nanmean(tmp.trial(1,replace(1,iarti)-3:replace(1,iarti)+3));
                tmp.trial(1,replace(1,iarti)) = mean_of;
            end
        end
        data.trial{1,1}(itrial,:) =  tmp.trial;
        
    end


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

    %% create onset matrices (remove last 250ms to ensure ripple event in trial)

    onsetmat = [onsets; onsets+(RT'.*data.fsample)-(settings.time_to_excl_RT*data.fsample)]';

    %% Extract ripples (optional to select long and short duration ripples)

    trl_ripple  = [];
    cnt         = 0;

    for ichannel = 1:numel(alldat)

        evs     = alldat{ichannel}.evtIndiv.maxTime;
        envSum  = alldat{ichannel}.evtIndiv.envSum;

        ripple_dur = alldat{ichannel}.evtIndiv.duration;

        duration_sel = logical(ones(1,numel(ripple_dur)));
       
        evs     = evs(duration_sel);
        envSum  = envSum(duration_sel);

        %% for each detected ripple, find the corresponding trial

        for iripple = 1:numel(evs)

            this_event = evs(iripple) >= onsetmat(:,1) & evs(iripple) <= onsetmat(:,2);

            if any(this_event)
                cnt=cnt+1;

                trl_ripple(cnt,1) = find(this_event);   % note down corresponding event number
                trl_ripple(cnt,2) = evs(iripple);       % note down ripple sample
                trl_ripple(cnt,3) = (evs(iripple) - onsetmat(this_event,1))/data.fsample; % note down time of ripple in trial
                trl_ripple(cnt,4) = ichannel; % note down channel
                trl_ripple(cnt,5) = envSum(iripple);

            end
        end

    end

    trl_ripple = sortrows(trl_ripple,1);

    %% restrict ripple data to retrieval

    trl_ripple = trl_ripple(ismember(trl_ripple(:,1),find(trls_ret)),:);

    %% [optional] do surrogates by taking time of ripple from other trial (for all or only for correct trials)


    correct_mem         = cell2mat(Memory(ismember(Memory(:,1),'SourceCorrect'),2));
    idx_correct         = ismember(trl_ripple(:,1), correct_mem);
    trl_ripple_correct  = trl_ripple(idx_correct,:);

    if settings.do_surrogates == 1 % 1 for all trials

        randtrials      = circshift(1:size(trl_ripple,1),1);
        tmp             = trl_ripple(:,2)-round(trl_ripple(:,3).*data.fsample); % find cue onset
        tmp             = tmp+round(trl_ripple(randtrials,3).*data.fsample); % add another ripple's event time
        trl_ripple(:,2) = tmp;
        trl_ripple(:,3) = trl_ripple(randtrials,3);

    elseif settings.do_surrogates == 2 % 2 for only correct trials

        randtrials_correct      = circshift(1:size(trl_ripple_correct,1),-1);
        tmp_correct             = trl_ripple_correct(:,2) - round(trl_ripple_correct(:,3) .* data.fsample); % find cue onset
        tmp_correct             = tmp_correct + round(trl_ripple_correct(randtrials_correct,3) .* data.fsample); % add shuffled ripple event time
        trl_ripple_correct(:,2) = tmp_correct;
        trl_ripple_correct(:,3) = trl_ripple_correct(randtrials_correct,3);

        trl_ripple(idx_correct,:) = trl_ripple_correct; % add to original structure with only correct trials swapped

    end

    %% Create trial structure around ripples

    pretrig      = round(abs(settings.TOI_test(1)) * data.fsample); % enough time to baseline correct later
    posttrig     = round(abs(settings.TOI_test(2)) * data.fsample);

    cfg          = [];
    cfg.trl      = [trl_ripple(:,2)-pretrig trl_ripple(:,2)+posttrig -pretrig*ones(size(trl_ripple,1),1)];

    data_ripples = ft_redefinetrial(cfg,data);

    % add trial info for each ripple trial, accounting for multiple ripples
    % per trial
    data_ripples.trialinfo = [];

    for itrial = 1:numel(data_ripples.trial)

        corresponding_trialinfo = find([trialinfo.EventNumber]==trl_ripple(itrial,1));

        trl_info                         = trialinfo(corresponding_trialinfo);
        trl_info.sample_ripple           = trl_ripple(itrial,2);
        trl_info.time_ripple             = trl_ripple(itrial,3);
        trl_info.channel                 = trl_ripple(itrial,4);
        trl_info.name_channel            = {alldat{trl_ripple(itrial,4)}.evtIndiv.label};
        trl_info.envSum                  = trl_ripple(itrial,5);

        data_ripples.trialinfo{itrial,1} = trl_info;
    end

    %% [optional] if there are multiple ripples per trial - pick the one ripple with highest activity (captured in the summed envelope metric)

    if settings.solo_ripple

        tmp_trl_info = data_ripples.trialinfo;
        % find the trials in which there were more than one ripple
        trlinfo         = cell2mat(data_ripples.trialinfo);
        ripple_trial    = [trlinfo.EventNumber];
        envSum          = [trlinfo.envSum];
        idx_unique      = unique(ripple_trial);

        sel     = [];
        counter = 1;

        for itrial = 1:numel(idx_unique)

            idx                 = find(ripple_trial==idx_unique(itrial));
            [~,max_effect]      = max(envSum(idx));
            sel(counter)        = idx(max_effect);

            counter             = counter+1;
        end

        cfg             = [];
        cfg.trials      = sel;
        data_ripples    = ft_selectdata(cfg, data_ripples);

    end

    %% [optional] pick ripples in a specific time window
    if any(settings.ripple_latency)

        trlinfo         = cell2mat(data_ripples.trialinfo);
        sel             = [trlinfo.time_ripple] > settings.ripple_latency(1) & [trlinfo.time_ripple] < settings.ripple_latency(2);

        cfg             = [];
        cfg.trials      = sel;
        data_ripples    = ft_selectdata(cfg, data_ripples);

    end

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
        data_ripples    = ft_resampledata(cfg, data_ripples);
        data_stimuli    = ft_resampledata(cfg, data_stimuli);
    end

    %% Time-lock data

    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.removemean  = 'no';

    data_ripples    = ft_timelockanalysis(cfg, data_ripples);
    data_stimuli    = ft_timelockanalysis(cfg, data_stimuli);

    fsample         = round(1/(data_ripples.time(2)-data_ripples.time(1)));

    %% [optional] Running mean filter

    if ~isempty(settings.smooth_win)
        data_ripples.trial = smoothdata(data_ripples.trial,3,'movmean',settings.smooth_win*fsample);
        data_stimuli.trial = smoothdata(data_stimuli.trial,3,'movmean',settings.smooth_win*fsample);
    end

    %% [optional] BL correct (us pre-stim baseline also for ripple-locked data)

    if settings.bs_correct == 1

        % training data
        bl_idx                  = nearest(data_stimuli.time,settings.bs_period(1)):nearest(data_stimuli.time,settings.bs_period(2));
        bldat                   = trimmean(data_stimuli.trial(:,:,bl_idx),settings.bs_trim,'round',3);
        blmat                   = repmat(bldat,[1 1 size(data_stimuli.trial,3)]);
        data_stimuli.trial      = data_stimuli.trial - blmat;

        % testing data
        trlinfo     = cell2mat(data_ripples.trialinfo);
        orig_events = [trlinfo.EventNumber];

        bl4ripples  = nan(size(data_ripples.trial,1),size(bldat,2));

        for iripple = 1:size(data_ripples.trial,1)
            bl4ripples(iripple,:) = bldat(orig_events(iripple),:);
        end
        blmat = repmat(bl4ripples,[1 1 size(data_ripples.trial,3)]);

        data_ripples.trial = data_ripples.trial - blmat;

    end


    %% time-lock data

    cfg             = [];
    cfg.keeptrials  = 'yes';
    data_stimuli    = ft_timelockanalysis(cfg,data_stimuli);
    data_ripples    = ft_timelockanalysis(cfg,data_ripples);

    %% calculate time of ripple and RT for those trials

    ripple_tmp  = [];
    RT_tmp      = [];
    for itrial = 1:numel(data_ripples.trialinfo)

        ripple_tmp(itrial)  = data_ripples.trialinfo{itrial, 1}.time_ripple;
        RT_tmp(itrial)      = data_ripples.trialinfo{itrial, 1}.RT;

    end
    trlinfo = cell2mat(data_ripples.trialinfo);
    sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.Memory],'SourceCorrect');
    sel2 = ismember([trlinfo.ExpPhase],'retrieval') &(ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));

    RT_all_subj_correct{isubject} = RT_tmp(sel1);
    RT_all_subj_incorrect{isubject} = RT_tmp(sel2);
    ripple_time_correct{isubject} = ripple_tmp(sel1);
    ripple_time_incorrect{isubject} = ripple_tmp(sel2);


    %%

    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% PCA %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%

    trlinfo = cell2mat(data_stimuli.trialinfo);

    sel1    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'color');
    sel2    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'scene');

    dataToClassifyTraining      = cat(1,data_stimuli.trial(sel1,:,:),data_stimuli.trial(sel2,:,:));
    clabelTraining              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
    samples_train               = nearest(data_stimuli.time,settings.TOI_train(1)):settings.timesteps_train*fsample:nearest(data_stimuli.time,settings.TOI_train(2));

    dataToClassifyTraining      = dataToClassifyTraining(:,:,samples_train);

    %% 1. category PCA [colours vs. scenes; coarse]

    trlinfo = cell2mat(data_ripples.trialinfo);

    %% PCA on the data to get the eigenvalues that explain 90% of the variance.
    % Do it with a sliding window of 60ms with 100Hz)





    for icond = 1:2
        if icond == 1
            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Memory],'SourceCorrect');
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & ismember([trlinfo.Memory],'SourceCorrect');
            perf{isubject}.correct.trl_num_test = [sum(sel1) sum(sel2)];
        elseif icond == 2
            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color')& (ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene')& (ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));
            perf{isubject} .incorrect.trl_num_test = [sum(sel1) sum(sel2)];
        end

        clabelTest = cat(1, 1*ones(sum(sel1), 1), 2*ones(sum(sel2), 1));
        dataToClassifyTest = cat(1, data_ripples.trial(sel1,:,:), data_ripples.trial(sel2,:,:));
        fs = 1 / (data_ripples.time(2) - data_ripples.time(1));
        chunk_size = round(settings.nu_time_points / ((1 / fs) * 1000)); % Size of each chunk
        overlap_size = round(chunk_size * settings.prc_overlap);
        num_chunks = floor((size(dataToClassifyTest, 3) - overlap_size) / (chunk_size - overlap_size));
        TOI_ripple = linspace(data_ripples.time(1), data_ripples.time(end), num_chunks);
        explained_variances = [];
        accuracy = [];
        ripple_to_decode = [];

        for i = 1:num_chunks
            % Calculate the start and end indices of the current chunk
            start_idx = (i - 1) * (chunk_size - overlap_size) + 1;
            end_idx = start_idx + chunk_size - 1;

            % Extract data for the current chunk and reshape
            chunk_data = dataToClassifyTest(:, :, start_idx:end_idx);
            chunk_data = reshape(chunk_data, size(dataToClassifyTest, 1), []);

            if settings.smooth_before_dim == 1
                % Smooth NaNs through linear interpolation
                for ismooth = 1:size(chunk_data, 1)
                    valid_indices = ~isnan(chunk_data(ismooth, :));
                    chunk_data(ismooth, :) = interp1(find(valid_indices), chunk_data(ismooth, valid_indices), 1:size(chunk_data, 2), 'linear', 'extrap');
                end
            end

            % Perform PCA and compute explained variance
            [coefficients, ~, latent, ~, explained] = pca(chunk_data);
            explained_variance_pca = latent / sum(latent);

            % Fit a power-law distribution
            pl_fit = fit((1:length(explained_variance_pca))', explained_variance_pca, 'power1');

            % Calculate distance to the elbow
            curvature = diff(diff(explained_variance_pca));
            [~, elbow_index] = max(curvature);
            elbow_component = elbow_index + 1; % Add 1 because of diff operation
            distance_to_elbow = sum((explained_variance_pca(elbow_component) - feval(pl_fit, elbow_component)).^2);


            signed_distance_to_elbow = explained_variance_pca(elbow_component) - feval(pl_fit, elbow_component);
            distance_to_elbow = sum(signed_distance_to_elbow);  % Sum of signed differences

            explained_variances(i, :) = distance_to_elbow;

        end

        accuracy = explained_variances;
        

        if icond == 1
            perf{isubject}.correct.accuracy = accuracy;
           
        elseif icond == 2
            perf{isubject}.incorrect.accuracy = accuracy;
           
        end

        % Add info
        perf{isubject}.channelcount_test = size(data_ripples.trial, 2);
        perf{isubject}.channels_test = data_ripples.label;
        perf{isubject}.time_train = TOI_ripple;
        perf{isubject}.time_test = TOI_ripple;
    end

    %         clearvars -except perf settings isubject subjects SubjectIDs RT_all_subj_correct RT_all_subj_incorrect
end
toc

delete(gcp);

return
%% Stats and plots

RT_correct          = cellfun(@mean, RT_all_subj_correct);
RT_max_correct      = cellfun(@max, RT_all_subj_correct);
RT_min_correct      = cellfun(@min, RT_all_subj_correct);

RT_incorrect        = cellfun(@mean, RT_all_subj_incorrect);
RT_max_incorrect    = cellfun(@max, RT_all_subj_incorrect);
RT_min_incorrect    = cellfun(@min, RT_all_subj_incorrect);

perf{1}.RT.correct      = RT_all_subj_correct;
perf{1}.RT.incorrect    = RT_all_subj_incorrect;


RT = RT_correct;

ripple_time_mean    = cellfun(@mean, ripple_time_correct);
ripple_time_max     = cellfun(@max, ripple_time_correct);
ripple_time_min     = cellfun(@min, ripple_time_correct);

perf{1}.ripple_times.correct    = ripple_time_correct;
perf{1}.ripple_times.incorrect  = ripple_time_incorrect;

ripples_to_plot = [];
rt_to_plot      = [];
for participant = 1:numel(subjects)

    ripples_to_plot = [ripples_to_plot, ripple_time_correct{participant}];
    rt_to_plot      = [rt_to_plot, RT_all_subj_correct{participant}];

end

delay_ripple_rt = rt_to_plot-ripples_to_plot;

figure;
subplot(3,1,1)
nhist(ripples_to_plot','proportion','color',settings.colour_scheme(8,:))
title('Time of ripples')
xlabel('time of ripples')
ylabel('proportion')
set(gca,'FontSize',14)
set(gca,'TickDir','out')
title(sprintf('Time of ripples, median = %.2fms',median(ripple_time_mean*1000)),'interpreter','none')
subplot(3,1,2)
nhist(rt_to_plot','proportion','color',settings.colour_scheme(8,:))
title('Reaction time in trials of ripples')
xlabel('Reaction time')
ylabel('proportion')
set(gca,'FontSize',14)
set(gca,'TickDir','out')
title(sprintf('Reaction time in trials of ripples, median = %.2fms',median(RT*1000)),'interpreter','none')
subplot(3,1,3)
nhist(delay_ripple_rt','proportion','color',settings.colour_scheme(8,:))
xlabel('Delay ripple RT')
ylabel('proportion')
set(gca,'FontSize',14)
set(gca,'TickDir','out')
title(sprintf('Delay ripples RT, median = %.2fms',median(delay_ripple_rt*1000)),'interpreter','none')


for participant = 1:numel(subjects)
    trl_num_test                = perf{participant}.correct.trl_num_test;
    trl_correct(participant)    = sum(trl_num_test);
    trl_num_test                = perf{participant}.incorrect.trl_num_test;
    trl_incorrect(participant)  = sum(trl_num_test);
end

data_nu_trl = {};
data_nu_trl{1,1} = trl_correct;
data_nu_trl{2,1} = trl_incorrect;

[~,p_val,~,stats] = ttest(trl_correct,trl_incorrect)

figure;
h = rm_raincloud(data_nu_trl, [settings.colour_scheme(6,:)],0,'ks',[],settings.colour_scheme);

h.p{1, 1}.FaceColor         = settings.colour_scheme(1,:);
h.s{1, 1}.MarkerFaceColor   = settings.colour_scheme(1,:);
h.m(1, 1).MarkerFaceColor    = settings.colour_scheme(1,:);
h.p{2, 1}.FaceColor         = settings.colour_scheme(10,:);
h.s{2, 1}.MarkerFaceColor   = settings.colour_scheme(10,:);
h.m(2, 1).MarkerFaceColor    = settings.colour_scheme(10,:);


hold on
title(sprintf('number of trials for the two conditions, t-stat = %.2f', stats.tstat))
set(gca,'TickDir','out')
xlabel('number of trials')
yticklabels({sprintf('incorrect %d',mean(trl_incorrect)), sprintf('correct %d',floor(mean(trl_correct)))})
ylabel('condition')
set(gca,'FontSize',20)
axis tight

[~,p_val,~,stats] = ttest(trl_correct,trl_incorrect)

correct_incorrect = {};
correct_baseline = {};

for isubject = 1:numel(subjects)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

     mean_zscore = mean([perf{isubject}.correct.accuracy;perf{isubject}.incorrect.accuracy]);
    stdzscore = std([perf{isubject}.correct.accuracy;perf{isubject}.incorrect.accuracy]);

    tmp                                                     = (perf{isubject}.correct.accuracy-mean_zscore)./stdzscore;


    correct_incorrect{1}.individual(isubject,1,:)   =tmp;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

    correct_incorrect{2}.label               = {'Channels'};
    correct_incorrect{2}.time                = perf{1,1}.time_train;

    tmp                                                     = (perf{isubject}.incorrect.accuracy-mean_zscore)./stdzscore;


     correct_incorrect{2}.individual(isubject,1,:)           = tmp;
    correct_incorrect{2}.dimord              = 'subj_chan_time';

     correct_baseline               = correct_incorrect;

     correct_baseline{2}.label               = {'Channels'};
    correct_baseline{2}.time                = perf{1,1}.time_train;
    correct_baseline{2}.individual(isubject,1,:)   =zeros(size(perf{isubject}.correct.accuracy));
    correct_baseline{2}.dimord              = 'subj_chan_time';


end

correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);

correct_baseline{1}.avg            = squeeze(correct_baseline{1}.individual);
correct_baseline{2}.avg            = squeeze(correct_baseline{2}.individual);



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

d = squeeze(correct_incorrect{1}.individual)-squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

figure;
boundedline(perf{1,1}.time_train,m,s,'b');
plot(perf{1,1}.time_train,m,'k','linewidth',2);
hold on


stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time,sigline,'r','linewidth',4);

set(gca,'FontSize',16,'FontName','Arial')
xlabel('ripple time (s)')
ylabel('signed distance to power-law distribution')
set(gca,'TickDir','out')

axis tight
vline(0)
hline(0)


xlim([cfg.latency(1), cfg.latency(end)])


%% relate sigline to reaction time on a group level
d = squeeze(correct_incorrect{1}.individual);

RT = RT_correct;

dimensionality_change = mean(d(:,sigline==0),2);

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



%% FT stats (decoding)

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];
cfg.frequency           = [-.2 freqaxis(end)]; % encoding time
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
cfg.computecritval      = 'yes';

cfg.numrandomization    = 1000;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = numel(subjects);
% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

[Fieldtripstats] = ft_freqstatistics(cfg, correct{:}, incorrect{:});

length(find(Fieldtripstats.mask==1))

% plot (significant) t vals

stats_time = nearest(perf{1, 1}.time_test,cfg.latency(1)):nearest(perf{1, 1}.time_test,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);

figure;
imagesc(...
    perf{1, 1}.time_test(stats_time),...
    freqaxis(stats_freq),...
    tvals);
colormap(jet)
caxis([-4 4])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',12)
xlabel('ripple time (sec)')
ylabel('encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title(sprintf('coarse-grained classifier, correct vs. incorrect \n%s correction, %s',cfg.correctm,cfg.method));

if any(Fieldtripstats.mask(:))
    plot_contour(perf{1, 1}.time_test(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats.mask)))
end


%% Correlate with original data (still debating whether we want to include this or not)

% get size of the different decoding analyses

xlimits = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);

% limit decoding from reconstructed data
correct_dec_rec     = correct_dec(:, :, xlimits);
incorrect_dec_rec   = incorrect_dec(:, :, xlimits);
% load original data
load('tvals_coarse_enc_ripple.mat')
load('mask_t_vals_coarse_enc_ripple.mat')

tvals_single_ripple_coarse  = tvals_coarse_enc_ripple;
mask_t_vals_coarse          = squeeze(mask_t_vals_coarse_enc_ripple);

% Load data from decoding analysis
load('correct_dec_all.mat')
load('incorrect_dec_all.mat')

correct_dec_all = correct_dec_all;
incorrect_dec_all = incorrect_dec_all;

[rows_large, cols_large] = size(tvals_single_ripple_coarse);

% Get the size of the smaller matrix
[rows_small, cols_small] = size(tvals);

% Calculate the number of columns to add
cols_to_add = cols_large - cols_small;

% Pad the smaller matrix with zeros to match the size of the larger one
% (only a few data points)
tvals_reconstructed_resized                 = [tvals, zeros(rows_small, cols_to_add)];
correct_dec_rec(:,:,end:end+cols_to_add)    = 0;
incorrect_dec_rec(:,:,end:end+cols_to_add)  = 0;
TOI_to_plot                                 = linspace(-1,1, size(tvals_single_ripple_coarse,2));
FOI_to_plot                                 = correct{1, 1}.freq;


[r, p] = corr(tvals_single_ripple_coarse(logical(mask_t_vals_coarse)),tvals_reconstructed_resized(logical(mask_t_vals_coarse)));

r
p



for isubject = 1:size(settings.subjects,1)
    tmp_sub = [];
    tmp_sub = squeeze(correct_dec_rec(isubject, :,:));
    tmp_sub(~mask_t_vals_coarse) = NaN;
    correct_dec_rec(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(correct_dec_all(isubject, :,:));
    tmp_sub(~mask_t_vals_coarse) = NaN;
    correct_dec_all(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(incorrect_dec_rec(isubject, :,:));
    tmp_sub(~mask_t_vals_coarse) = NaN;
    incorrect_dec_rec(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(incorrect_dec_all(isubject, :,:));
    tmp_sub(~mask_t_vals_coarse) = NaN;
    incorrect_dec_all(isubject, :,:) = tmp_sub;
end

correct_dec_mean    = nanmean(nanmean(correct_dec_all,3),2);
incorrect_dec_mean  = nanmean(nanmean(incorrect_dec_all,3),2);

correct_dec_rec_mean    = nanmean(nanmean(correct_dec_rec,3),2);
incorrect_dec_rec_mean  = nanmean(nanmean(incorrect_dec_rec,3),2);

tvals_single_ripple_coarse(~mask_t_vals_coarse)     = 0;
tvals_reconstructed_resized(~mask_t_vals_coarse)    = 0;

figure;
subplot(2,1,1)
imagesc(TOI_to_plot, FOI_to_plot, tvals_single_ripple_coarse)
axis xy
title('original')
set(gca,'FontSize',12)
xlabel('ripple time (sec)')
ylabel('encoding time (sec)')
set(gca,'TickDir','out')
clim([0 4])
vline(0)
hcb = colorbar('Location','EastOutside');
set(gca,'FontSize',20)
hold on
subplot(2,1,2)
imagesc(TOI_to_plot, FOI_to_plot, tvals_reconstructed_resized)
axis xy
title('reconstructed')
vline(0)
set(gca,'FontSize',12)
xlabel('ripple time (sec)')
ylabel('encoding time (sec)')
set(gca,'TickDir','out')
clim([0 4])
hcb = colorbar('Location','EastOutside');
% ax = text(.8,3.8, sprintf('correlation = %0.2f', r))
% ax.FontSize = 16;
% ax = text(.8,3.5, sprintf('p-value = %0.3f', p))
% ax.FontSize = 20;

set(gca,'FontSize',20)

% Number of data points
n = tvals_single_ripple_coarse(:);
n = numel(find(n));

% Convert r to Fisher's z-score
z = atanh(r);

% Calculate effect size (Cohen's d)
effect_size = z * sqrt(n - 3);

% Display effect size
disp(['Effect size (Cohen''s d): ' num2str(effect_size)]);


%% correlate fine-grained decoding with dimensionality

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

% Calculate effect size (Cohen's d)
effect_size = z * sqrt(n - 3);

% Display effect size
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









