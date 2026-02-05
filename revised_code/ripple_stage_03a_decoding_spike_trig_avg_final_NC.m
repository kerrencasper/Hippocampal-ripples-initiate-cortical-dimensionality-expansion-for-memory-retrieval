%%
% [ripple_stage_03a_decoding] find ripples in hippocampal channels,
% extract and realign data based on ripple events.
% Use encoding data to train LDA and test at ripple aligned
% data.
%                  Bernhard Staresina [bernhard.staresina@psy.ox.ac.uk]
%                  Casper Kerren      [kerren@cbs.mpg.de]


clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

[~,ftpath]=ft_version;

%% path settings
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects            = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW','09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

settings.data_dir           = [settings.base_path_castle,'preprocessing/channel_removal_all_channels/common_trimmed_average/'];
settings.data_dir_art       = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];

settings.save_dir           = [settings.base_path_castle,'output_data/decoding/'];
settings.data_dir_channels  = [settings.base_path_castle,'ripple_project_publication_for_replication/templates'];
settings.anatomy_dir        = [settings.base_path_castle,'ripple_project_publication_for_replication/additional_analyses/visualisation/'];
settings.AAL_dir            = fullfile(settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions/AAL3');
settings.SPM_dir            = fullfile(settings.base_path_castle,'/ripple_project_publication_for_replication/subfunctions/spm12');

settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};

subjects                    = {'CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK'};
SubjectIDs                  = {'01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK'};
settings.nu_rep             = [1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1];

settings.healthyhemi        = {'R' 'LR' 'R' 'L' 'L' 'R' 'R' 'R' 'R' 'R' 'R' 'LR'};
addpath('/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/scripts_to_publish')
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/MVPA-Light-master'))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))
%% pre-decoding settings
% ripple extraction

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1; % remove co-occuring ripples
settings.time_to_excl_RT            = .25; % exclude last 250 ms of trials, to make sure ripple event was in trial (Plotting -1 to 1 around ripple. Max ripple is 500ms, meaning it starts 250ms before 0.
settings.solo_ripple                = 1; % pick one ripple per trial if multiple ripple events are found
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur, eg [.5 1.5]

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

settings.TOI_test               = [-1 1]; % time around ripple
settings.timesteps_test         = ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample); % in s

settings.classifier             = 'lda';
settings.metric                 = 'auc';

%% channel settings

settings.channel_sel            = 1; % 1 exclude hippo, 2 only hippo, 3 all channels

%% start for loop

timeaxis = settings.TOI_test(1):settings.timesteps_test:settings.TOI_test(2);
freqaxis = settings.TOI_train(1):settings.timesteps_train:settings.TOI_train(2);

perf    = cell(1,numel(subjects));
channs  = cell(1,numel(subjects));

numWorkers = 6; %

parpool('local', numWorkers);

tic
parfor isubject = 1:numel(subjects)


    fprintf('processing subject %01d/%02d\n',isubject,numel(subjects));

    %% find spikes

    data_in         = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_',num2str(settings.subjects(isubject,:))],'data','onsets_session');
    data            = data_in.data;
    onsets_session  = data_in.onsets_session;
    around          = [5 10]; % how much time to include before first onset and after last onset


    %% restrict to relevant data (note that if you want the actual onset in ms after having cut around the relevant data onsets will have to be adjusted henceforth, eg onsets  = onsets_session - data.time{1}(1)*data.fsample;

    cfg         = [];
    cfg.latency = [onsets_session(1)/data.fsample-around(1) onsets_session(end)/data.fsample+around(2)];
    data        = ft_selectdata(cfg,data);

    sampleinfo_tmp = data.sampleinfo(1);

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

    if isubject == 8
        cfg             = [];
        cfg.prewindow   = 3;
        cfg.postwindow  = 3;
        data            = ft_interpolatenan(cfg,data);
    end

    X = data.trial{1}';      % Nsamples x Nchannels
    fs = data.fsample;
    exclusion_range = round(0.01 * fs);


    [xnew, spkindx] = despike_dilate_peak_spike_old(X,2, exclusion_range, false, 0);


    fs             = data.fsample;
    [Nsamp, Nchan] = size(X);
    [spk_samp, spk_chan] = ind2sub(size(X), spkindx);


    refrac_ms = 100;                 % e.g., 50 ms
    refrac_samp = round(refrac_ms/1000 * fs);

    keep = true(size(spk_samp));
    s = [];
    for c_idx = 1:Nchan
        idx = find(spk_chan == c_idx);
        [~,ord] = sort(spk_samp(idx));
        idx = idx(ord);

        last = -inf;
        for j = 1:numel(idx)
            s = spk_samp(idx(j));
            if s - last < refrac_samp
                keep(idx(j)) = false;
            else
                last = s;
            end
        end
    end

    spk_samp = spk_samp(keep);
    spk_chan = spk_chan(keep);
    spkindx  = spkindx(keep);

    % counts = accumarray(spk_chan, 1, [size(X,2) 1]);
    % disp(counts')


    win_sec   = [-0.1 0.1];
    win_samp  = round(win_sec * fs);
    lags      = win_samp(1):win_samp(2);
    nLags     = numel(lags);
    t         = lags / fs;

    % baseline window (pre-spike)
    baseline_sec  = [-0.1 -0.02];
    baseline_samp = round(baseline_sec * fs);
    base_idx      = baseline_samp(1):baseline_samp(2);
    base_mask     = ismember(lags, base_idx);

    STA_chan = zeros(nLags, Nchan);
    count    = zeros(Nchan, 1);
    c = [];
    s = [];
    for k = 1:numel(spk_samp)
        s = spk_samp(k);
        c = spk_chan(k);

        if s + lags(1) < 1 || s + lags(end) > Nsamp
            continue
        end

        seg = X(s + lags, c);

        [b,a] = butter(2, 5/(fs/2), 'high');   % 5 Hz HP
        seg = filtfilt(b,a, seg);

        % epoch-wise baseline correction
        b   = mean(seg(base_mask), 'omitnan');
        seg = seg - b;

        %        seg = detrend(seg, 'linear');

        STA_chan(:, c) = STA_chan(:, c) + seg;
        count(c)       = count(c) + 1;
    end

    for c_idx = 1:Nchan
        if count(c_idx) > 0
            STA_chan(:, c_idx) = abs(STA_chan(:, c_idx) ./ count(c_idx));
        else
            STA_chan(:, c_idx) = NaN;
        end
    end

    valid_ch = count > 0;
    STA_mean = nanmean(STA_chan(:, valid_ch), 2);

    perf{isubject}.STA              = STA_chan(:, valid_ch);
    perf{isubject}.STA_mean         = STA_mean;
    perf{isubject}.STA_mean_time    = t;

    %% Load data to realign based on cue onset encoding and based on spikes


    data_in                 = load([settings.data_dir_art,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',subjects{isubject}],'data','onsets_session');
    data                    = data_in.data;
    data_in                 = [];
    onsets                  = onsets_session - data.time{1}(1)*data.fsample;
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

    spkindx = spkindx-sampleinfo_tmp;

    fs                 = data.fsample;
    [Nsamp, Nchan]     = size(X);

    % Convert linear indices → (sample, channel)
    [spk_samp, spk_chan] = ind2sub(size(X), spkindx);

    unique_chan = unique(spk_chan);

    for ichannel = 1:numel(unique_chan)

        this_chan = spk_chan==ichannel;
        this_spike_chan = spk_samp(this_chan);


        %% for each detected ripple, find the corresponding trial

        for ispike = 1:numel(this_spike_chan)

            this_event = this_spike_chan(ispike) >= onsetmat(:,1) & this_spike_chan(ispike) <= onsetmat(:,2);

            if any(this_event)
                cnt=cnt+1;

                trl_ripple(cnt,1) = find(this_event);   % note down corresponding event number
                trl_ripple(cnt,2) = this_spike_chan(ispike);       % note down spike sample
                trl_ripple(cnt,3) = (this_spike_chan(ispike) - onsetmat(this_event,1))/data.fsample; % note down time of spike in trial
                trl_ripple(cnt,4) = ichannel; % note down channel


            end
        end

    end

    if isempty(trl_ripple)
        continue;
    end

    trl_ripple = sortrows(trl_ripple,1);

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
        trl_info.name_channel            = {channels{trl_ripple(itrial,4)}};

        data_ripples.trialinfo{itrial,1} = trl_info;
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



    %%%%%%%%%%%%%%%%%%%%%%
    %%% CLASSIFICATION %%%
    %%%%%%%%%%%%%%%%%%%%%%

    %% 1. category classification [colours vs. scenes; coarse]

    %--- 1a encoding

    trlinfo = cell2mat(data_stimuli.trialinfo);

    sel1    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'color');
    sel2    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'scene');

    dataToClassifyTraining      = cat(1,data_stimuli.trial(sel1,:,:),data_stimuli.trial(sel2,:,:));
    clabelTraining              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
    samples_train               = nearest(data_stimuli.time,settings.TOI_train(1)):settings.timesteps_train*fsample:nearest(data_stimuli.time,settings.TOI_train(2));

    dataToClassifyTraining      = dataToClassifyTraining(:,:,samples_train);

    %--- 1b retrieval (separately for correct and incorrect source memory)

    trlinfo = cell2mat(data_ripples.trialinfo);

    for icond = 1:2

        if icond == 1

            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Memory],'SourceCorrect');
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & ismember([trlinfo.Memory],'SourceCorrect') ;

            perf{isubject}.correct.trl_num_test = [sum(sel1) sum(sel2)];

        elseif icond == 2

            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color')& (ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene')& (ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));

            perf{isubject}.incorrect.trl_num_test         = [sum(sel1) sum(sel2)];
        end

        if (sum(sel1) < 2 | sum(sel2) < 2)
            continue;
        end

        dataToClassifyTest      = cat(1,data_ripples.trial(sel1,:,:),data_ripples.trial(sel2,:,:));
        clabelTest              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
        samples_test            = nearest(data_ripples.time,settings.TOI_test(1)):settings.timesteps_test*fsample:nearest(data_ripples.time,settings.TOI_test(2));

        dataToClassifyTest      = dataToClassifyTest(:,:,samples_test);

        %--- [optional] Preprocess parameters (z-scoring)

        if settings.zscore_data4class
            dataToClassifyTraining  = zscore(dataToClassifyTraining);
            dataToClassifyTest      = zscore(dataToClassifyTest);
        end

        %--- Running classifier

        cfg                         = [];
        cfg.classifier              = settings.classifier;
        cfg.metric                  = settings.metric;
        [accuracy, ~]               = mv_classify_timextime(cfg, dataToClassifyTraining, clabelTraining, dataToClassifyTest, clabelTest);

        if icond == 1
            perf{isubject}.correct.accuracy   = accuracy;
        elseif icond == 2
            perf{isubject}.incorrect.accuracy = accuracy;
        end

        % add info
        perf{isubject}.channelcount_train  = size(data_stimuli.trial,2);
        perf{isubject}.channelcount_test   = size(data_ripples.trial,2);
        perf{isubject}.channels_test       = data_ripples.label;
        perf{isubject}.trl_num_train       = size(dataToClassifyTraining,1);
        perf{isubject}.time_train          = data_stimuli.time(samples_train);
        perf{isubject}.time_test           = data_ripples.time(samples_test);

    end

    %% 2. exemplar classification (fine)

    %% 2.1 colours [blue vs. red]

    %--- 2.1a encoding

    trlinfo = cell2mat(data_stimuli.trialinfo);

    sel1    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'blue');
    sel2    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'red');

    dataToClassifyTraining      = cat(1,data_stimuli.trial(sel1,:,:),data_stimuli.trial(sel2,:,:));
    clabelTraining              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
    samples_train               = nearest(data_stimuli.time,settings.TOI_train(1)):settings.timesteps_train*fsample:nearest(data_stimuli.time,settings.TOI_train(2));

    dataToClassifyTraining = dataToClassifyTraining(:,:,samples_train);

    %--- 2.1b retrieval (separately for correct and incorrect source memory)

    trlinfo = cell2mat(data_ripples.trialinfo);

    for icond = 1:2

        if icond == 1

            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'blue') & ismember([trlinfo.Memory],'SourceCorrect');
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'red')  & ismember([trlinfo.Memory],'SourceCorrect');

            perf{isubject}.colour.correct.trl_num_test    = [sum(sel1) sum(sel2)];

        elseif icond == 2

            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'blue') & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'red')  & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));

            perf{isubject}.colour.incorrect.trl_num_test  = [sum(sel1) sum(sel2)];
        end

        if (sum(sel1) < 2 | sum(sel2) < 2)
            continue;
        end

        dataToClassifyTest      = cat(1,data_ripples.trial(sel1,:,:),data_ripples.trial(sel2,:,:));
        clabelTest              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
        samples_test            = nearest(data_ripples.time,settings.TOI_test(1)):settings.timesteps_test*fsample:nearest(data_ripples.time,settings.TOI_test(2));

        dataToClassifyTest = dataToClassifyTest(:,:,samples_test);

        %--- [optional] Preprocess parameters (z-scoring)

        if settings.zscore_data4class
            dataToClassifyTraining  = zscore(dataToClassifyTraining);
            dataToClassifyTest      = zscore(dataToClassifyTest);
        end

        %--- Running classifier

        cfg                         = [];
        cfg.classifier              = settings.classifier;
        cfg.metric                  = settings.metric;
        [accuracy, ~]               = mv_classify_timextime(cfg, dataToClassifyTraining, clabelTraining, dataToClassifyTest, clabelTest);

        if icond == 1
            perf{isubject}.colour.correct.accuracy         = accuracy;
            perf{isubject}.colour.correct.trl_num_train    = [numel(find(clabelTraining==1)),numel(find(clabelTraining==2))];
        elseif icond == 2
            perf{isubject}.colour.incorrect.accuracy       = accuracy;
            perf{isubject}.colour.incorrect.trl_num_train  = [numel(find(clabelTraining==1)),numel(find(clabelTraining==2))];
        end

    end

    %% 2.2 scenes [indoor vs. outdoor]

    %--- 2.2a encoding

    trlinfo = cell2mat(data_stimuli.trialinfo);

    sel1    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'indoor','office'}));
    sel2    = ismember([trlinfo.ExpPhase],'encoding') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'outdoor','nature'}));

    dataToClassifyTraining      = cat(1,data_stimuli.trial(sel1,:,:),data_stimuli.trial(sel2,:,:));
    clabelTraining              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
    samples_train               = nearest(data_stimuli.time,settings.TOI_train(1)):settings.timesteps_train*fsample:nearest(data_stimuli.time,settings.TOI_train(2));
    dataToClassifyTraining      = dataToClassifyTraining(:,:,samples_train);


    %--- 2.2b retrieval (separately for correct and incorrect source memory)

    trlinfo     = cell2mat(data_ripples.trialinfo);

    for icond = 1:2

        if icond == 1

            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'indoor','office'}))  & ismember([trlinfo.Memory],'SourceCorrect');
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'outdoor','nature'})) & ismember([trlinfo.Memory],'SourceCorrect');

            perf{isubject}.scene.correct.trl_num_test     = [sum(sel1) sum(sel2)];

        elseif icond == 2

            sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'indoor','office'}))  & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));
            sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'outdoor','nature'})) & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));

            perf{isubject}.scene.incorrect.trl_num_test   = [sum(sel1) sum(sel2)];

        end

        if (sum(sel1) < 2 | sum(sel2) < 2)
            continue;
        end

        dataToClassifyTest      = cat(1,data_ripples.trial(sel1,:,:),data_ripples.trial(sel2,:,:));
        clabelTest              = cat(1,1*ones(sum(sel1),1),2*ones(sum(sel2),1));
        samples_test            = nearest(data_ripples.time,settings.TOI_test(1)):settings.timesteps_test*fsample:nearest(data_ripples.time,settings.TOI_test(2));
        dataToClassifyTest      = dataToClassifyTest(:,:,samples_test);

        %--- [optional] Preprocess parameters (z-scoring)
        if settings.zscore_data4class
            dataToClassifyTraining  = zscore(dataToClassifyTraining);
            dataToClassifyTest      = zscore(dataToClassifyTest);
        end

        %--- Running classifier

        cfg                         = [];
        cfg.classifier              = settings.classifier;
        cfg.metric                  = settings.metric;
        [accuracy, ~]               = mv_classify_timextime(cfg, dataToClassifyTraining, clabelTraining, dataToClassifyTest, clabelTest);

        if icond == 1

            perf{isubject}.scene.correct.accuracy         = accuracy;
            perf{isubject}.scene.correct.trl_num_train    = [numel(find(clabelTraining==1)),numel(find(clabelTraining==2))];

        elseif icond == 2

            perf{isubject}.scene.incorrect.accuracy       = accuracy;
            perf{isubject}.scene.incorrect.trl_num_train  = [numel(find(clabelTraining==1)),numel(find(clabelTraining==2))];

        end

    end

%         clearvars -except perf settings isubject subjects SubjectIDs
end
toc

delete(gcp);

return
%% Stats and plots

correct                 = cell(1);
incorrect              = cell(1);
correct_colour         = cell(1);
incorrect_colour        = cell(1);
correct_scene            = cell(1);
incorrect_scene         = cell(1);
baseline_all            = cell(1);
correct_scene_colour     = cell(1);
incorrect_scene_colour   = cell(1);

correct_dec_all     = [];
incorrect_dec_all   = [];

sigma = .5; % STD of 2D gaussian smoothing

count_trls = 1;
for isubject = 1:numel(subjects)

    if ~isfield(perf{isubject}, 'correct') | ~isfield(perf{isubject}, 'incorrect') |...
            ~isfield(perf{isubject}.correct,'accuracy') | ~isfield(perf{isubject}.incorrect,'accuracy') |...
            ~isfield(perf{isubject}.colour.correct,'accuracy') | ~isfield(perf{isubject}.colour.incorrect,'accuracy') |...
            ~isfield(perf{isubject}.scene.correct,'accuracy') | ~isfield(perf{isubject}.scene.incorrect,'accuracy')
        continue;
    end

    nu_trials(count_trls,:) = sum(perf{isubject}.correct.trl_num_test);

    % correct
    correct{1,count_trls}                                 = struct;
    correct{1,count_trls}.label                           = {'chan'};
    correct{1,count_trls}.dimord                          = 'chan_freq_time';
    correct{1,count_trls}.freq                            = perf{1,1}.time_train;
    correct{1,count_trls}.time                            = perf{1,1}.time_test;
    correct{1,count_trls}.powspctrm(1,:,:)                = perf{isubject}.correct.accuracy;
   

    correct_dec_all(count_trls,:,:) = imgaussfilt(perf{isubject}.correct.accuracy, sigma);

    % incorrect
    incorrect{1,count_trls}                               = correct{1,count_trls};
    incorrect{1,count_trls}.powspctrm(1,:,:)          = perf{isubject}.incorrect.accuracy;
   
    incorrect_dec_all(count_trls,:,:) = imgaussfilt(perf{isubject}.incorrect.accuracy, sigma);

    % correct colour
    correct_colour{1,count_trls}                          = correct{1,count_trls};
    correct_colour{1,count_trls}.powspctrm(1,:,:)     = perf{isubject}.colour.correct.accuracy;
  

    % incorrect colour
    incorrect_colour{1,count_trls}                        = correct{1,count_trls};
    incorrect_colour{1,count_trls}.powspctrm(1,:,:)   = perf{isubject}.colour.incorrect.accuracy;
   
    % correct scene
    correct_scene{1,count_trls}                           = correct{1,count_trls};
    correct_scene{1,count_trls}.powspctrm(1,:,:)     = perf{isubject}.scene.correct.accuracy;
   
    % incorrect scene
    incorrect_scene{1,count_trls}                        = correct{1,count_trls};
    incorrect_scene{1,count_trls}.powspctrm(1,:,:)   = perf{isubject}.scene.incorrect.accuracy;
   
    % correct scene and colour collapsed
    correct_scene_colour{1,count_trls}                    = correct{1,count_trls};
    correct_scene_colour{1,count_trls}.powspctrm(1,:,:) = mean(cat(1,correct_colour{1,count_trls}.powspctrm,correct_scene{1,count_trls}.powspctrm));
   
    correct_dec_fine_all(count_trls,:,:) = imgaussfilt(mean(cat(1,correct_colour{1,count_trls}.powspctrm,correct_scene{1,count_trls}.powspctrm)), sigma);

    % incorrect scene and colour collapsed
    incorrect_scene_colour{1,count_trls}                  = correct{1,count_trls};
    incorrect_scene_colour{1,count_trls}.powspctrm(1,:,:) = mean(cat(1,incorrect_colour{1,count_trls}.powspctrm,incorrect_scene{1,count_trls}.powspctrm));
  
    incorrect_dec_fine_all(count_trls,:,:) = imgaussfilt(mean(cat(1,incorrect_colour{1,count_trls}.powspctrm,incorrect_scene{1,count_trls}.powspctrm)), sigma);


    % baseline
    baseline_all{1,count_trls}                            = correct{1,count_trls};
    baseline_all{1,count_trls}.powspctrm(1,:,:)           = .5*ones(size(perf{isubject}.correct.accuracy));

    count_trls = count_trls+1;
end

%% Run stats scene and colour collapsed (fine-grained)

timeaxis = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);
timeaxis = correct{1, 1}.time(timeaxis);



contrast = 'incorrect'; % baseline incorrect

cfg                     = [];
cfg.latency             = [timeaxis(1) timeaxis(end)]; % ripple time
cfg.frequency           = [-.2 freqaxis(end)]; % encoding time  [freqaxis(1) freqaxis(end)]
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

cfg.numrandomization    = 500;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = numel(correct);
% set up design matrix

% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];
cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
if strcmp(contrast,'baseline')
    [Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, baseline_all{:});
elseif strcmp(contrast,'incorrect')
    [Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, incorrect_scene_colour{:});
end
length(find(Fieldtripstats.mask==1))

%% plot (significant) t vals

stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);

figure;
imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    tvals);
colormap(jet)
caxis([-5 5])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',12)
xlabel('ripple time (sec)')
ylabel('encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title(sprintf('fine-grained classifier, correct vs. %s\n%s correction, %s',contrast,cfg.correctm,cfg.method));

if any(Fieldtripstats.mask(:))
    plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats.mask)))
end

%% plot average, with significant pixels highlighted, followed by time course

this_TOI_enc = [0.5 3];

stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

this_TOI_idx = nearest(freqaxis,this_TOI_enc(1)):nearest(freqaxis,this_TOI_enc(2));


d = [];
for isubject = 1:numel(correct)
    d(:,:,isubject) = squeeze(correct_scene_colour{isubject}.powspctrm)-squeeze(incorrect_scene_colour{isubject}.powspctrm);
end

m = squeeze(nanmean(nanmean(d(this_TOI_idx,stats_time,:),1),3));
s = squeeze(nanstd(nanmean(d(this_TOI_idx,stats_time,:),1),0,3))./sqrt(size(d,3));


figure;
imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    nanmean(d,3));
colormap(jet)
caxis([0 .1])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',12)
xlabel('ripple time (sec)')
ylabel('encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title(sprintf('coarse-grained classifier, correct vs. %s\n%s correction, %s',contrast,cfg.correctm,cfg.method));

if any(Fieldtripstats.mask(:))
    pause(.5)
    plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats.mask)))
end

figure
boundedline(timeaxis,m,s,'k')
vline(0)
hline(0)
set(gca,'FontSize',12)
xlabel('ripple time (sec)')
ylabel('reinstatement')
set(gca,'TickDir','out')
title(sprintf('encoding averaged from %1.1f-%1.1fs',this_TOI_enc))

%% visualise included channels

elec_size   = 15;
transp      = 0.25;
extracolor  = .5;

epos_all    = [];

for isubject = 1:numel(subjects)
    epos_all = [epos_all;cell2mat({channs{isubject}.coords}')];
end

colvec = ones(1,size(epos_all,1));

views = [-90 0;0 0;180 -90];

for iview=1:size(views,1)
    figure
    plot_ecog(colvec, ...
        fullfile(ftpath,'template/anatomy/'),...
        epos_all,[-max(colvec) max(colvec)+extracolor], transp, views(iview,:), elec_size);
    colorbar off
end

%%