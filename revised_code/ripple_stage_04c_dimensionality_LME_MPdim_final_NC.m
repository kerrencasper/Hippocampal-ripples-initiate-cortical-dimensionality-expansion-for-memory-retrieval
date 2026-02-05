%%
% [ripple_stage_04d_dimensionality_LME] find ripples in hippocampal channels,
% extract and realign data based on ripple events.
% Use encoding data to train LDA and test at ripple aligned
% data. Do PCA to estimate dimensionality of correct and incorrect,
% perform linear-mixed effects model.
%                  Bernhard Staresina [bernhard.staresina@psy.ox.ac.uk]
%                  Casper Kerren      [kerren@cbs.mpg.de]

clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

% [~,ftpath]=ft_version;

%% path settings
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
addpath('/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/scripts_to_publish')

load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;

settings.data_dir           = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];

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

addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/MVPA-Light-master'))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/help_functions'))
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))


settings.colour_scheme_1 = brewermap(30,'RdBu');

%% pre-decoding settings
% ripple extraction

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1; % remove co-occuring ripples
settings.time_to_excl_RT            = .25; % exclude last 250 ms of trials, to make sure ripple event was in trial
settings.solo_ripple                = 1; % pick one (maxEnv) ripple per trial if multiple ripple events are found
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur, eg [.5 1.5]
settings.do_surrogates              = 0; % switch time of ripples between trials.

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
settings.prc_overlap            = .9; % perceptage overlap sliding window

%% start for loop

timeaxis = settings.TOI_test(1):settings.timesteps_test:settings.TOI_test(2);
freqaxis = settings.TOI_train(1):settings.timesteps_train:settings.TOI_train(2);

perf    = cell(1,numel(subjects));
channs  = cell(1,numel(subjects));

RT_all_subj_correct     = cell(1,numel(subjects));
RT_all_subj_incorrect   = cell(1,numel(subjects));
ripple_time = cell(1,numel(subjects));

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

    %% [optional] do surrogates by taking time of ripple from other trial
    if settings.do_surrogates

        randtrials      = randperm(size(trl_ripple,1));
        tmp             = trl_ripple(:,2)-round(trl_ripple(:,3).*data.fsample); % find cue onset
        tmp             = tmp+round(trl_ripple(randtrials,3).*data.fsample); % add another ripple's event time
        trl_ripple(:,2) = tmp;
        trl_ripple(:,3) = trl_ripple(randtrials,3);

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
    ripple_time{isubject} = ripple_tmp;

    %% Do simensionality for each block instead, using a linear mixed effect modelling
    trlinfo     = cell2mat(data_ripples.trialinfo);
    these_runs  = [trlinfo.RunNumber];
    avail_runs  = unique(these_runs);


    accuracy_correct    = [];
    accuracy_incorrect  = [];
    for irun = 1:numel(avail_runs)

        trlinfo     = cell2mat(data_ripples.trialinfo);

        idx         = these_runs == avail_runs(irun);
        data_tmp    = data_ripples.trial(idx,:,:);

        trlinfo = trlinfo(idx,:);


        %%

        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% PCA %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%


        %% PCA on the data to get the eigenvalues that explain 90% of the variance.
        % Do it with a sliding window of 60ms with 100Hz)


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

            dataToClassifyTest      = cat(1,data_tmp(sel1,:,:),data_tmp(sel2,:,:));
            if size(dataToClassifyTest,1) < 2
                continue
            end

            fs              = 1/(data_ripples.time(2)-data_ripples.time(1));
            chunk_size      = round(settings.nu_time_points/((1/fs)*1000)); % Size of each chunk
            overlap_size    = round(chunk_size * settings.prc_overlap);
            num_chunks      = floor((size(dataToClassifyTest, 3) - overlap_size) / (chunk_size - overlap_size));
            TOI_ripple      = linspace(data_ripples.time(1), data_ripples.time(end),num_chunks);

            explained_variances     = zeros(num_chunks,1);
            for i = 1:num_chunks
                % Calculate the start and end indices of the current chunk
                start_idx   = (i - 1) * (chunk_size - overlap_size) + 1;
                end_idx     = start_idx + chunk_size - 1;

                % Extract data for the current chunk and reshape
                chunk_data = dataToClassifyTest(:, :, start_idx:end_idx);
                chunk_data = reshape(chunk_data, size(dataToClassifyTest, 1), []);

                if settings.smooth_before_dim == 1
                    for ch = 1:size(chunk_data, 1)
                        valid_idx = ~isnan(chunk_data(ch,:));
                        if sum(valid_idx) < 2
                            continue
                        end
                        chunk_data(ch,:) = interp1( find(valid_idx), ...
                            chunk_data(ch,valid_idx), ...
                            1:size(chunk_data,2), ...
                            'linear','extrap');
                    end
                end

                % optionally drop completely NaN channels
                good_ch = ~all(isnan(chunk_data),2);
                chunk_data = chunk_data(good_ch,:);

                % 2) z-score each channel across time (mean 0, var 1)
                chunk_z = zscore(chunk_data, 0, 2);    % channels x time

                [nChan, nTime] = size(chunk_z);

                % 3) channel × channel correlation matrix
                %    (z-scored → covariance == correlation)
                C = (chunk_z * chunk_z.') / nTime;     % nChan x nChan

                % 4) eigenvalues of C
                eigvals = eig(C);                      % unsorted
                eigvals = sort(eigvals, 'descend');

                % 5) Marčenko–Pastur upper bound
                q = nTime / nChan;                     % B/n  (should be ≥ 1)
                if q < 1
                    q = 1/q;                           % MP formula assumes q ≥ 1
                end
                lambda_max = (1 + sqrt(1/q))^2;        % σ^2 = 1 after z-scoring

                % 6) number of "signal" components
                nSignal = sum(eigvals > lambda_max);

                % You can also store the spectrum and λ_max if you want to plot them
                %     eig_spectra{i}  = eigvals;
                %     lambda_max_all(i) = lambda_max;

                % This is your dimensionality estimate for this chunk
                explained_variances(i) = nSignal;

            end

            if icond == 1
                accuracy_correct(irun,:) = explained_variances;

            elseif icond == 2
                accuracy_incorrect(irun,:) = explained_variances;

            end

            % add info

        end

    end
    perf{isubject}.channelcount_test   = size(data_ripples.trial,2);
    perf{isubject}.channels_test       = data_ripples.label;
    perf{isubject}.time_train          = TOI_ripple;
    perf{isubject}.time_test           = TOI_ripple;

    perf{isubject}.correct.accuracy     = accuracy_correct;
    perf{isubject}.incorrect.accuracy   = accuracy_incorrect;


    %     clearvars -except perf settings isubject subjects SubjectIDs RT_all_subj_correct RT_all_subj_incorrect
end
toc

delete(gcp);

return
%% create table for LME

RT_correct = cellfun(@mean, RT_all_subj_correct);
RT_max_correct = cellfun(@max, RT_all_subj_correct);
RT_min_correct = cellfun(@min, RT_all_subj_correct);

RT_incorrect = cellfun(@mean, RT_all_subj_incorrect);
RT_max_incorrect = cellfun(@max, RT_all_subj_incorrect);
RT_min_incorrect = cellfun(@min, RT_all_subj_incorrect);


% first for accuracy and condition

subjectIDs = [];
blockIDs = [];
conditions = [];
accuracies = [];

TOI = nearest(perf{1, 1}.time_test,-1):nearest(perf{1, 1}.time_test,1);

for isubject = 1:numel(perf)
    perf_correct     = zscore(perf{isubject}.correct.accuracy(:,TOI),[],2);
    perf_incorrect   = zscore(perf{isubject}.incorrect.accuracy(:,TOI),[],2);

    idx = (perf_correct(:,1) ~= 0 & perf_incorrect(:,1) ~= 0);

    perf_correct    = perf_correct(idx,:);
    perf_incorrect  = perf_incorrect(idx,:);

    % Correct trials
    num_blocks_correct = size(perf_correct, 1);

    for iblock = 1:num_blocks_correct
        num_samples = size(perf_correct, 2);
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples, 1)];
        conditions  = [conditions; repmat(1, num_samples, 1)];  % 1 for correct
        accuracies  = [accuracies; perf_correct(iblock, :)'];
    end

    % Incorrect trials
    num_blocks_incorrect = size(perf_incorrect, 1);

    for iblock = 1:num_blocks_incorrect
        num_samples = size(perf_incorrect, 2);
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples, 1)];
        conditions  = [conditions; repmat(2, num_samples, 1)];  % 2 for incorrect
        accuracies  = [accuracies; perf_incorrect(iblock, :)'];
    end
end


dataTable = table(subjectIDs, blockIDs, conditions, accuracies, ...
    'VariableNames', {'SubjectID', 'BlockID', 'Condition', 'Accuracy'});


% Define the formula for the linear mixed-effects model (SubjectID and BlockID as random effects)
lme = fitlme(dataTable, 'Accuracy ~ Condition + (1|SubjectID) + (1|BlockID)');

disp(lme)

%% now accuracy, conditions and before/after ripple

RT_correct = cellfun(@mean, RT_all_subj_correct);
RT_incorrect = cellfun(@mean, RT_all_subj_incorrect);


subjectIDs  = [];
blockIDs    = [];
conditions  = [];
halves      = [];
accuracies  = [];

for isubject = 1:numel(perf)
    perf_correct     = perf{isubject}.correct.accuracy(:,TOI);
    perf_incorrect   = perf{isubject}.incorrect.accuracy(:,TOI);

    idx = (perf_correct(:,1) ~= 0 & perf_incorrect(:,1) ~= 0);

    perf_correct    = perf_correct(idx,:);
    perf_incorrect  = perf_incorrect(idx,:);

    % Correct trials
    num_blocks_correct = size(perf_correct, 1);
    for iblock = 1:num_blocks_correct
        num_samples = size(perf_correct, 2);
        half_point  = floor(num_samples / 2);
        % First half
        subjectIDs  = [subjectIDs; repmat(isubject, half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, half_point, 1)];
        conditions  = [conditions; repmat(1, half_point, 1)];  % 1 for correct
        halves      = [halves; repmat(1, half_point, 1)];  % 1 for first half
        accuracies  = [accuracies; perf_correct(iblock, 1:half_point)'];
        % Second half
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples - half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples - half_point, 1)];
        conditions  = [conditions; repmat(1, num_samples - half_point, 1)];  % 1 for correct
        halves      = [halves; repmat(2, num_samples - half_point, 1)];  % 2 for second half
        accuracies  = [accuracies; perf_correct(iblock, half_point+1:end)'];
    end

    % Incorrect trials
    num_blocks_incorrect = size(perf_incorrect, 1);
    for iblock = 1:num_blocks_incorrect
        num_samples = size(perf_incorrect, 2);
        half_point  = floor(num_samples / 2);
        % First half
        subjectIDs  = [subjectIDs; repmat(isubject, half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, half_point, 1)];
        conditions  = [conditions; repmat(2, half_point, 1)];  % 2 for incorrect
        halves      = [halves; repmat(1, half_point, 1)];  % 1 for first half
        accuracies  = [accuracies; perf_incorrect(iblock, 1:half_point)'];
        % Second half
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples - half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples - half_point, 1)];
        conditions  = [conditions; repmat(2, num_samples - half_point, 1)];  % 2 for incorrect
        halves      = [halves; repmat(2, num_samples - half_point, 1)];  % 2 for second half
        accuracies  = [accuracies; perf_incorrect(iblock, half_point+1:end)'];
    end
end

% Create a table
dataTable = table(subjectIDs, blockIDs, conditions, halves, accuracies, ...
    'VariableNames', {'SubjectID', 'BlockID', 'Condition', 'Half', 'Accuracy'});


% Include interaction between Condition and Half
lme = fitlme(dataTable, 'Accuracy ~ Condition * Half + (1|SubjectID) + (1|BlockID)');

% Display the results
disp(lme)

%% PLOT
% Extract unique subject IDs
uniqueSubjects = unique(dataTable.SubjectID);

% Initialize arrays to store the means and SEMs
means = zeros(numel(uniqueSubjects), 2, 2);  % Dimensions: subject x condition x half
sems = zeros(numel(uniqueSubjects), 2, 2);

% Loop through each subject, condition, and half to calculate means and SEMs
for isubject = 1:numel(uniqueSubjects)
    for condition = 1:2
        for half = 1:2
            % Filter data for the current subject, condition, and half
            subset = dataTable(dataTable.SubjectID == uniqueSubjects(isubject) & ...
                dataTable.Condition == condition & ...
                dataTable.Half == half, :);
            % Calculate the mean and SEM
            means(isubject, condition, half) = mean(subset.Accuracy);
            sems(isubject, condition, half) = std(subset.Accuracy) / sqrt(height(subset));
        end
    end
end


% Create labels for the plot
conditionLabels = {'Correct', 'Incorrect'};
halfLabels = {'First Half', 'Second Half'};

% Plot the interaction effect for each participant

mean_correct    = squeeze(means(:,1,:));
mean_incorrect  = squeeze(means(:,2,:));

data_lme = {};
data_lme{1, 1} = mean_correct(:,2)-mean_correct(:,1);
data_lme{2, 1} = mean_incorrect(:,2)-mean_incorrect(:,1);


glm_fit = [];
yfit = [];

for isubject = 1:numel(perf)
    x = [data_lme{1, 1}(isubject,1),1;data_lme{2, 1}(isubject,1),2];
    [bb,dev,~] = glmfit(x(:,2),x(:,1),'normal');
    glm_fit(isubject,:) = bb(2,1);
    yfit(isubject,:) = polyval([bb(2,1),bb(1,1)],[1,2]);
end


h   = rm_raincloud(data_lme, [settings.colour_scheme_1(7,:)],0,'ks',[],settings.colour_scheme_1,yfit);
h.p{1, 1}.FaceColor = settings.colour_scheme_1(24,:);
h.p{2, 1}.FaceColor = settings.colour_scheme_1(3,:);




hold on
set(gca,'TickDir','out')
xlabel('dimensionality')
set(gca,'FontSize',24)
% yticks([.3 10 30 50 70 90 110 130])
yticklabels({'incorrect 2nd - incorrect 1st', 'correct 2nd - correct 1st'})
ylabel('condition')
axis tight
range=axis;


data_lme = {};
data_lme{1, 1} = mean_correct(:,1);
data_lme{2, 1} = mean_correct(:,2);
data_lme{3, 1} = mean_incorrect(:,1);
data_lme{4, 1} = mean_incorrect(:,2);


figure;
h   = rm_raincloud(data_lme, [settings.colour_scheme_1(7,:)],0,'ks',[],settings.colour_scheme_1);
h.p{1, 1}.FaceColor = settings.colour_scheme_1(24,:);
h.p{2, 1}.FaceColor = settings.colour_scheme_1(21,:);
h.p{3, 1}.FaceColor = settings.colour_scheme_1(7,:);
h.p{4, 1}.FaceColor = settings.colour_scheme_1(3,:);



hold on
set(gca,'TickDir','out')
xlabel('dimensionality')
set(gca,'FontSize',24)
% yticks([.3 10 30 50 70 90 110 130])
yticklabels({'incorrect second', 'incorrect first','correct second', 'correct first'})
ylabel('condition')
axis tight
range=axis;










data_to_plot = [mean_correct(:,2)-mean_correct(:,1),mean_incorrect(:,2)-mean_incorrect(:,1)];



N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector

scattercorrect = scatter(.9, mean_correct(:,2)-mean_correct(:,1), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;


scattercorrect = scatter(1.1, mean_incorrect(:,2)-mean_incorrect(:,1), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;





hold on;
[f, xi] = ksdensity(mean_correct(:,2)-mean_correct(:,1),'Bandwidth',.1); 
fill(.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2)

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(mean_incorrect(:,2)-mean_incorrect(:,1),'Bandwidth',.1); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;

    line([.9 1.1], [mean_correct(:,2)-mean_correct(:,1), mean_incorrect(:,2)-mean_incorrect(:,1)], 'Color', [0.5 0.5 0.5]);

    
end




[a,b,c,d] = ttest(mean_correct(:,2)-mean_correct(:,1),mean_incorrect(:,2)-mean_incorrect(:,1))

[a,b,c,d] = ttest(mean_correct(:,2)-mean_correct(:,1))

[a,b,c,d] = ttest(mean_incorrect(:,2)-mean_incorrect(:,1))




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