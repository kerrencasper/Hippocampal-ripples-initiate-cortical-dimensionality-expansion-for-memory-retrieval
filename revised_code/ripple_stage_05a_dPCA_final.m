

%                  Casper Kerren      [kerren@cbs.mpg.de]

clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

% [~,ftpath]=ft_version;

%% path settings
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';

% FC has no trials in blue correct
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');



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
addpath(genpath('/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/dPCA-master'))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/help_functions'))
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))


load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;

settings.colour_scheme_1 = brewermap(30,'RdBu');
settings.colour_scheme_1 = settings.colour_scheme_1;
%% pre-decoding settings
% ripple extraction

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1; % remove co-occuring ripples
settings.time_to_excl_RT            = .25; % exclude last 250 ms of trials, to make sure ripple event was in trial
settings.solo_ripple                = 1; % pick one (maxEnv) ripple per trial if multiple ripple events are found
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur
settings.do_surrogates              = 0; % switch time of ripples between trials.

%% decoding settings

% data preprocessing
settings.ori_fsample            = 1000; % original sample frequency
settings.do_resample            = 100; % [] or sample frequency

% baseline and zscoring settings
settings.zscore                 = 1;
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

%% start for loop

timeaxis = settings.TOI_test(1):settings.timesteps_test:settings.TOI_test(2);
freqaxis = settings.TOI_train(1):settings.timesteps_train:settings.TOI_train(2);

% the maximum number of trials for one condition is 46, so I will keep all
% trials in the structure, meaning a matrix will be (trials, channels,
% time). When concatenating them over channels (2nd dim), I need to have made the matrix
% nan(46,channels, time).

data_all_blue_correct       = [];
data_all_red_correct        = [];
data_all_indoor_correct     = [];
data_all_outdoor_correct    = [];

data_all_blue_incorrect     = [];
data_all_red_incorrect      = [];
data_all_indoor_incorrect   = [];
data_all_outdoor_incorrect  = [];


trialinfo_all     = {};
max_all = [];


numWorkers = 4; %

parpool('local', numWorkers);

tic
for isubject = 1:numel(subjects)

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


    if settings.zscore == 1
        data_ripples.trial  = zscore(data_ripples.trial,[],1);
    end

    % trialNum: N x S x D
    % firingRates: N x S x D x T x maxTrialNum
    % firingRatesAverage: N x S x D x T
    %
    % N is the number of neurons (in my case channels)
    % S is the number of stimuli conditions (F1 frequencies in Romo's task) (in
    % my case indoor outdoor red blue)
    % D is the number of decisions (D=2) (in my case correct incorrect)
    % T is the number of time-points (note that all the trials should have the
    % same length in time!)
    % E is the maximum number of trials

    trlinfo = cell2mat(data_ripples.trialinfo);

    TOI = nearest(data_ripples.time, -1):nearest(data_ripples.time, 1);

    sel1 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'blue') & ismember([trlinfo.Memory],'SourceCorrect');
    sel2 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'red')  & ismember([trlinfo.Memory],'SourceCorrect');

    sel3 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'indoor','office'}))  & ismember([trlinfo.Memory],'SourceCorrect');
    sel4 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'outdoor','nature'})) & ismember([trlinfo.Memory],'SourceCorrect');

    sel5 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'blue') & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));
    sel6 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'color') & ismember([trlinfo.Subcat],'red')  & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));

    sel7 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'indoor','office'}))  & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));
    sel8 = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.BlockType],'scene') & (ismember([trlinfo.Subcat],{'outdoor','nature'})) & (ismember([trlinfo.Memory],{'SourceIncorrect','SourceDunno','dunno'}));



    max_all(isubject,:) = max(sum([sel1;sel2;sel3;sel4;sel5;sel6;sel7;sel8],2));
    min_all(isubject,:) = min(sum([sel1;sel2;sel3;sel4;sel5;sel6;sel7;sel8],2));


    correct_trials      = sel1|sel2|sel3|sel4;
    incorrect_trials    = sel5|sel6|sel7|sel8;



    % create a matrix so that I have the maximum number of trials (46) in the
    % first. The rest will be nans.
    dataToClassifyBlueCorrect       = nan(46, length(data_ripples.label), length(data_ripples.time));
    dataToClassifyRedCorrect        = nan(46, length(data_ripples.label), length(data_ripples.time));
    dataToClassifyIndoorCorrect     = nan(46, length(data_ripples.label), length(data_ripples.time));
    dataToClassifyOutdoorCorrect    = nan(46, length(data_ripples.label), length(data_ripples.time));

    dataToClassifyBlueIncorrect     = nan(46, length(data_ripples.label), length(data_ripples.time));
    dataToClassifyRedIncorrect      = nan(46, length(data_ripples.label), length(data_ripples.time));
    dataToClassifyIndoorIncorrect   = nan(46, length(data_ripples.label), length(data_ripples.time));
    dataToClassifyOutdoorIncorrect  = nan(46, length(data_ripples.label), length(data_ripples.time));

    dataToClassifyBlueCorrect(1:size(data_ripples.trial(sel1,:,:),1),:,:)       = data_ripples.trial(sel1,:,:);
    dataToClassifyRedCorrect(1:size(data_ripples.trial(sel2,:,:),1),:,:)        = data_ripples.trial(sel2,:,:);
    dataToClassifyIndoorCorrect(1:size(data_ripples.trial(sel3,:,:),1),:,:)     = data_ripples.trial(sel3,:,:);
    dataToClassifyOutdoorCorrect(1:size(data_ripples.trial(sel4,:,:),1),:,:)    = data_ripples.trial(sel4,:,:);

    dataToClassifyBlueIncorrect(1:size(data_ripples.trial(sel5,:,:),1),:,:)     = data_ripples.trial(sel5,:,:);
    dataToClassifyRedIncorrect(1:size(data_ripples.trial(sel6,:,:),1),:,:)      = data_ripples.trial(sel6,:,:);
    dataToClassifyIndoorIncorrect(1:size(data_ripples.trial(sel7,:,:),1),:,:)   = data_ripples.trial(sel7,:,:);
    dataToClassifyOutdoorIncorrect(1:size(data_ripples.trial(sel8,:,:),1),:,:)  = data_ripples.trial(sel8,:,:);

    data_all_blue_correct       = cat(2,data_all_blue_correct,dataToClassifyBlueCorrect);
    data_all_red_correct        = cat(2,data_all_red_correct,dataToClassifyRedCorrect);
    data_all_indoor_correct     = cat(2,data_all_indoor_correct,dataToClassifyIndoorCorrect);
    data_all_outdoor_correct    = cat(2,data_all_outdoor_correct,dataToClassifyOutdoorCorrect);

    data_all_blue_incorrect     = cat(2,data_all_blue_incorrect,dataToClassifyBlueIncorrect);
    data_all_red_incorrect      = cat(2,data_all_red_incorrect,dataToClassifyRedIncorrect);
    data_all_indoor_incorrect   = cat(2,data_all_indoor_incorrect,dataToClassifyIndoorIncorrect);
    data_all_outdoor_incorrect  = cat(2,data_all_outdoor_incorrect,dataToClassifyOutdoorIncorrect);

    trialinfo_all{isubject} = trlinfo;

    clearvars -except data_all_blue_correct data_all_red_correct data_all_indoor_correct...
        data_all_outdoor_correct data_all_blue_incorrect data_all_red_incorrect...
        data_all_indoor_incorrect data_all_outdoor_incorrect isubject settings subjects SubjectIDs trialinfo_all max_all min_all

end

delete(gcp);

return

%%

% maximum number of trials is 46. The rest I have set to NanN.

fs          = 100;
timeEvents  = -1.2:1/fs:1.2;
TOI         = nearest(timeEvents,-1):nearest(timeEvents,1);
timeEvents  = timeEvents(TOI);
time        = 0;


dataToClassifyTest_correct      = cat(1, data_all_blue_correct, data_all_red_correct, data_all_indoor_correct,data_all_outdoor_correct);
dataToClassifyTest_incorrect    = cat(1, data_all_blue_incorrect, data_all_red_incorrect, data_all_indoor_incorrect,data_all_outdoor_incorrect);

dataToClassifyTest_correct      = reshape(dataToClassifyTest_correct, 1, size(dataToClassifyTest_correct,1), size(dataToClassifyTest_correct,2), size(dataToClassifyTest_correct,3));
dataToClassifyTest_incorrect    = reshape(dataToClassifyTest_incorrect, 1, size(dataToClassifyTest_incorrect,1), size(dataToClassifyTest_incorrect,2), size(dataToClassifyTest_incorrect,3));

dataToClassifyTest = cat(1, dataToClassifyTest_correct,dataToClassifyTest_incorrect);

dataToClassifyTest = dataToClassifyTest(:,:,:,TOI);

save dataToClassifyTest dataToClassifyTest

% use plot_dPCA to plot the results.

