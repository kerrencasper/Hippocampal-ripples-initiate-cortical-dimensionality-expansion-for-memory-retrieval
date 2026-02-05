clear
close all


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

addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/MVPA-Light-master'))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/help_functions'))
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))

%% ripple settings

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1; % remove co-occuring ripples
settings.time_to_excl_RT            = .25; % exclude last 250 ms of trials, to make sure ripple event was in trial
settings.rippleselection            = 0; % use all (0), short (1) or long (2) ripples
settings.ripples_most_chan          = 0; % use only ripples from the channel with greatest number of ripples
settings.solo_ripple                = 1; % pick one ripple per trial if multiple ripple events are found
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur, eg [.5 1.5]

%% TFR settings

% define peak spectral power

settings.FOI_theta              = .5:10; % for good resolution
settings.cycles_theta           = ceil(settings.FOI_theta * 0.5); % ~ 500 ms for each frequency
% settings.cycles_theta(settings.cycles_theta < 5)    = 5;

settings.FOI_gamma              = 30:5:150;
settings.cycles_gamma           = ceil(settings.FOI_gamma * 0.5); % ~ 500 ms for each frequency
settings.cycles_gamma(settings.cycles_gamma < 5)    = 5;
settings.peak_frequency         = 1; % 1 == PAC based on theta and gamma peak, 2 == for given frequency range
settings.peak_method            = 2; % 1 == IRASA, 2 == FOOOF

% for TF after peak estimate

settings.bs_correct             = 0; % 0 = none | 1 = % baseline change | 2 = zscore each trial
settings.bs_period              = [-.5 -.1]; % for option 1 above

settings.TOI                    = [-1 1]; % time around ripple
settings.timesteps              = .020; % for TFR
settings.FOI                    = [1:1:10, 30:5:150];

settings.cycles                         = ceil(settings.FOI * 0.5); % ~ 500 ms for each frequency
settings.cycles(settings.cycles < 5)    = 5;

settings.TFR_padding = ceil(max((1./settings.FOI).*settings.cycles)/2); % to ensure sufficient time for spectral resolution is available

%% PAC settings

settings.PAC_or_MI       = 1; % 1 for PAC (mean vector length) (Bragin et al, 1995; Canolty et al 2006; Jensen and Colgin, 2007; Lakatos et al, 2005) and 2 for MI (Tort et al).
settings.PAC_phase       = [1:1:10]; % frequencies for phase data
settings.PAC_power       = 30:5:150; % frequencies for phase data
settings.PAC_TOI         = [-1 1]; % time of interest for MI calculation
settings.del_erp         = 0; % if 1, erp is subtracted from each trial before TF
settings.PAC_bins        = 18; % nu bins for MI
settings.PAC_binEdges    = linspace(-pi,pi,settings.PAC_bins+1); % -pi to pi
settings.PAC_binCenters  = settings.PAC_binEdges(1:end-1)-diff(settings.PAC_binEdges)/2;
settings.PAC_nusurro     = 250; % number of surrogates for baseline

%% anatomical channel selection for TFR analysis

settings.channel_sel              = 1; % 1 exclude hippo, 2 only hippo, 3 all channels
settings.channel_posteriors       = 0; % only include posterior temporal and occipital channels y < 0 and z < 25
settings.channel_AAL              = 0; % only contacts that fall within predefined AAL regions
settings.pick_one_chan_per_ROI    = 0; % pick one channel within each ROI from AAL.

%% start for loop

numWorkers = 6; %

parpool('local', numWorkers);

tic

PAC = cell(1,numel(subjects));

parfor isubject = 1:numel(subjects)

    fprintf('processing subject %01d/%02d\n',isubject,numel(subjects));

    %% LOAD
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

    perf{isubject}.nu_HC_contacts = numel(channels);

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

    %% Load data to realign based on ripples

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



    %% Load subject file and change RT

    [numbers,strings] = xlsread([settings.base_path_castle,'well01_behavior_all.xls']);

    strings = strings(2:end,:);

    if isnan(numbers(1,1))
        numbers = numbers(2:end,:);
    end

    sel         = find(strcmp(strings(:,2),SubjectID));
    sel         = sel(1:numel(onsets));

    trls_enc    = strcmp(strings(sel,4),'encoding');
    trls_ret    = strcmp(strings(sel,4),'retrieval');

    RT                          = numbers(sel,11);
    RT(RT==-1 & trls_enc==1)    = 3; % -1 no press in time - set to 3s at encoding
    RT(RT==-1 & trls_ret==1)    = 5; % -1 no press in time - set to 5s at encoding and 5s at retrieval

    if settings.full_enc_trial
        RT(trls_enc==1) = 3; % [optional] set all encoding to 3s
    end


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

    %% create onset matrices (remove last xxx ms to ensure ripple event in trial)

    onsetmat = [onsets; onsets+(RT'.*data.fsample)-(settings.time_to_excl_RT*data.fsample)]';


    %% Extract artifacts
    % These are the timing of the artifacts in the data.
    % This is important later on when I want to divide the number of
    % ripples found in one trial with the amount of data we had in that
    % trial to find a ripple. It is not possible to find a ripple during
    % A) and artifact, and B) when we pad around the artifact.
    % I here add the padding immidiately, so that I can later exclude it.
    % I also want to make sure that if there is less than 39ms between
    % end of one artifact and beginning of next, I set this to artifact too.

    disp('finding all artifacts...')

    datectable_per_trial = nan(size(onsetmat,1),numel(alldat));

    for ichannel = 1:numel(alldat)

        datectable = ones(1,numel(data.time{1}));

        % mark all artifact including padding period
        for iartifact = 1:size(data.artifact{ichannel},1)

            idx             = data.artifact{ichannel}(iartifact,1)+alldat{ichannel}.param.artfctPad(1)*data.fsample : data.artifact{ichannel}(iartifact,2)+alldat{ichannel}.param.artfctPad(2)*data.fsample;
            idx             = idx(idx>0 & idx <= numel(data.time{1}));
            datectable(idx) = 0;
        end

        % mark all clean segments shorter than minimum ripple duration as
        % unavailable
        transitions      = [datectable(1) diff(datectable)];
        clean_beginnings = find(transitions == 1);
        clean_endings    = find(transitions == -1);

        for iblock = 1:numel(clean_beginnings)-1
            if clean_endings(iblock)-clean_beginnings(iblock) < alldat{ichannel}.criterion.len(1) * data.fsample
                datectable(clean_beginnings(iblock):clean_endings(iblock)) = 0;
            end
        end

        % write out how many detectable seconds there are per trial and channel
        for itrial = 1:size(onsetmat,1)
            datectable_per_trial(itrial,ichannel) = sum(datectable(round(onsetmat(itrial,1)):round(onsetmat(itrial,2))))/data.fsample;
        end
        detectable = [];

    end


    %% extract relevant info

    ripples_trial_channel               = struct;
    ripples_trial_channel.number        = nan(size(onsetmat,1),numel(alldat));
    ripples_trial_channel.duration      = nan(size(onsetmat,1),numel(alldat));
    ripples_trial_channel.frequency     = nan(size(onsetmat,1),numel(alldat));
    ripples_trial_channel.amplitude     = nan(size(onsetmat,1),numel(alldat));
    ripples_trial_channel.activity      = nan(size(onsetmat,1),numel(alldat));
    ripples_trial_channel.latency_first = nan(size(onsetmat,1),numel(alldat));
    ripples_trial_channel.latency_all   = nan(size(onsetmat,1),numel(alldat));

    for ichannel = 1:numel(alldat)

        evs = alldat{ichannel}.evtIndiv.maxTime;

        for itrial = 1:size(onsetmat,1)

            these_ripples = evs >= onsetmat(itrial,1) & evs <= onsetmat(itrial,2);
            ripples_trial_channel.number(itrial,ichannel) = sum(these_ripples);

            if any(these_ripples)

                latencies = (evs(these_ripples) - onsetmat(itrial,1))/data.fsample;

                ripples_trial_channel.duration(itrial,ichannel)     = mean(alldat{ichannel}.evtIndiv.duration(these_ripples));
                ripples_trial_channel.frequency(itrial,ichannel)    = mean(alldat{ichannel}.evtIndiv.freq(these_ripples));
                ripples_trial_channel.amplitude(itrial,ichannel)    = mean(alldat{ichannel}.evtIndiv.maxAmp(these_ripples));
                ripples_trial_channel.activity(itrial,ichannel)     = mean(alldat{ichannel}.evtIndiv.envSum(these_ripples));

                ripples_trial_channel.latency_first(itrial,ichannel) = latencies(1);
                ripples_trial_channel.latency_all(itrial,ichannel)   = mean(latencies);

            end
        end
    end

    ripples_trial_channel.density = ripples_trial_channel.number ./ datectable_per_trial;



    %% channel selection for TFR data

    % 1 exclude hippo, 2 only hippo, 3 all channels
    tmp                             = load([settings.data_dir_channels,'/channels_to_exclude_all_hipp_both_hem.mat']);
    channels_to_exclude_all_hipp    = tmp.channels_to_exclude_all_hipp;
    tmp                             = load([settings.data_dir_channels,'/channels_hipp_ripples.mat']);
    channels_hipp_ripples           = tmp.channels_hipp_ripples;

    % for later tfr
    cfg         = [];
    cfg.channel = data.label;

    switch settings.channel_sel
        case 1
            cfg.channel = setdiff(setdiff([data.label],char(channels_to_exclude_all_hipp{isubject,:})),settings.scalp_channels);
        case 2 % ensure that hipp data are locked to ripple-providing channels only
            cfg.channel = intersect(cellstr(setdiff([data.label],settings.scalp_channels)),char(channels_hipp_ripples{isubject,:}));
            tmp         = cell2mat(alldat);
            cfg.channel = intersect(cfg.channel,[tmp.detectCh]);
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


    if settings.channel_posteriors == 1 % only include posterior temporal and occipital channels abs(x) > 40 | y < 0 | z < 25

        post_chan       = these_labels(these_y < 0 & these_z < 25);
        %         post_chan       = these_labels(these_y > -21);
        cfg.channel     = intersect(cfg.channel,post_chan);
    end

    if settings.channel_AAL == 1

        addpath(settings.AAL_dir);
        addpath(settings.SPM_dir);

        aal     = spm_read_vols(spm_vol(fullfile(settings.AAL_dir,'AAL3v1_1mm.nii')));
        origin  = [91 127 73]; % this is the [0 0 0] origin in voxel space

        tmp     = load(fullfile(settings.AAL_dir,'ROI_MNI_V7_List.mat'));
        names   = {tmp.ROI.Nom_L};
        tmp     = [];

        %         ROIs    = {'Temporal_Sup','Temporal_Mid','Temporal_Inf','Fusiform','ParaHippocampal','Lingual','Calcarine','Precuneus'};
        %         ROIs    = {'Temporal_Sup','Temporal_Mid','Temporal_Inf','Fusiform','ParaHippocampal'};
        ROIs    = {'Temporal_Sup','Temporal_Mid','Temporal_Inf','Fusiform'};
        %         ROIs    = {'Hippocampus'};

        if settings.pick_one_chan_per_ROI == 1  % pick one channel in each ROI

            counter     = 1;
            keep_these  = [];
            for iroi = 1:numel(ROIs)

                ROI_brain = ismember(aal,find(contains(names,ROIs(iroi))));

                sel = nan(size(cfg.channel,1),1);

                for ichannel = 1:size(cfg.channel,1)

                    idx             = strcmp(cfg.channel{ichannel},these_labels);
                    sel(ichannel)   = ROI_brain(origin(1)+these_x(idx),origin(2)+these_y(idx),origin(3)+these_z(idx));
                end
                if any(sel)

                    idx_chan            = find(sel);
                    keep_these(counter) = randsample(idx_chan,1);
                    counter = counter + 1;

                end

            end
            sel             = zeros(size(cfg.channel,1),1);
            sel(keep_these) = 1;

        else % take all channels that fall within ROIs

            ROI_brain = ismember(aal,find(contains(names,ROIs)));

            sel = nan(size(cfg.channel,1),1);

            for ichannel = 1:size(cfg.channel,1)

                idx              = strcmp(cfg.channel{ichannel},these_labels);
                sel(ichannel)    = ROI_brain(origin(1)+these_x(idx),origin(2)+these_y(idx),origin(3)+these_z(idx));

            end
        end

        cfg.channel = cfg.channel(sel==1);

    end

    data = ft_selectdata(cfg, data);

    % write out info of retained channels

    for ichannel = 1:numel(data.label)

        idx = strcmp(data.label{ichannel},these_labels);

        channs{isubject}(ichannel).names  = data.label{ichannel};
        channs{isubject}(ichannel).coords = [these_x(idx) these_y(idx) these_z(idx)];
    end

    perf{isubject}.nu_extra_HC_contacts = numel(cfg.channel);


    trl_ripple  = [];
    cnt         = 0;
    
    for ichannel = 1:numel(alldat)
        
        evs     = alldat{ichannel}.evtIndiv.maxTime;
        envSum  = alldat{ichannel}.evtIndiv.envSum;

        amplitude_ripples = alldat{ichannel}.evtIndiv.maxAmp;
        
        ripple_dur = alldat{ichannel}.evtIndiv.duration;


        
        if settings.rippleselection == 0 % all ripples
            duration_sel = logical(ones(1,numel(ripple_dur)));
        elseif settings.rippleselection == 1 % short ripples
            duration_sel = ripple_dur < median(ripple_dur);
        elseif settings.rippleselection == 2 % long ripples
            duration_sel = ripple_dur > median(ripple_dur);
        end
        
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
                trl_ripple(cnt,6) = ripples_trial_channel.number(this_event,ichannel); % number of ripples in this channel and trial
                trl_ripple(cnt,7) = amplitude_ripples(iripple); % amplitude of ripple
                trl_ripple(cnt,8) = ripple_dur(iripple); % duration of ripple
            end
        end
        
    end
    
    trl_ripple = sortrows(trl_ripple,1);
    
    %% Create trial structure around ripples
    
    pretrig      = round((abs(settings.TOI(1))+settings.TFR_padding) * data.fsample); % enough time to baseline correct later
    posttrig     = round((abs(settings.TOI(2))+settings.TFR_padding) * data.fsample);
    
    cfg          = [];
    cfg.trl      = [trl_ripple(:,2)-pretrig trl_ripple(:,2)+posttrig -pretrig*ones(size(trl_ripple,1),1)];
    
    data_ripples = ft_redefinetrial(cfg,data);
    % add trial info for each ripple trial, accounting for multiple ripples
    % per trial
    data_ripples.trialinfo = [];
    
    for itrial = 1:numel(data_ripples.trial)
        
        corresponding_trialinfo = find([trialinfo.EventNumber]==trl_ripple(itrial,1));
        
        trl_info                        = trialinfo(corresponding_trialinfo);
        trl_info.sample_ripple          = trl_ripple(itrial,2);
        trl_info.time_ripple            = trl_ripple(itrial,3);
        trl_info.channel                = trl_ripple(itrial,4);
        trl_info.name_channel           = {alldat{trl_ripple(itrial,4)}.evtIndiv.label};
        trl_info.envSum                 = trl_ripple(itrial,5);
        trl_info.duration               = trl_ripple(itrial,8);
        trl_info.amplitude              = trl_ripple(itrial,7);
        trl_info.nu_ripples             = trl_ripple(itrial,6);
        
        data_ripples.trialinfo{itrial,1}    = trl_info;
        
    end
    
    %% [optional] if there are multiple ripples per trial - pick the one ripple with highest activity (captured in the summed envelope metric)
    
    if settings.solo_ripple
        
        tmp_trl_info    = data_ripples.trialinfo;
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
        
        cfg                     = [];
        cfg.trials              = sel;
        data_ripples            = ft_selectdata(cfg, data_ripples);
    end
    
    %% [optional] pick ripples in a specific time window
    
    if any(settings.ripple_latency)

        trlinfo         = cell2mat(data_ripples.trialinfo);
        sel             = [trlinfo.time_ripple] > settings.ripple_latency(1) & [trlinfo.time_ripple] < settings.ripple_latency(2);
        trlinfo         = data_ripples.trialinfo;

        cfg                     = [];
        cfg.trials              = sel;
        data_ripples            = ft_selectdata(cfg, data_ripples);

    end

    trlinfo         = cell2mat(data_ripples.trialinfo);

    sel1 = (ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.Memory],'SourceCorrect'));
    sel2 = (ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));

    perf{isubject}.nu_correct_trials = sum(sel1);
    perf{isubject}.nu_incorrect_trials = sum(sel2);

end

delete(gcp);

return


for isubject = 1:size(settings.subjects,1)

hc_contacts(isubject,:) = perf{isubject}.nu_HC_contacts;  
ctx_contacts(isubject,:) = perf{isubject}.nu_extra_HC_contacts;  
corr_trials(isubject,:) = perf{isubject}.nu_correct_trials;  
incorr_trials(isubject,:) = perf{isubject}.nu_incorrect_trials;  

end