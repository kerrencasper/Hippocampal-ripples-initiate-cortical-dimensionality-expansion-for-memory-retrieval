%%
% [ripple_stage_05a_PAC] find ripples in hippocampal channels,
% extract and realign data based on ripple events.
% estimate phase-amplitude coupling between HC and cortex
%                  Bernhard Staresina [bernhard.staresina@psy.ox.ac.uk]
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear
close all


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

load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;

settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
subjects                    = {'CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK'};
SubjectIDs                  = {'01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK'};
settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};

addpath(genpath([paths.MVPA_Light_master]))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([paths.subfunctions])
addpath(genpath([paths.help_functions]))
addpath(genpath([paths.plotting]))

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

settings.FOI_theta              = 1:30; % for good resolution
settings.cycles_theta           = ceil(settings.FOI_theta * 0.5); % ~ 500 ms for each frequency
settings.cycles_theta(settings.cycles_theta < 5)    = 5;

settings.FOI_gamma              = 30:5:150;
settings.cycles_gamma           = ceil(settings.FOI_gamma * 0.5); % ~ 500 ms for each frequency
settings.cycles_gamma(settings.cycles_gamma < 5)    = 5;
settings.peak_frequency         = 1; % 1 == PAC based on theta and gamma peak, 2 == for given frequency range
settings.peak_method            = 2; % 1 == IRASA, 2 == FOOOF

% for TF after peak estimate

settings.bs_correct             = 1; % 0 = none | 1 = % baseline change | 2 = zscore each trial
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

    %% add back hippocampal data

    data_in     = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',subjects{isubject}],'data','onsets_session');
    data_hipp   = data_in.data;
    data_in     = [];

    data_hipp.sampleinfo = 1+data_hipp.sampleinfo - data_hipp.sampleinfo(1);
    data_hipp.time{1}    = 1/data_hipp.fsample+data_hipp.time{1}-data_hipp.time{1}(1);


    cfg.channel = intersect(cellstr(setdiff([data_hipp.label],settings.scalp_channels)),char(channels_hipp_ripples{isubject,:}));
    tmp         = cell2mat(alldat);
    cfg.channel = intersect(cfg.channel,[tmp.detectCh]);
    data_hipp   = ft_selectdata(cfg, data_hipp);

    %     data = ft_appenddata([],data,data_hipp); % doesn't work on mine.
    tmp = [];
    tmp.trial   = [data.trial{1,1};data_hipp.trial{1,1}];
    tmp.label   = [data.label(:);data_hipp.label(:)];

    data.trial  = [];
    data.label  = [];

    data.trial{1,1} = tmp.trial;
    data.label      = tmp.label;
    data            = rmfield(data, {'artifact', 'artifactdef'});
    data_hipp       = [];
    tmp             = [];


    %% create onset matrices (remove last xxx ms to ensure ripple event in trial)

    onsetmat = [onsets; onsets+(RT'.*data.fsample)-(settings.time_to_excl_RT*data.fsample)]';

    %% pick trials with ripples (optional to select long and short duration ripples)

    trl_ripple  = [];
    cnt         = 0;

    for ichannel = 1:numel(alldat)

        evs     = alldat{ichannel}.evtIndiv.maxTime;
        envSum  = alldat{ichannel}.evtIndiv.envSum;

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


    %% [optional] subtract erp from data (to not mix erp phase-reset with PAC).
    if settings.del_erp == 1

        cfg             = [];
        cfg.keeptrials  = 'yes';
        cfg.removemean  = 'no';
        data_ripples    = ft_timelockanalysis(cfg,data_ripples);

        for ichan = 1:numel(data_ripples.label)

            % compute ERP
            erp = squeeze(mean(data_ripples.trial(:,ichan,:),1));

            % compute induced power by subtracting ERP from each trial
            data_ripples.trial(:,ichan,:) = squeeze(data_ripples.trial(:,ichan,:)) - repmat(erp',[size(data_ripples.trial,1),1]);

        end

    end

    %% peak frequency

    % select only correct trials
    trlinfo = cell2mat(data_ripples.trialinfo);
    sel = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.Memory],'SourceCorrect');

    cfg                     = [];
    cfg.trials              = sel;
    data_ripples_correct    = ft_selectdata(cfg,data_ripples);

    % for theta
    tfr_cfg             = [];
    tfr_ripples = [];
    tfr_cfg.method      = 'mtmfft';
    tfr_cfg.output      = 'fooof_peaks';
    tfr_cfg.keeptrials  = 'no';

    tfr_cfg.taper       = 'hanning';
    tfr_cfg.foi         = settings.FOI_theta;
    tfr_cfg.pad         = 'nextpow2';
    tfr_cfg.width       = settings.cycles_theta;


    tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
    freq_theta_all         = ft_freqanalysis(tfr_cfg,data_ripples_correct);

    % for gamma

    tfr_cfg             = [];

    tfr_cfg.method      = 'mtmfft';
    tfr_cfg.output      = 'fooof_peaks';
    tfr_cfg.keeptrials  = 'no';

    tfr_cfg.taper       = 'hanning';
    tfr_cfg.foi         = settings.FOI_gamma;
    tfr_cfg.pad         = 'nextpow2';
    tfr_cfg.width       = settings.cycles_gamma;


    tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
    freq_gamma_all    = ft_freqanalysis(tfr_cfg,data_ripples_correct);

    %     cfg = [];
    % cfg.parameter = 'powspctrm';
    % cfg.appenddim = 'rpt';   % "repeat" dimension = trials
    %
    % freq_theta_all = ft_appendfreq(cfg, subtracted_osc_theta{:});
    % freq_gamma_all = ft_appendfreq(cfg, subtracted_osc_gamma{:});


    tfr_ripples.trialinfo   = data_ripples.trialinfo;

    tfr_ripples.freq    = [settings.FOI_theta,settings.FOI_gamma];
    tfr_ripples.label   = data_ripples.label;


    pow_ripples     = cat(2,freq_theta_all.powspctrm,freq_gamma_all.powspctrm);



    % pick cortical and hippocampal channels to find peak power for hipp and cortex
    these_hipp = intersect(tfr_ripples.label,intersect(cellstr(setdiff([data.label],settings.scalp_channels)),char(channels_hipp_ripples{isubject,:})));
    these_cort = setdiff(tfr_ripples.label,these_hipp);

    hipp_chann = ismember(tfr_ripples.label,these_hipp);
    cort_chann = ismember(tfr_ripples.label,these_cort);

    theta_FOI       = nearest(tfr_ripples.freq,3):nearest(tfr_ripples.freq,8);
    gamma_FOI       = nearest(tfr_ripples.freq,40):nearest(tfr_ripples.freq,140);

    theta_freq = settings.PAC_phase(theta_FOI);



    [~,max_theta] = max(pow_ripples(hipp_chann,theta_FOI)');
    [~,max_gamma] = max(pow_ripples(cort_chann,gamma_FOI)');



    PAC{isubject}.max_theta = theta_freq(max_theta);
    PAC{isubject}.max_gamma = settings.PAC_power(max_gamma);

    [~,max_theta] = max(nanmean(pow_ripples(hipp_chann,theta_FOI)));
    [~,max_gamma] = max(nanmean(pow_ripples(cort_chann,gamma_FOI)));

    PAC{isubject}.max_theta_average = theta_freq(max_theta);
    PAC{isubject}.max_gamma_average = settings.PAC_power(max_gamma);


    these_hipp = intersect(data_ripples.label,intersect(cellstr(setdiff([data.label],settings.scalp_channels)),char(channels_hipp_ripples{isubject,:})));
    these_cort = setdiff(data_ripples.label,these_hipp);


    fs = data_ripples_correct.fsample;

    % --- identify cortical channels (same logic you used for PAC) ---

    all_labels   = data_ripples_correct.label;

    % you already had these earlier when finding peaks:
    % these_hipp = intersect(...);
    % these_cort = setdiff(tfr_ripples.label,these_hipp);

    is_cort = ismember(all_labels, these_cort);
    cort_labels = all_labels(is_cort);

    % gamma peaks per cortical channel (MUST be aligned to cort_labels order)
    gamma_peaks_chan = PAC{isubject}.max_gamma(:);   % one per cortical channel

    gamma_shape = [];
    gamma_shape.chan_label      = cort_labels;
    PAC{isubject}.chan_label = cort_labels;
    gamma_shape.peak_trough_ratio = [];
    gamma_shape.rise_decay_ratio  = [];
    gamma_shape.skew             = [];
    gamma_shape.kurt             = [];
    gamma_shape.wave_t = [];           % will fill once we know window length
    gamma_shape.wave_mean = cell(1,numel(cort_labels));
    gamma_shape.wave_sem  = cell(1,numel(cort_labels));


    cfg = [];
    cfg.latency = [-1 1];
    data_ripples_correct = ft_selectdata(cfg,data_ripples_correct);

    k = 3;  % samples around extremum for sharpness; adjust if you want

    for ic = 1:numel(cort_labels)

        this_label = cort_labels{ic};
        f0 = gamma_peaks_chan(ic);           % centre gamma frequency for this channel

        if isnan(f0)
            continue; % skip if no peak
        end

        % --- design band-pass filter around peak gamma (±10 Hz) ---
        f_low  = max(30,  f0 - 10);          % keep within sensible gamma range
        f_high = min(150, f0 + 10);
        if f_low >= f_high
            continue;
        end

        % simple FIR or IIR; here: 4th order Butterworth
        [b, a] = butter(4, [f_low f_high] / (fs/2), 'bandpass');

        % --- concatenate trials in time for this channel ---
        sig_all = [];

        for itrial = 1:numel(data_ripples_correct.trial)
            trial_data = data_ripples_correct.trial{itrial};   % [nChan x nTime]
            chan_idx   = find(strcmp(all_labels, this_label));

            if isempty(chan_idx), continue; end

            sig = trial_data(chan_idx, :);     % 1 x nTime
            sig_filt = filtfilt(b, a, double(sig));   % zero-phase
            sig_all  = [sig_all sig_filt];     % concatenate across trials
        end

        if isempty(sig_all)
            continue;
        end

        % --- waveform-shape metrics on concatenated gamma signal ---

        sig = sig_all(:)';   % row vector

        % peaks & troughs
        [peak_amp, peak_idx]   = findpeaks(sig);
        [trough_amp, trough_idx] = findpeaks(-sig);
        trough_amp = -trough_amp;

        % peak & trough sharpness
        peak_sharp   = nan(size(peak_idx));
        trough_sharp = nan(size(trough_idx));

        for i = 1:numel(peak_idx)
            idx = peak_idx(i);
            if idx <= k || idx >= numel(sig)-k, continue; end
            neigh = sig([idx-k:idx-1 idx+1:idx+k]);
            peak_sharp(i) = peak_amp(i) - mean(neigh);
        end

        for i = 1:numel(trough_idx)
            idx = trough_idx(i);
            if idx <= k || idx >= numel(sig)-k, continue; end
            neigh = sig([idx-k:idx-1 idx+1:idx+k]);
            trough_sharp(i) = mean(neigh) - trough_amp(i);
        end

        peak_trough_ratio = nanmean(peak_sharp) / nanmean(trough_sharp);

        % rise–decay asymmetry
        rise_times  = [];
        decay_times = [];

        for i = 1:numel(trough_idx)
            next_peak = peak_idx(find(peak_idx > trough_idx(i), 1, 'first'));
            if ~isempty(next_peak)
                rise_times(end+1) = next_peak - trough_idx(i);
            end
        end

        for i = 1:numel(peak_idx)
            next_trough = trough_idx(find(trough_idx > peak_idx(i), 1, 'first'));
            if ~isempty(next_trough)
                decay_times(end+1) = next_trough - peak_idx(i);
            end
        end

        rise_decay_ratio = nanmean(rise_times) / nanmean(decay_times);

        % skew / kurtosis
        gamma_skew = skewness(sig);
        gamma_kurt = kurtosis(sig);

        % store
        gamma_shape.peak_trough_ratio(ic) = peak_trough_ratio;
        gamma_shape.rise_decay_ratio(ic)  = rise_decay_ratio;
        gamma_shape.skew(ic)              = gamma_skew;
        gamma_shape.kurt(ic)              = gamma_kurt;


        % now average gamma

        min_peak_dist = round((fs / f0) * 0.5);   % half a cycle in samples
        [peak_amp, peak_idx] = findpeaks(sig, 'MinPeakDistance', min_peak_dist);


        n_cycles = 2;  % window size in cycles (e.g., +/- 1 cycle around peak)
        half_win = round((fs / f0) * (n_cycles/2));   % samples either side
        win = -half_win:half_win;

        % time axis (seconds)
        if isempty(gamma_shape.wave_t)
            gamma_shape.wave_t = win / fs;
        end

        % collect snippets around peaks
        valid = peak_idx > half_win & peak_idx <= (numel(sig) - half_win);
        peak_idx_valid = peak_idx(valid);

        snips = nan(numel(peak_idx_valid), numel(win));
        for ii = 1:numel(peak_idx_valid)
            idx = peak_idx_valid(ii);
            snips(ii,:) = sig(idx + win);
        end

        % normalise each snippet if you want shape only (optional)
        % snips = snips ./ max(abs(snips),[],2);


        gamma_shape.wave_mean{ic} = mean(snips, 1, 'omitnan');
        gamma_shape.wave_sem{ic}  = std(snips, [], 1, 'omitnan') ./ sqrt(size(snips,1));

        n_phase = 100;  % number of points per cycle (choose once)

        phase_old = linspace(-pi, pi, numel(gamma_shape.wave_mean{ic}));
        phase_new = linspace(-pi, pi, n_phase);

        gamma_shape.wave_mean_phase{ic} = interp1( ...
            phase_old, gamma_shape.wave_mean{ic}, phase_new, 'linear');

        gamma_shape.wave_sem_phase{ic} = interp1( ...
            phase_old, gamma_shape.wave_sem{ic}, phase_new, 'linear');




    end

    PAC{isubject}.peak_trough_ratio = gamma_shape.peak_trough_ratio;
    PAC{isubject}.rise_decay_ratio = gamma_shape.rise_decay_ratio;
    PAC{isubject}.skew = gamma_shape.skew;
    PAC{isubject}.kurt = gamma_shape.kurt;

    wave_mat = [];
    for ic = 1:numel(cort_labels)
        if ~isempty(gamma_shape.wave_mean_phase{ic})
            wave_mat(end+1,:) = gamma_shape.wave_mean_phase{ic};
        end
    end

    gamma_shape.wave_grand_mean = mean(wave_mat, 1, 'omitnan');
    gamma_shape.wave_grand_sem  = std(wave_mat, [], 1, 'omitnan') ./ sqrt(size(wave_mat,1));


    PAC{isubject}.gamma_wave_t          = gamma_shape.wave_t;
    PAC{isubject}.gamma_wave_mean       = gamma_shape.wave_mean;
    PAC{isubject}.gamma_wave_sem        = gamma_shape.wave_sem;
    PAC{isubject}.gamma_wave_grand_mean = gamma_shape.wave_grand_mean;
    PAC{isubject}.gamma_wave_grand_sem  = gamma_shape.wave_grand_sem;


end

delete(gcp);

return


%% stats

phase_to_plot       = [];
pow_to_plot         = [];
phase_to_plot_avg       = [];
pow_to_plot_avg         = [];

gamma_to_plot = [];
% figure;
for participant = 1:numel(subjects)
    %     hold on
    %     subplot(4,3,participant)

    phase_to_plot{participant}       = PAC{participant}.max_theta;
    pow_to_plot{participant}         = PAC{participant}.max_gamma;

    phase_to_plot_avg(participant,:)       = PAC{participant}.max_theta_average;
    pow_to_plot_avg(participant,:)         = PAC{participant}.max_gamma_average;

    %     hist(PAC{participant}.max_theta)

    gamma_to_plot(participant,:)         = PAC{participant}.gamma_wave_grand_mean;

end

m = mean(gamma_to_plot);
s = std(gamma_to_plot)/sqrt(12);

x = linspace(-pi, pi, numel(m));
x_label = 'Gamma phase (rad)';
set(gca,'XTick',[-pi 0 pi],'XTickLabel',{'-\pi','0','\pi'});
figure; hold on

fill([x fliplr(x)], ...
     [m+s fliplr(m-s)], ...
     [0.7 0.7 0.7], ...
     'EdgeColor','none', ...
     'FaceAlpha',0.5);

plot(x, m, 'k', 'LineWidth', 2);

xlabel(x_label);
ylabel('Gamma amplitude (a.u.)');
title('Average cortical gamma waveform');

box off

plot(x, gamma_to_plot', 'Color', [0.6 0.6 0.6 0.3]);
hold on
plot(x, m, 'k', 'LineWidth', 2);




theta_dispersion   = [];
gamma_dispersion   = [];

for participant = 1:numel(subjects)

    theta_vals = PAC{participant}.max_theta;
    theta_avg  = PAC{participant}.max_theta_average;

    gamma_vals = PAC{participant}.max_gamma;
    gamma_avg  = PAC{participant}.max_gamma_average;

    % --- Theta dispersion ---
    theta_dispersion.std(participant)  = std(theta_vals);
    theta_dispersion.mad(participant)  = mad(theta_vals,1);   % robust measure
    theta_dispersion.range(participant,:) = [min(theta_vals), max(theta_vals)];
    theta_dispersion.within2Hz(participant) = mean(abs(theta_vals - theta_avg) <= 2);

    % --- Gamma dispersion ---
    gamma_dispersion.std(participant)  = std(gamma_vals);
    gamma_dispersion.mad(participant)  = mad(gamma_vals,1);
    gamma_dispersion.range(participant,:) = [min(gamma_vals), max(gamma_vals)];
    gamma_dispersion.within10Hz(participant) = mean(abs(gamma_vals - gamma_avg) <= 10);

end


mean(theta_dispersion.std)
mean(theta_dispersion.mad)
mean(theta_dispersion.within2Hz)



mean(gamma_dispersion.std)
mean(gamma_dispersion.mad)
mean(gamma_dispersion.within10Hz)




%% waveform metrics

all_peaktrough = [];
all_risedecay = [];
all_skew = [];
all_kurt = [];
all_gammafreq = [];

for s = 1:numel(PAC)
    all_peaktrough = [all_peaktrough; PAC{s}.peak_trough_ratio(:)];
    all_risedecay = [all_risedecay; PAC{s}.rise_decay_ratio(:)];
    all_skew = [all_skew; PAC{s}.skew(:)];
    all_kurt = [all_kurt; PAC{s}.kurt(:)];
    all_gammafreq = [all_gammafreq; PAC{s}.max_gamma(:)];
end



mean_peaktr = mean(all_peaktrough,'omitnan')
std_peaktr  = std(all_peaktrough,'omitnan')

mean_risedecay = mean(all_risedecay,'omitnan')
std_risedecay  = std(all_risedecay,'omitnan')

mean_skew = mean(all_skew,'omitnan')
std_skew  = std(all_skew,'omitnan')

mean_kurt = mean(all_kurt,'omitnan')
std_kurt  = std(all_kurt,'omitnan')



[r_pt, p_pt] = corr(all_peaktrough, all_gammafreq, 'rows','complete');
[r_rd, p_rd] = corr(all_risedecay, all_gammafreq, 'rows','complete');
[r_sk, p_sk] = corr(all_skew, all_gammafreq, 'rows','complete');
[r_ku, p_ku] = corr(all_kurt, all_gammafreq, 'rows','complete');

fprintf('Peak–trough corr = %.3f, p=%.3g\n', r_pt, p_pt);
fprintf('Rise–decay corr = %.3f, p=%.3g\n', r_rd, p_rd);
fprintf('Skew corr = %.3f, p=%.3g\n', r_sk, p_sk);
fprintf('Kurtosis corr = %.3f, p=%.3g\n', r_ku, p_ku);



figure;
subplot(2,2,1); histogram(all_peaktrough); title('Peak–Trough Ratio');
subplot(2,2,2); histogram(all_risedecay); title('Rise–Decay Ratio');
subplot(2,2,3); histogram(all_skew); title('Skewness');
subplot(2,2,4); histogram(all_kurt); title('Kurtosis');
