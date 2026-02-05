%%
% [ripple_stage_05a_PAC] find ripples in hippocampal channels,
% extract and realign data based on ripple events.
% estimate phase-amplitude coupling between HC and cortex
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear
close all


%% path settings

settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

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

numWorkers = 4; %

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
if settings.peak_frequency == 1
    

    % select only correct trials
    trlinfo = cell2mat(data_ripples.trialinfo);
    sel = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.Memory],'SourceCorrect');

    cfg                     = [];
    cfg.trials              = sel;
    data_ripples_correct    = ft_selectdata(cfg,data_ripples);

    % for theta
    tfr_cfg             = [];
    if settings.peak_method == 1
        tfr_cfg.output      = 'fractal';
        tfr_cfg.method      = 'irasa';
        tfr_cfg.keeptrials  = 'yes';


        tfr_cfg.taper       = 'hanning';
        tfr_cfg.foi         = settings.FOI_theta;
        tfr_cfg.pad         = 'nextpow2';
        tfr_cfg.width       = settings.cycles_theta;


        tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
        tfr_ripples_fraq    = ft_freqanalysis(tfr_cfg,data_ripples_correct);

        tfr_cfg.output      = 'original';
        tfr_ripples         = ft_freqanalysis(tfr_cfg, data_ripples_correct);

        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2-x1';
        subtracted_osc_theta    = ft_math(cfg, tfr_ripples_fraq, tfr_ripples);

    elseif settings.peak_method == 2
        tfr_cfg.method      = 'mtmfft';
        tfr_cfg.output      = 'fooof_peaks';
        tfr_cfg.keeptrials  = 'no';

        tfr_cfg.taper       = 'hanning';
        tfr_cfg.foi         = settings.FOI_theta;
        tfr_cfg.pad         = 'nextpow2';
        tfr_cfg.width       = settings.cycles_theta;


        tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
        subtracted_osc_theta         = ft_freqanalysis(tfr_cfg,data_ripples_correct);

    end

    % for gamma

    tfr_cfg             = [];
    if settings.peak_method == 1
        tfr_cfg.output      = 'fractal';
        tfr_cfg.method      = 'irasa';
        tfr_cfg.keeptrials  = 'yes';

        tfr_cfg.taper       = 'hanning';
        tfr_cfg.foi         = settings.FOI_gamma;
        tfr_cfg.pad         = 'nextpow2';
        tfr_cfg.width       = settings.cycles_gamma;


        tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
        tfr_ripples_fraq    = ft_freqanalysis(tfr_cfg,data_ripples_correct);

        tfr_cfg.output      = 'original';
        tfr_ripples         = ft_freqanalysis(tfr_cfg, data_ripples_correct);

        cfg               = [];
        cfg.parameter     = 'powspctrm';
        cfg.operation     = 'x2-x1';
        subtracted_osc_gamma    = ft_math(cfg, tfr_ripples_fraq, tfr_ripples);

    elseif settings.peak_method == 2
        tfr_cfg.method      = 'mtmfft';
        tfr_cfg.output      = 'fooof_peaks';
        tfr_cfg.keeptrials  = 'no';

        tfr_cfg.taper       = 'hanning';
        tfr_cfg.foi         = settings.FOI_gamma;
        tfr_cfg.pad         = 'nextpow2';
        tfr_cfg.width       = settings.cycles_gamma;


        tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
        subtracted_osc_gamma    = ft_freqanalysis(tfr_cfg,data_ripples_correct);
    end


    tfr_ripples.trialinfo   = data_ripples.trialinfo;

    tfr_ripples.freq    = [settings.FOI_theta,settings.FOI_gamma];
    tfr_ripples.label   = data_ripples.label;

    if settings.peak_method == 1
        pow_ripples     = cat(3,subtracted_osc_theta.powspctrm,subtracted_osc_gamma.powspctrm);
    elseif settings.peak_method == 2
        pow_ripples     = cat(2,subtracted_osc_theta.powspctrm,subtracted_osc_gamma.powspctrm);
    end


    % pick cortical and hippocampal channels to find peak power for hipp and cortex
    these_hipp = intersect(tfr_ripples.label,intersect(cellstr(setdiff([data.label],settings.scalp_channels)),char(channels_hipp_ripples{isubject,:})));
    these_cort = setdiff(tfr_ripples.label,these_hipp);

    hipp_chann = ismember(tfr_ripples.label,these_hipp);
    cort_chann = ismember(tfr_ripples.label,these_cort);

    theta_FOI       = nearest(tfr_ripples.freq,3):nearest(tfr_ripples.freq,8);
    gamma_FOI       = nearest(tfr_ripples.freq,40):nearest(tfr_ripples.freq,140);

    theta_freq = settings.PAC_phase(theta_FOI);


    if settings.peak_method == 1
        [~,max_theta] = max(squeeze(nanmean(nanmean(pow_ripples(:,hipp_chann,theta_FOI)))));
        [~,max_gamma] = max(squeeze(nanmean(nanmean(pow_ripples(:,cort_chann,gamma_FOI)))));
    elseif settings.peak_method == 2
        [~,max_theta] = max(squeeze(nanmean(pow_ripples(hipp_chann,theta_FOI))));
        [~,max_gamma] = max(squeeze(nanmean(pow_ripples(cort_chann,gamma_FOI))));
    end


    PAC{isubject}.max_theta = theta_freq(max_theta);
    PAC{isubject}.max_gamma = settings.PAC_power(max_gamma);

    freq_theta = [theta_freq(max_theta)-2:theta_freq(max_theta)+2];
    freq_gamma = [settings.PAC_power(max_gamma)-10:settings.PAC_power(max_gamma)+10];

    cycles                = ceil([freq_theta,freq_gamma] * 0.5);
    cycles(cycles < 5)    = 5;


    tfr_cfg             = [];
    tfr_cfg.output      = 'fourier';
    tfr_cfg.method      = 'wavelet';
    tfr_cfg.taper       = 'hanning';
    tfr_cfg.foi         = freq_theta;
    tfr_cfg.pad         = 'nextpow2';
    tfr_cfg.width       = cycles(1:numel(freq_theta));
    tfr_cfg.keeptrials  = 'yes';
    tfr_cfg.polyremoval = 1;
    
    tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
    tfr_ripples_theta   = ft_freqanalysis(tfr_cfg,data_ripples);

    tfr_cfg             = [];
    tfr_cfg.tapsmofrq   = freq_gamma/4;
    tfr_cfg.output      = 'pow';
    tfr_cfg.method      = 'mtmconvol';
    tfr_cfg.taper       = 'dpss';
    tfr_cfg.foi         = freq_gamma;
    tfr_cfg.t_ftimwin   =  0.2 * ones(length(tfr_cfg.foi), 1)
    tfr_cfg.pad         = 'nextpow2';
    tfr_cfg.keeptrials  = 'yes';
    tfr_cfg.polyremoval = 1;
    
    tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
    tfr_ripples_gamma   = ft_freqanalysis(tfr_cfg,data_ripples);


elseif settings.peak_frequency == 2
    

    tfr_cfg             = [];
    tfr_cfg.output      = 'fourier';
    tfr_cfg.method      = 'wavelet';
    tfr_cfg.taper       = 'hanning';
    tfr_cfg.foi         = [settings.PAC_phase];
    tfr_cfg.pad         = 'nextpow2';
    tfr_cfg.width       = settings.cycles(1:numel(settings.PAC_phase));
    tfr_cfg.keeptrials  = 'yes';
    tfr_cfg.polyremoval = 1;
    
    tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
    tfr_ripples_theta   = ft_freqanalysis(tfr_cfg,data_ripples);


    tfr_cfg             = [];
    tfr_cfg.tapsmofrq   = settings.PAC_power/4;
    tfr_cfg.output      = 'pow';
    tfr_cfg.method      = 'mtmconvol';
    tfr_cfg.taper       = 'dpss';
    tfr_cfg.foi         = [settings.PAC_power];
    tfr_cfg.t_ftimwin   =  0.2 * ones(length(tfr_cfg.foi), 1)
    tfr_cfg.pad         = 'nextpow2';
    tfr_cfg.keeptrials  = 'yes';
    tfr_cfg.polyremoval = 1;
    
    tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
    tfr_ripples_gamma   = ft_freqanalysis(tfr_cfg,data_ripples);


end

    pow_ripples     = [];
    plv_ripples     = [];
    pow_ripples     = tfr_ripples_gamma.powspctrm;
    plv_ripples     = angle(tfr_ripples_theta.fourierspctrm);

    
    %% [optional] baseline correction
    
    if settings.bs_correct == 1 % based on prestim baselin
        
        pretrig      = round((abs(settings.bs_period(1))+settings.TFR_padding) * data.fsample);
        posttrig     = round((abs(settings.bs_period(2))+settings.TFR_padding) * data.fsample);
        
        cfg          = [];
        cfg.trl      = [onsets'-pretrig onsets'+posttrig -pretrig*ones(numel(onsets),1)];
        
        data_BL      = ft_redefinetrial(cfg,data);
        
        data_BL.trialinfo = [];
        
        for itrial = 1:numel(data_BL.trial)
            
            data_BL.trialinfo{itrial,1} = trialinfo(itrial);
        end
        
        %--- match BL events and channels with ripple-locked events
        
        trlinfo_BL              = cell2mat(data_BL.trialinfo);
        trlinfo_ripplelocked    = cell2mat(data_ripples.trialinfo);
        
        sel_trl                 = intersect([trlinfo_BL.EventNumber],[trlinfo_ripplelocked.EventNumber]);
        sel_chan                = intersect(data_BL.label,data_ripples.label);
        
        cfg                     = [];
        cfg.trials              = sel_trl;
        cfg.channel             = sel_chan;
        data_BL                 = ft_selectdata(cfg, data_BL);
        
        %--- calculate TFR for BL
        disp('frequency analysis of baseline data...')
        tfr_cfg.toi      = settings.bs_period(1):settings.timesteps:settings.bs_period(2);
        tfr_BL           = ft_freqanalysis(tfr_cfg,data_BL);
        
        tfr_BL.trialinfo = data_ripples.trialinfo;
        
        %--- power for baseline period
        
%         pow_BL   = abs(tfr_BL.fourierspctrm) .^2;
        pow_BL   = tfr_BL.powspctrm;
        
        
        dims            = size(pow_ripples);
        
        blmat           = repmat(mean(pow_BL,numel(dims)),[ones(1,numel(dims)-1) dims(end)]);
        pow_ripples     = (pow_ripples - blmat)./blmat;
        
        blmat  = [];
        pow_BL = [];
        
    elseif settings.bs_correct == 2 % zscore within each trial (that is, dim 4), independent of channel or frequency
        
        pow_ripples     = zscore(pow_ripples,[],4);
        
    end


    tfr_BL      = [];
    
    %%

    these_hipp = intersect(data_ripples.label,intersect(cellstr(setdiff([data.label],settings.scalp_channels)),char(channels_hipp_ripples{isubject,:})));
    these_cort = setdiff(data_ripples.label,these_hipp);
    
    hipp_chann = ismember(data_ripples.label,these_hipp);
    cort_chann = ismember(data_ripples.label,these_cort);
    
    pow_ripples = pow_ripples(:,cort_chann,:,:);
    plv_ripples = plv_ripples(:,hipp_chann,:,:);

    %% [optional] pick ripples from a single channel
    if settings.ripples_most_chan
        
        n_ripples = [];
        
        for ichannel = 1:numel(alldat)
            n_ripples(ichannel) = alldat{ichannel}.evtIndiv.numEvt;
        end
        
        [~,maxchannel]  = max(n_ripples);
        
        plv_ripples     = plv_ripples(:,maxchannel,:,:);
    end
    
    
    %% calculate PAC

    trlinfo = cell2mat(data_ripples.trialinfo);

    for icond = 1:1

        if icond == 1

            sel = ismember([trlinfo.ExpPhase],'retrieval') & ismember([trlinfo.Memory],'SourceCorrect');
          
            perf{isubject}.correct.trl_num_test = sum(sel);

        elseif icond == 2

            sel = ismember([trlinfo.ExpPhase],'retrieval') & (ismember([trlinfo.Memory],{'SourceIncorrect' ,'SourceDunno' 'dunno'}));

            perf{isubject}.incorrect.trl_num_test         = sum(sel);
        end

        pow_ripples_PAC = pow_ripples(sel,:,:,:);
        plv_ripples_PAC = plv_ripples(sel,:,:,:);

        TOI_PAC   = nearest(tfr_ripples_theta.time,settings.PAC_TOI(1)):nearest(tfr_ripples_theta.time,settings.PAC_TOI(2));

        if settings.peak_frequency == 1

            FOI_phase = nearest(tfr_ripples_theta.freq,freq_theta(1)):nearest(tfr_ripples_theta.freq,freq_theta(end));
            FOI_power = nearest(tfr_ripples_gamma.freq,freq_gamma(1)):nearest(tfr_ripples_gamma.freq,freq_gamma(end));

        else

            FOI_phase = nearest(tfr_ripples_theta.freq,1):nearest(tfr_ripples_theta.freq,8);
            FOI_power = nearest(tfr_ripples_gamma.freq,30):nearest(tfr_ripples_gamma.freq,150);

        end

        PAC_tmp_hipp        = [];
        PAC_tmp_hipp_perm   = [];
        PAC_cort            = [];
        PAC_cort_perm       = [];
        PAC_pow             = [];
        PAC_pow_perm        = [];
        PAC_all             = [];
        PAC_perm_all        = [];

        if settings.PAC_or_MI == 1
        PAC_perm_all = [];
         PAC_all = zeros(numel(FOI_power),numel(FOI_phase));
            for ifreq = 1:numel(FOI_phase)
                ifreq
                PAC_pow = zeros(1,numel(FOI_power));
                for ipow = 1:numel(FOI_power)
                    
                    PAC_cort = zeros(size(pow_ripples_PAC,2),1);
                    for ichan = 1:size(pow_ripples_PAC,2)
                        PAC_tmp_hipp = zeros(size(plv_ripples_PAC,2),1);
                        for iichan = 1:size(plv_ripples_PAC,2)

                            % extract temporally localized power and phase
                            ampl = squeeze(pow_ripples_PAC(:,ichan,FOI_power(ipow),TOI_PAC));
                            phas = squeeze(plv_ripples_PAC(:,iichan,FOI_phase(ifreq),TOI_PAC));

                            obsPAC = abs(mean(ampl(:).*exp(1i*phas(:))));
                            PAC_perm = [];
                            if any(settings.PAC_nusurro)

                                for i=1:settings.PAC_nusurro

                                    rand_trial  = randperm(size(ampl,1));
                                    ampl        = ampl(rand_trial,:);

                                    PAC_perm(i,:) = abs(mean(ampl(:).*exp(1i*phas(:))));
                                end
                            end

                            PAC_tmp_hipp(iichan,:)       = obsPAC;
                            if any(settings.PAC_nusurro)
                                PAC_tmp_hipp_perm(iichan,:)  = mean(PAC_perm);
                            end

                        end
                        PAC_cort(ichan,:)        = mean(PAC_tmp_hipp);
                        if any(settings.PAC_nusurro)
                            PAC_cort_perm(ichan,:)   = mean(PAC_tmp_hipp_perm);
                        end
                    end
                    PAC_pow(:,ipow)      = mean(PAC_cort);
                    if any(settings.PAC_nusurro)
                        PAC_pow_perm(:,ipow) = mean(PAC_cort_perm);
                    end
                end
                PAC_all(:,ifreq)         = PAC_pow;
                if any(settings.PAC_nusurro)
                    PAC_perm_all(:,ifreq)    = PAC_pow_perm;
                end

            end

        elseif settings.PAC_or_MI == 2
            PAC_perm_all = [];
            PAC_all = zeros(numel(FOI_power),numel(FOI_phase));
            for ifreq = 1:numel(FOI_phase)
                ifreq
                PAC_pow = zeros(1,numel(FOI_power));
                for ipow = 1:numel(FOI_power)
                    PAC_cort = zeros(size(pow_ripples_PAC,2),1);
                    for ichan = 1:size(pow_ripples_PAC,2)
                        PAC_tmp_hipp = zeros(size(plv_ripples_PAC,2),1);
                        for iichan = 1:size(plv_ripples_PAC,2)
                            % extract temporally localized power and phase

                            ampl = squeeze(pow_ripples_PAC(:,ichan,FOI_power(ipow),TOI_PAC));
                            phas = squeeze(exp(1i*plv_ripples_PAC(:,iichan,FOI_phase(ifreq),TOI_PAC)));

                            [~,binIdx]  = histc(angle(phas),settings.PAC_binEdges);

                            amplBin     = zeros(1,settings.PAC_bins);
                            for bin=1:settings.PAC_bins
                                if sum(sum(binIdx==bin))>0
                                    amplBin(bin) = mean(ampl(binIdx==bin));
                                end
                            end

                            % avoid negative numbers for later KL calculation. Add abs
                            % of minimum to all bins.
                            if any(amplBin<0)
                                amplBin = amplBin+abs(min(amplBin));
                            end
                            amplP = amplBin/sum(amplBin);

                            % compute Kullback-Leibler distance and modulation index (MI)
                            amplQ = ones(1,settings.PAC_bins)./settings.PAC_bins;

                            % in the special case where observed probability in a bin is 0, this tweak
                            % allows computing a meaningful KL distance nonetheless
                            if any(amplP==0)
                                amplP(amplP==0)=eps;
                            end

                            distKL  = sum(amplP.*log(amplP./amplQ));
                            PAC_tmp = distKL./log(settings.PAC_bins);

                            % surrogates
                            PAC_perm = [];
                            if any(settings.PAC_nusurro)

                                for i=1:settings.PAC_nusurro

                                    rand_trial  = randperm(size(ampl,1));
                                    ampl        = ampl(rand_trial,:);

                                    amplBin = zeros(1,settings.PAC_bins);
                                    for bin=1:settings.PAC_bins
                                        if sum(sum(binIdx==bin))>0
                                            amplBin(bin)=mean(ampl(binIdx==bin));
                                        end
                                    end

                                    if any(amplBin<0)
                                        amplBin = amplBin+abs(min(amplBin));
                                    end

                                    amplP = amplBin/sum(amplBin);

                                    % compute Kullback-Leibler distance and modulation index (MI)
                                    amplQ = ones(1,settings.PAC_bins)./settings.PAC_bins;

                                    if any(amplP==0)
                                        amplP(amplP==0)=eps;
                                    end

                                    distKL          = sum(amplP.*log(amplP./amplQ));
                                    PAC_perm(i,:)   = distKL./log(settings.PAC_bins);
                                end
                            end
                            PAC_tmp_hipp(iichan,:)       = PAC_tmp;
                            if any(settings.PAC_nusurro)
                                PAC_tmp_hipp_perm(iichan,:)  = mean(PAC_perm);
                            end

                        end
                        PAC_cort(ichan,:)        = mean(PAC_tmp_hipp);
                        if any(settings.PAC_nusurro)
                            PAC_cort_perm(ichan,:)   = mean(PAC_tmp_hipp_perm);
                        end
                    end
                    PAC_pow(:,ipow)      = mean(PAC_cort);
                    if any(settings.PAC_nusurro)
                        PAC_pow_perm(:,ipow) = mean(PAC_cort_perm);
                    end
                end
                PAC_all(:,ifreq)         = PAC_pow;
                if any(settings.PAC_nusurro)
                    PAC_perm_all(:,ifreq)    = PAC_pow_perm;
                end

            end
        end

        if icond == 1
            PAC{isubject}.correct   = PAC_all;
        elseif icond == 2
            PAC{isubject}.incorrect   = PAC_all;
        end
        if any(settings.PAC_nusurro)
            PAC{isubject}.permuted  = PAC_perm_all;
        end
    end

end

delete(gcp);

return


%% stats

if settings.peak_frequency == 1
timeaxis = -2:2;
freqaxis = -10:10;
end


phase_to_plot       = [];
pow_to_plot         = [];
for participant = 1:numel(subjects)

    phase_to_plot       = [phase_to_plot, PAC{participant}.max_theta];
    pow_to_plot         = [pow_to_plot, PAC{participant}.max_gamma];

end

boxplot([phase_to_plot',pow_to_plot'])

timeaxis = settings.PAC_phase;
freqaxis = settings.PAC_power;

correct_PAC = {};
baseline_PAC = {};
baseline_all = [];
correct_all = [];

phase_to_plot       = [];
pow_to_plot         = [];
for participant = 1:numel(subjects)

    phase_to_plot       = [phase_to_plot, PAC{participant}.max_theta];
    pow_to_plot         = [pow_to_plot, PAC{participant}.max_gamma];

end

boxplot([phase_to_plot',pow_to_plot'])

for isubject = 1:numel(subjects)
    
    % correct
    correct_PAC{1,isubject}                                 = struct;
    correct_PAC{1,isubject}.label                           = {'chan'};
    correct_PAC{1,isubject}.dimord                          = 'chan_freq_time';
    correct_PAC{1,isubject}.freq                            = freqaxis;
    correct_PAC{1,isubject}.time                            = timeaxis;
    

    correct_PAC{1,isubject}.powspctrm(1,:,:)                = PAC{1,isubject}.correct;
    correct_all(isubject,:,:) = PAC{1,isubject}.correct;
    
    % baseline
    baseline_PAC{1,isubject}                                = correct_PAC{1,isubject};
    baseline_PAC{1,isubject}.powspctrm(1,:,:)               = PAC{1,isubject}.permuted;
    baseline_all(isubject,:,:) = PAC{1,isubject}.permuted; 
    
end

% Run stats

cfg                     = [];
cfg.latency             = [timeaxis(1) timeaxis(end)];%[2 6]%[4 8];
cfg.frequency           = [freqaxis(1) freqaxis(end)]; %[60 80]%[10 20]
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

cfg.numrandomization    = 250;%1000;%'all';

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

% run stats
[Fieldtripstats] = ft_freqstatistics(cfg, correct_PAC{:}, baseline_PAC{:});

length(find(Fieldtripstats.mask==1))

%% plot (significant) t vals

stats_time = nearest(timeaxis,cfg.latency(1)):nearest(timeaxis,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

% tvalues
tvals = squeeze(Fieldtripstats.stat);

figure;
imagesc(...
    timeaxis,...
    freqaxis,...
    tvals.*squeeze(Fieldtripstats.mask));
colormap(hot)
caxis([0 max(abs(tvals(:)))])


axis xy
set(gca,'FontSize',12)
xlabel('frequency phase (Hz)')
ylabel('frequency power (Hz)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title(sprintf('MI %s correction, %s',cfg.correctm,cfg.method));
ylim([min(freqaxis(stats_freq)) max(freqaxis(stats_freq))])
xlim([min(timeaxis(stats_time)) max(timeaxis(stats_time))])

