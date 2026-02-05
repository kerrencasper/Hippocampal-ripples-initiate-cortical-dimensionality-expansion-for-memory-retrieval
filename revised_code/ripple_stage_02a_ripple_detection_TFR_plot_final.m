%% detect ripples and plot TFR plot of detected ripples
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults



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

settings.healthyhemi        = {'R' 'LR' 'R' 'L' 'L' 'R' 'R' 'R' 'R' 'R' 'R' 'LR'};

addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/MVPA-Light-master'))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/help_functions'))
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))

%% pre-decoding settings
% ripple extraction

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1; % remove co-occuring ripples
settings.time_to_excl_RT            = .25; % exclude last 250 ms of trials, to make sure ripple event was in trial
settings.solo_ripple                = 1; % pick one (maxEnv) ripple per trial if multiple ripple events are found
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur, eg [.5 1.5]

%% decoding settings

% data preprocessing
settings.ori_fsample            = 1000; % original sample frequency
settings.do_resample            = 500; % [] or sample frequency

settings.TOI_train              = [-.5 3];
settings.timesteps_train        = ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample); % in s. if you want it to take less sample points multiply [e.g., ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample)*2

settings.TOI_test               = [-6 6]; % time around ripple
settings.timesteps_test         = ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample); % in s

settings.classifier             = 'lda';
settings.metric                 = 'auc';

%% channel settings

settings.channel_sel            = 2; % 1 exclude hippo, 2 only hippo, 3 all channels

%% start for loop


%% TFR settings

% define peak spectral power


settings.FOI       = 60:2:120;
settings.cycles                         = ceil(settings.FOI * 0.5); % ~ 500 ms for each frequency
settings.cycles(settings.cycles < 5)    = 5;
% for TF after peak estimate

settings.bs_correct             = 0; % 0 = none | 1 = % baseline change | 2 = zscore each trial
settings.bs_period              = [-.5 -.1]; % for option 1 above
settings.bs_period_ripple       = [-1.25 -1; 1.25 1]; % [-.5 -.1]


settings.TOI                    = [-1 1]; % time around ripple
settings.timesteps              = .020; % for TFR

settings.TFR_padding = ceil(max((1./settings.FOI).*settings.cycles)/2); % to ensure sufficient time for spectral resolution is available


%% start for loop

timeaxis = settings.TOI_test(1):settings.timesteps_test:settings.TOI_test(2);
freqaxis = settings.TOI_train(1):settings.timesteps_train:settings.TOI_train(2);


perf    = cell(1,numel(subjects));
channs  = cell(1,numel(subjects));

numWorkers = 6; %

parpool('local', numWorkers);

tic

TF = cell(1,numel(subjects));

parfor isubject = 1:numel(subjects)

    fprintf('processing subject %01d/%02d\n',isubject,numel(subjects));

    %% LOAD data

    data_tmp_1 = [];
    data_tmp_2 = [];

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

    trlinfo = cell2mat(data_ripples.trialinfo);

    sel1 = ismember([trlinfo.ExpPhase],'retrieval') 
   
    cfg             = [];
    cfg.trials      = sel1;
    data_ripples    = ft_selectdata(cfg,data_ripples)

    trlinfo = cell2mat(data_stimuli.trialinfo);

    sel1 = ismember([trlinfo.ExpPhase],'retrieval') 
   
    cfg             = [];
    cfg.trials      = sel1;
    data_stimuli    = ft_selectdata(cfg,data_stimuli)
    %% power estimate

    tfr_cfg             = [];
    tfr_cfg.output      = 'pow';
    tfr_cfg.method      = 'wavelet';
    tfr_cfg.taper       = 'hanning';
    tfr_cfg.foi         = settings.FOI;
    tfr_cfg.pad         = 'nextpow2';
    tfr_cfg.width       = settings.cycles;
    tfr_cfg.keeptrials  = 'yes';
    tfr_cfg.polyremoval = 1;


    tfr_cfg.toi         = settings.TOI(1):settings.timesteps:settings.TOI(2);
    tfr_ripples_gamma   = ft_freqanalysis(tfr_cfg,data_ripples);

    cfg                 = [];
    cfg.baseline        = [-1 1];
    cfg.baselinetype    = 'relchange';
    tfr_ripples_gamma   = ft_freqbaseline(cfg,tfr_ripples_gamma);


    pow_ripples     = [];
    pow_ripples     = tfr_ripples_gamma.powspctrm;





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

    hipp_chann = ismember(data_ripples.label,these_hipp);

    pow_ripples = squeeze(mean(pow_ripples(:,hipp_chann,:,:)));

  
    TF{isubject}.all   = pow_ripples;
    TF{isubject}.time   = tfr_ripples_gamma.time;
    TF{isubject}.freq   = tfr_ripples_gamma.freq;


end

delete(gcp);

return


%% plot data
to_plot_ripple = [];
peak_power = [];
for isubject = 1:numel(TF)
    to_plot_ripple(isubject,:,:) = squeeze(mean(TF{isubject}.all));
    [idx,peak_pow] = max(squeeze(mean(mean(TF{isubject}.all(:,:,51),3),1)));
    peak_power(isubject,:) = peak_pow;
end

peak_power = TF{1}.freq(peak_power);
mean_peak = mean(peak_power);
std_peak = std(peak_power);
to_plot_ripple = squeeze(mean(to_plot_ripple));

timeaxis   = TF{1}.time;
freqaxis    = TF{1}.freq;

figure;
imagesc(...
    timeaxis,...
    freqaxis,...
    to_plot_ripple);
colormap(hot)
caxis([0 .7])


axis xy
set(gca,'FontSize',20)
xlabel('ripple time')
ylabel('frequency power (Hz)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');

%%
figure;
for isubject = 1:numel(TF)
    subplot(4,3,isubject)
    to_plot_ripple = [];
    to_plot_ripple = squeeze(mean(TF{isubject}.all));

    imagesc(...
        timeaxis,...
        freqaxis,...
        to_plot_ripple);
    colormap(hot)
    caxis([0 .7])


    axis xy
    set(gca,'FontSize',20)
    xlabel('ripple time')
    ylabel('frequency power (Hz)')
    set(gca,'TickDir','out')

    hcb = colorbar('Location','EastOutside');
end