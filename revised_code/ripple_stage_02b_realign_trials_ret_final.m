%% Find ripples and make not of the time points then realign the artifact rejected data based on the time points.
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear

%% path settings
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW','09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');


settings.data_dir           = [settings.base_path_castle,'output_data/ripple_detection/'];
settings.data_dir_arti      = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];
settings.save_dir           = [settings.base_path_castle,'output_data/realigned_data/'];

mean_all = [];
for isubject = 1:size(settings.subjects,1)

    fprintf('realigning subject %01d/%02d\n',isubject,size(settings.subjects,1));

    settings.SubjectID   = num2str(settings.SubjectIDs(isubject,:));

    %% Load data


    data_in                 = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_ripple_',num2str(settings.subjects(isubject,:))],'alldat','onsets_session');
    data_in_arti            = load([settings.data_dir_arti,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',num2str(settings.subjects(isubject,:))],'data','onsets_session');
    alldat                  = data_in.alldat;
    onsets_session          = data_in.onsets_session;
    data                    = data_in_arti.data;
    time_to_excl_RT         = .25; % in sec


    %% Find onsets

    onsets              = onsets_session - data.time{1}(1)*data.fsample;
    [numbers,strings]   = xlsread([settings.base_path_castle,'well01_behavior_all.xls']);
    strings             = strings(2:end,:);
    if isnan(numbers(1,1))
        numbers         = numbers(2:end,:);
    end
    %% Find RT



    sel                 = find(strcmp(strings(:,2),settings.SubjectID));
    sel                 = sel(1:numel(onsets));
    RT                  = numbers(sel,11);
    trls_enc            = strcmp(strings(sel,4),'encoding');
    trls_ret            = strcmp(strings(sel,4),'retrieval');



    RT                          = numbers(sel,11);

    RT(RT==-1 & trls_enc==1)    = 3; % -1 no press in time - set to 3s at encoding and 5s at retrieval
    RT(RT==-1 & trls_ret==1)    = 5; % -1 no press in time - set to 5s at encoding and 5s at retrieval

    %% Create a matrix with onsets for retrieval, ending at RT minus 250ms (max ripple = 500)/2 = 250ms.


    onsetmat_ret = zeros(numel(onsets),5001);
    for itrl=1:numel(onsets)
        to_add = onsets(1,itrl):onsets(1,itrl)+((RT(itrl,1)*data.fsample))-(time_to_excl_RT*data.fsample);
        onsetmat_ret(itrl,1:size(to_add,2)) = to_add;
        clear to_add
    end

    %% Extract ripples

    mat = zeros(numel(alldat),numel(data.time{1}));

    for ichannel = 1:numel(alldat)
        evs  = alldat{ichannel}.evtIndiv.maxTime;
        mat(ichannel,evs) = 1;
        clear evs
    end


    %% When ripples occured per channel


    for ichannel = 1:numel(alldat)
        tmp = mat(ichannel,:);
        allripples_samples_chan(ichannel,1:size(find(tmp),2)) = find(tmp);
        clear allripples tmp
    end
    clear tmp


    %% All relevant information about each trial.

    OverallTrial    = numbers(sel,1);
    OverallTrial    = OverallTrial-(OverallTrial(1,1)-1);
    SubjectID_all   = strings(sel,2);
    ExpPhase        = strings(sel,4);
    TrialNumber     = numbers(sel,5);
    BlockType       = strings(sel,6);
    Subcat          = strings(sel,7);
    Word            = strings(sel,8);
    OldNew          = strings(sel,9);
    Response        = strings(sel,10);
    Memory          = strings(sel,12);
    RunNumber       = numbers(sel,3);

    sel_ret         = ismember(ExpPhase,'retrieval');

    OverallTrial_ret= OverallTrial(sel_ret,1);

    %% Reduce matrices to only retrieval trials
    onsets_ret      = onsets(:,sel_ret);
    onsetmat_ret    = onsetmat_ret(sel_ret,:);

    %% Time stamps and duration of ripples
    rippletime_ret_tmp              = [];
    rippletime_ret_tmp_with_onset   = [];
    rippledur_ret_tmp               = [];

    rippletime_ret              = cell(numel(alldat),1);
    rippletime_ret_with_onset   = cell(numel(alldat),1);
    rippledur_ret               = cell(numel(alldat),1);


    for ichannel = 1:numel(alldat)

        fprintf('channel %01d / %01d\n',ichannel,numel(alldat));

        tmp = allripples_samples_chan(ichannel,:); % All timestamps per channel
        tmp(tmp==0)=[];
        for nucol = 1:size(tmp,2)
            time_rip = tmp(1,nucol); % single timestamp
            for itrl=1:numel(onsets_ret)
                [~,b] = find(onsetmat_ret(itrl,:)==time_rip); % Find trial of timestamp
                if isempty(b)
                    time_stamp  = NaN; % If not that trial, NaN it.
                    dur         = NaN;
                else
                    [~,tmp_b]   = find(alldat{ichannel}.evtIndiv.maxTime == time_rip); % If a trial, find the idx of that ripple in the data.
                    time_stamp  = alldat{ichannel}.evtIndiv.maxTime(1,tmp_b); % Add timestamp of ripple.
                    dur         = alldat{ichannel}.evtIndiv.duration(1,tmp_b); % Use idx to find the duration of that particular ripple
                end
                %
                [rippletime_ret_tmp(itrl,:)]                = time_stamp-onsets_ret(1,itrl); % Subtract onset
                [rippletime_ret_tmp_with_onset(itrl,:)]     = time_stamp+data.sampleinfo(1); % Leave onset in and add sampleinfo
                [rippledur_ret_tmp(itrl,:)]                 = dur; % Duration
                clear dur b time_stamp
            end
            rippletime_ret_col(:,nucol)                 = rippletime_ret_tmp;
            rippletime_ret_with_onset_col(:,nucol)      = rippletime_ret_tmp_with_onset;
            rippledur_ret_col(:,nucol)                  = rippledur_ret_tmp;

            clear rippletime_ret_tmp rippledur_ret_tmp...
                time_rip rippletime_ret_tmp_with_onset
        end
        clear tmp
        rippletime_with_onset_col   = [];
        rippletime_col              = [];
        rippledur_col               = [];

        for nurows = 1:size(rippletime_ret_col,1) % Reduce matrix
            [~,b] = find(rippletime_ret_col(nurows,:)>0);
            if isempty(b)
                [rippletime_col(nurows,1)]              = NaN;
                [rippletime_with_onset_col(nurows,1)]   = NaN;
                [rippledur_col(nurows,1)]               = NaN;

            else
                [rippletime_col(nurows,1:size(b,2))]                = rippletime_ret_col(nurows,b);
                [rippletime_with_onset_col(nurows,1:size(b,2))]     = rippletime_ret_with_onset_col(nurows,b);
                [rippledur_col(nurows,1:size(b,2))]                 = rippledur_ret_col(nurows,b);
            end
            clear b
        end
        clear rippletime_ret_col rippledur_ret_col rippleamp_ret_col ripplefreq_ret_col rippletime_ret_with_onset_col
        rippletime_col(rippletime_col==0)                       = NaN;
        rippletime_ret{ichannel,:}                              = [rippletime_col,OverallTrial_ret];
        rippletime_col                                          = [];
        rippletime_with_onset_col(rippletime_with_onset_col==0) = NaN;
        rippletime_ret_with_onset{ichannel,:}                   = [rippletime_with_onset_col,OverallTrial_ret];
        rippletime_with_onset_col                               = [];
        rippledur_col(rippledur_col==0)                         = NaN;
        rippledur_ret{ichannel,:}                               = [rippledur_col,OverallTrial_ret];
        rippledur_col                                           = [];
        clear tmp

    end


    %% for each column, make a mXn matrix with time of ripple in relation
    % to onset, time of ripple, duration of ripple and trial number.

    tmp_and_trial_all_col   = [];
    ripples_trials          = [];
    tmp_and_trial_all_chan  = [];


    for ichannel = 1:numel(alldat)
        tmp                 = rippletime_ret_with_onset{ichannel,:};
        tmp_dur             = rippledur_ret{ichannel,:};
        tmp_without_onset   = rippletime_ret{ichannel,:};

        for nucol = 1:size(tmp,2)-1
            tmp_col               = tmp(:,nucol);
            tmp_dur_col           = tmp_dur(:,nucol);
            tmp_dup_col           = tmp_without_onset(:,nucol);
            tmp_and_trial         = [tmp_dup_col,tmp_col,tmp_dur_col,tmp(:,end)];
            tmp_and_trial_all_col = [tmp_and_trial_all_col;tmp_and_trial];
            clear tmp_col tmp_and_trial tmp_dur_col
        end

        tmp_and_trial_all_col_duplicate = [];
        ripples_trials                  = [];

        for nucol = 1:size(rippletime_ret,2)-1
            tmp_col                         = rippletime_ret(:,nucol);
            tmp_and_trial                   = [tmp_col,rippletime_ret(:,end)];
            tmp_and_trial_all_col_duplicate = [tmp_and_trial_all_col_duplicate;tmp_and_trial];
            clear tmp_col tmp_and_trial
        end

        % At this point I can add all channels togehter, as they are not important
        % anymore.
        tmp_and_trial_all_chan  = [tmp_and_trial_all_chan;tmp_and_trial_all_col];
        tmp_and_trial_all_col   = [];
    end

    %% Delete non-ripple elements and make a matrix with the remaining ripples.

    [tmp,~]                         = find(tmp_and_trial_all_chan(:,2)>0);
    rippletime_ret_sorted           = tmp_and_trial_all_chan(tmp,:);
    [~,tmpsort]                     = sort(rippletime_ret_sorted(:,2),'ascend');
    rippletime_ret_sorted           = rippletime_ret_sorted(tmpsort,:);

    %% Create TRL

    trl                     = [];
    pretrig                 = round(6 * data.fsample);
    posttrig                = round(4 * data.fsample);

    for i = 1:size(rippletime_ret_sorted,1)
        offset              = -pretrig;
        trlbegin            = rippletime_ret_sorted(i,2) - pretrig;
        trlend              = rippletime_ret_sorted(i,2) + posttrig;
        time_ripple         = rippletime_ret_sorted(i,1);
        duration_ripple     = rippletime_ret_sorted(i,3);
        trial_no            = rippletime_ret_sorted(i,4);
        newtrl              = [trlbegin trlend offset trial_no time_ripple duration_ripple];
        trl                 = [trl; newtrl;];
    end

    %% Find trial number and add all the relevant data for the trialinfo
    clear tmp tmp_
    for nucol = 1:size(trl,1)
        tmp = trl(nucol,4);
        tmp_.trialinfo{nucol,1}.Cnt             = OverallTrial(tmp,1);
        tmp_.trialinfo{nucol,1}.SubjectID       = SubjectID_all(tmp,1);
        tmp_.trialinfo{nucol,1}.RunNumber       = RunNumber(tmp,1);
        tmp_.trialinfo{nucol,1}.ExpPhase        = ExpPhase(tmp,1);
        tmp_.trialinfo{nucol,1}.TrialNumber     = TrialNumber(tmp,1);
        tmp_.trialinfo{nucol,1}.BlockType       = BlockType(tmp,1);
        tmp_.trialinfo{nucol,1}.Source          = Subcat(tmp,1);
        tmp_.trialinfo{nucol,1}.Word            = Word(tmp,1);
        tmp_.trialinfo{nucol,1}.OldNew          = OldNew(tmp,1);
        tmp_.trialinfo{nucol,1}.Response        = Response(tmp,1);
        tmp_.trialinfo{nucol,1}.RT              = RT(tmp,1);
        tmp_.trialinfo{nucol,1}.Memory          = Memory(tmp,1);
        tmp_.trialinfo{nucol,1}.ripple_time     = trl(nucol,5);
        tmp_.trialinfo{nucol,1}.duration        = trl(nucol,6);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% Test the realignment %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Make sure the ripples are at time zero in hippocampal channels


    load('/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/ripple_project_publication_for_replication/templates/channels_hipp_ripples.mat')

    channels = channels_hipp_ripples(isubject,:);
    channels = channels(~cellfun(@isempty,channels));


    cfg                 = [];
    cfg.channel         = channels;
    data_hipp           = ft_selectdata(cfg, data);

    cfg                 = [];
    cfg.trl             = trl;
    data_realigned      = ft_redefinetrial(cfg,data_hipp);

    data_realigned.trialinfo = [];
    data_realigned.trialinfo = tmp_.trialinfo;

    cfg                 = [];
    cfg.detrend         = 'yes';
    cfg.demean          = 'yes';
    data_realigned                 = ft_preprocessing(cfg,data_realigned);

    cfg             = [];
    cfg.keeptrials  = 'yes';
    GA              = ft_timelockanalysis(cfg,data_realigned);


    m = squeeze(nanmean(nanmean(squeeze(GA.trial))));
    s = squeeze(std(nanmean(squeeze(GA.trial))))/sqrt(size(GA.trial,1));

    mean_all(isubject,:) = m;

    TOI             = [-6 4];
    toi_idx         = TOI(1):1/round(data.fsample):TOI(2);


    figure;
    boundedline(toi_idx,m,s,'k','alpha')
    axis tight
    hline(0);
    vline(0);
    xlabel('Time (sec)')
    ylabel('\muV')
    %

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Realign trials and save
    % Exclude artifacts and interpolate the NaNs with 3 before and after.
    for itrial= 1:numel(data.artifact)
        itrial
        tmp = data;
        tmp.trial = [];
        tmp.trial = data.trial{1,1}(itrial,:);
        tmp.artifact = data.artifact{1,itrial};
        for iarti = 1:size(tmp.artifact,1)
            tmp_arti = tmp.artifact(iarti,:);
            tmp.trial(1,tmp_arti(1):tmp_arti(2)) = NaN;
            clear tmp_arti
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
        clear tmp
    end

    cfg               = [];
    cfg.trl           = trl;
    data_realigned    = ft_redefinetrial(cfg,data);

    data_realigned.trialinfo = [];
    data_realigned.trialinfo = tmp_.trialinfo;


    %     save([settings.save_dir,'data_realigned_',num2str(settings.subjects(isubject,:))],'data_realigned','-v7.3')


    clearvars -except isubject settings mean_all
end





%% plot

fs = 1000;

TOI = -6:1/fs:4;
TOI_idx = nearest(TOI,-1):nearest(TOI,1);


m = squeeze(nanmean(mean_all(:,TOI_idx)));
s = squeeze(std(mean_all(:,TOI_idx)))/sqrt(size(mean_all,1));

TOI = TOI(TOI_idx);
figure;
boundedline(TOI,m,s,'k','alpha')
axis tight
hline(0);
vline(0);
xlabel('Time (sec)')
ylabel('\muV')