%% realign based on trial onset
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear

%% path settings
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects            = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW','09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

settings.data_dir      = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];
settings.save_dir      = [settings.base_path_castle,'output_data/time_locked/'];

for isubject = 1:size(settings.subjects,1)


    fprintf('time-locking subject %01d/%02d\n',isubject,size(settings.subjects,1));

    settings.SubjectID   = num2str(settings.SubjectIDs(isubject,:));

    %% LOAD


    data_in                  = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',num2str(settings.subjects(isubject,:))],'data','onsets_session');
    data                     = data_in.data;
    onsets                   = data_in.onsets_session;
    data_in                  = [];



    %% Exclude artifacts and interpolate the NaNs with 3 before and after.
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


    %% Create TRL
    trl           = [];
    pretrig       = round(1 * data.fsample);
    posttrig      = round(5 * data.fsample);

    for i = 1:numel(onsets)
        offset    = -pretrig;
        trlbegin  = onsets(i) - pretrig;
        trlend    = onsets(i) + posttrig;
        newtrl    = [trlbegin trlend offset];
        trl       = [trl; newtrl];
    end

    cfg     = [];
    cfg.trl = trl;
    data    = ft_redefinetrial(cfg,data);
    clear trl

    cfg.demean  = 'no';
    data        = ft_preprocessing(cfg,data);

    %% Get info from excel sheet
    [numbers,strings] = xlsread([settings.base_path_castle,'well01_behavior_all.xls']);

    strings = strings(2:end,:);
    if isnan(numbers(1,1))
        numbers         = numbers(2:end,:);
    end
    sel = find(strcmp(strings(:,2),settings.SubjectID));
    sel = sel(1:numel(onsets));

    for ientry=1:numel(onsets)

        data.trialinfo{ientry,1}.Cnt             = numbers(sel(ientry),1);
        data.trialinfo{ientry,1}.SubjectID       = strings(sel(ientry),2);
        data.trialinfo{ientry,1}.RunNumber       = numbers(sel(ientry),3);
        data.trialinfo{ientry,1}.ExpPhase        = strings(sel(ientry),4);
        data.trialinfo{ientry,1}.TrialNumber     = numbers(sel(ientry),5);
        data.trialinfo{ientry,1}.BlockType       = strings(sel(ientry),6);
        data.trialinfo{ientry,1}.Source          = strings(sel(ientry),7);
        data.trialinfo{ientry,1}.Word            = strings(sel(ientry),8);
        data.trialinfo{ientry,1}.OldNew          = strings(sel(ientry),9);
        data.trialinfo{ientry,1}.Response        = strings(sel(ientry),10);
        data.trialinfo{ientry,1}.RT              = numbers(sel(ientry),11);
        data.trialinfo{ientry,1}.Memory          = strings(sel(ientry),12);
        data.trialinfo{ientry,1}.PrevMemory      = strings(sel(ientry),13);
    end
    clear sel strings numbers

    cfg             = [];
    cfg.keeptrials  = 'yes';
    data_timelocked = ft_timelockanalysis(cfg,data);


    save([settings.save_dir,'data_timelocked_',num2str(settings.subjects(isubject,:))],'data_timelocked','-v7.3')

    clearvars -except settings isubject
end
%%