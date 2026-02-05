%% detect ripples around reaction time
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

[~,ftpath]=ft_version;


%% path settings
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

settings.data_dir           = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];
settings.save_dir           = [settings.base_path_castle,'output_data/ripple_detection/'];

load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1;

settings.pre_win    = 1; % Window before trigger in secs
settings.post_win   = 1; % Window after trigger in secs
settings.smoothing  = .4;

settings.zscore_task = 1;


addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))

% for plotting
settings.fs         = 1000;
settings.conditions = {'ret_IO' 'ret_IS'};

%% start for loop
ripple_latency_correct_retrieval = cell(1,size(settings.subjects,1));
ripple_latency_incorrect_retrieval = cell(1,size(settings.subjects,1));

tic

numWorkers = 8; %

parpool('local', numWorkers);

parfor isubject = 1:size(settings.subjects,1)

    fprintf('processing subject %01d/%02d\n',isubject,size(settings.subjects,1));

    %% LOAD
    data_in             = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',settings.subjects(isubject,:)],'data','onsets_session');
    data                = data_in.data;
    onsets_session      = data_in.onsets_session;
    data_in             = [];

    %% load channels and restrict to channels that are in hippocampus, based on anatomy

    tmp = load([settings.base_path_castle,'ripple_project_publication_for_replication/templates/channels_hipp_ripples']);
    channels_hipp_ripples = tmp.channels_hipp_ripples;
    tmp = [];
    channels = channels_hipp_ripples(isubject,:);
    channels = channels(~cellfun(@isempty,channels));

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

    %% LOAD excel sheet for RTs

    onsets              = onsets_session - data.time{1}(1)*inData.fsample;
    [numbers,strings]   = xlsread([settings.base_path_castle,'/well01_behavior_all.xls']);
    strings             = strings(2:end,:);

    if isnan(numbers(1,1))
        numbers         = numbers(2:end,:);
    end


    %% delete co-occuring ripples

    if settings.remove_ripple_duplicates

        alldat = remove_ripple_duplicates(alldat);

    end


    %% Divide into conditions


    [numbers,strings] = xlsread([settings.base_path_castle,'/well01_behavior_all.xls']);

    strings     = strings(2:end,:);
    sel         = find(strcmp(strings(:,2),settings.SubjectIDs(isubject,:)));
    sel         = sel(1:numel(onsets));

    if isnan(numbers(1,1))
        numbers         = numbers(2:end,:);
    end

    ExpPhase        = strings(sel,4);
    Memory          = strings(sel,12);
    RT              = numbers(sel,11);

    trls_enc            = strcmp(strings(sel,4),'encoding');
    trls_ret            = strcmp(strings(sel,4),'retrieval');


    RT(RT==-1 & trls_enc==1)    = 3; % -1 no press in time - set to 3s at encoding

    if settings.full_enc_trial
        RT(trls_enc==1) = 3; % [optional] set all encoding to 3s
    end
    RT(RT==-1 & trls_ret==1)    = 5; % -1 no press in time - set to 5s at retrieval

    %% create onset matrix
    onsetmat = [onsets; onsets+(RT'.*data.fsample)]';

    %% create trialinfo
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


    sel_ret_IS = ismember(ExpPhase,'retrieval') & ismember(Memory,'SourceCorrect');
    sel_ret_IO = ismember(ExpPhase,'retrieval') & ismember(Memory,{'SourceIncorrect' 'SourceDunno' 'dunno'});

    % create trialinfo for later epoching
    trialinfo_both      = sel_ret_IS;
    trialinfo_both      = double(trialinfo_both);
    idx                 = find(sel_ret_IO);
    trialinfo_both(idx) = 2;
   

     %% Find ripples, zscore and smooth before epoching
    channel_tracker = [];
    data_writeout   = {};
    for ichannel = 1:numel(alldat)

    
        % find ripple events
        evs = alldat{ichannel}.evtIndiv.maxTime;

        first_onset =  min(onsets, [], 'all');
        final_onset =  max(onsets, [], 'all');

        % Define time of interest from first to last trial
        time_idx = first_onset-round(settings.pre_win*inData.fsample):final_onset+(RT(end)*inData.fsample)+round(settings.post_win*inData.fsample)+1;

        % Do hist count
        ripple_all = histcounts(evs,time_idx);

        % Remove one from time axis to account for hist count
        time_idx = time_idx(1:end-1);

        % Convert to Hz
        ripple_all = ripple_all .* inData.fsample;

        % Smooth data
        ripple_all = smoothdata(ripple_all,2,'movmean',settings.smoothing*inData.fsample);

        % chop into trials - account for RT (onset plus RT minus pre-window : onset plus RT plus post-window)
        ripple_trial = [];
        for itrial = 1:numel(onsets)
            time_tral = onsets(itrial) + (RT(itrial)*inData.fsample) - (settings.pre_win*inData.fsample) : onsets(itrial) + (RT(itrial)*inData.fsample) + (settings.pre_win*inData.fsample);
            ripple_trial(itrial,:) = ripple_all(time_idx>=time_tral(1) & time_idx < time_tral(end));
        end

        % only take retrieval trials

        idx_all = sel_ret_IO | sel_ret_IS;


        ripple_trial        = ripple_trial(idx_all,:);
        onsets_tmp          = onsets(idx_all);
        trialinfo_both_tmp  = trialinfo_both(idx_all);

        tmp         = [];
        time_tral   = -settings.pre_win:1/inData.fsample:settings.post_win;
        for itrial = 1:numel(onsets_tmp)

            tmp.trial{itrial}   = ripple_trial(itrial,:);
            tmp.time{itrial}    = time_tral(1:end-1);

        end
        tmp.trialinfo   = trialinfo_both_tmp;
        tmp.label       = {'channel'};

        %% Z-score [optional]
        if settings.zscore_task
            concat = [];
            for itrial = 1:numel(tmp.trial)
                concat = [concat tmp.trial{itrial}];
            end
            mu      = mean(concat);
            sigma   = std(concat);
            for itrial = 1:numel(tmp.trial)
                tmp.trial{itrial} = (tmp.trial{itrial} -mu) / sigma;
            end
        end

        %% Time-lock
        for icondition = 1:numel(settings.conditions)

            cfg                                 = [];
            cfg.keeptrials                      = 'no';
            cfg.trials                          = tmp.trialinfo == icondition;
            data_writeout{icondition}{ichannel} = ft_timelockanalysis(cfg,tmp);

        end

    end

    %% event means
    tmp_1 = [];
    tmp_2 = [];
     for ichannel = 1:numel(alldat)
        tmp_1 = [tmp_1;data_writeout{1, 1}{1, ichannel}.avg];
        tmp_2 = [tmp_2;data_writeout{1, 2}{1, ichannel}.avg];
      end

    ripple_latency_correct_retrieval{isubject}.rate      = mean(tmp_1);
    ripple_latency_incorrect_retrieval{isubject}.rate    = mean(tmp_2);

    ripple_latency_correct_retrieval{isubject}.time      = time_tral;


end
toc

delete(gcp);
return


ripple_rate_time_resolved_RT.retrieval.correct     = ripple_latency_correct_retrieval;
ripple_rate_time_resolved_RT.retrieval.incorrect   = ripple_latency_incorrect_retrieval;

save ripple_rate_time_resolved_RT ripple_rate_time_resolved_RT

%% Gather data and plot

ripple_rate_correct     = [];
ripple_rate_incorrect   = [];
for isubject = 1:size(settings.subjects,1)
    ripple_rate_correct(isubject,:)     = ripple_latency_correct_retrieval{1, isubject}.rate;
    ripple_rate_incorrect(isubject,:)   = ripple_latency_incorrect_retrieval{1, isubject}.rate;
end

plot(mean(ripple_rate_correct))
hold on
plot(mean(ripple_rate_incorrect))

%%
correct_incorrect = {};

for isubject = 1:size(settings.subjects,1)

    correct_incorrect{1}.label                      = {'Channels'};
    correct_incorrect{1}.time                       = ripple_latency_correct_retrieval{1, 1}.time(1:end-1);

    correct_incorrect{1}.individual(isubject,1,:)   = ripple_rate_correct(isubject,:);
    correct_incorrect{1}.dimord                     = 'subj_chan_time';



end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);


correct_incorrect{2}                = correct_incorrect{1};

for isubject =1:size(settings.subjects,1)

    correct_incorrect{2}.individual(isubject,1,:)           = ripple_rate_incorrect(isubject,:);

end



cfg                     = [];
cfg.latency             = [-1 1];

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

d = squeeze(correct_incorrect{1}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

figure;
boundedline(correct_incorrect{1, 1}.time  ,m,s,'b');
plot(correct_incorrect{1, 1}.time  ,m,'k','linewidth',2);
hold on

d = squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(correct_incorrect{1, 1}.time  ,m,s,'r');
plot(correct_incorrect{1, 1}.time  ,m,'k','linewidth',2);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time,sigline,'r','linewidth',4);

set(gca,'FontSize',16,'FontName','Arial')
xlabel('Time after cue (sec)')
ylabel('ripple rate (hz)')
set(gca,'TickDir','out')
title('difference ripple rate retrieval')

axis tight
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])

