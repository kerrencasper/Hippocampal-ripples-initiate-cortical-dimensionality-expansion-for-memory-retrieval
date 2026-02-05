%% detect ripples
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

% [~,ftpath]=ft_version;

%% path settings
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects            = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW','09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');


settings.data_dir           = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];
settings.save_dir           = [settings.base_path_castle,'output_data/ripple_detection/'];

load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1;

addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])

% for plotting
settings.fs                     = 1000;
settings.conditions             = {'enc_IO' 'enc_IS' 'ret_IO' 'ret_IS'};

%% start for loop

ripples_condition = cell(1,size(settings.subjects,1));

tic
parfor isubject = 1:size(settings.subjects,1)


    fprintf('processing subject %01d/%02d\n',isubject,size(settings.subjects,1));


    %% LOAD
    data_in             = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',settings.subjects(isubject,:)],'data','onsets_session');
    data                = data_in.data;
    onsets_session      = data_in.onsets_session;
    data_in             = [];

    %% load channels and restrict to channels that are in hippocampus, based on anatomy

    tmp = load([settings.base_path_castle,'/ripple_project_publication_for_replication/templates/channels_hipp_ripples']);
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
    %% Find RT

    sel                 = find(strcmp(strings(:,2),settings.SubjectIDs(isubject,:)));
    sel                 = sel(1:numel(onsets));
    trls_enc            = strcmp(strings(sel,4),'encoding');
    trls_ret            = strcmp(strings(sel,4),'retrieval');

    RT                          = numbers(sel,11);
    RT(RT==-1 & trls_enc==1)    = 3; % -1 no press in time - set to 3s at encoding

    if settings.full_enc_trial
        RT(trls_enc==1) = 3; % [optional] set all encoding to 3s
    end
    RT(RT==-1 & trls_ret==1)    = 5; % -1 no press in time - set to 5s at retrieval

    %% create onset matrices

    onsetmat = [onsets; onsets+(RT'.*data.fsample)]';

      %% delete co-occuring ripples

    if settings.remove_ripple_duplicates

        alldat = remove_ripple_duplicates(alldat);

    end
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

    %% CHARACTERISTICS OF RIPPLES per trial and channel

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
    ripples_trial_channel.channeldensity = mean(ripples_trial_channel.number>0,2); % poportion of contacts showing a ripple

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

    %% event means

    sel_enc_IO = ismember(ExpPhase,'encoding')  & ismember(Memory,{'SourceIncorrect' 'SourceDunno'});
    sel_enc_IS = ismember(ExpPhase,'encoding')  & ismember(Memory,'SourceCorrect');
    sel_ret_IO = ismember(ExpPhase,'retrieval') & ismember(Memory,{'SourceIncorrect' 'SourceDunno' 'dunno'});
    sel_ret_IS = ismember(ExpPhase,'retrieval') & ismember(Memory,'SourceCorrect');


    ripples_condition{isubject}.number        =  [...
        mean(mean(ripples_trial_channel.number(sel_enc_IO,:),2))
        mean(mean(ripples_trial_channel.number(sel_enc_IS,:),2))
        mean(mean(ripples_trial_channel.number(sel_ret_IO,:),2))
        mean(mean(ripples_trial_channel.number(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.density        =  [...
        mean(mean(ripples_trial_channel.density(sel_enc_IO,:),2))
        mean(mean(ripples_trial_channel.density(sel_enc_IS,:),2))
        mean(mean(ripples_trial_channel.density(sel_ret_IO,:),2))
        mean(mean(ripples_trial_channel.density(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.duration        =  [...
        nanmean(nanmean(ripples_trial_channel.duration(sel_enc_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.duration(sel_enc_IS,:),2))
        nanmean(nanmean(ripples_trial_channel.duration(sel_ret_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.duration(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.frequency        =  [...
        nanmean(nanmean(ripples_trial_channel.frequency(sel_enc_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.frequency(sel_enc_IS,:),2))
        nanmean(nanmean(ripples_trial_channel.frequency(sel_ret_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.frequency(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.amplitude        =  [...
        nanmean(nanmean(ripples_trial_channel.amplitude(sel_enc_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.amplitude(sel_enc_IS,:),2))
        nanmean(nanmean(ripples_trial_channel.amplitude(sel_ret_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.amplitude(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.activity        =  [...
        nanmean(nanmean(ripples_trial_channel.activity(sel_enc_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.activity(sel_enc_IS,:),2))
        nanmean(nanmean(ripples_trial_channel.activity(sel_ret_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.activity(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.latency_first        =  [...
        nanmean(nanmean(ripples_trial_channel.latency_first(sel_enc_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.latency_first(sel_enc_IS,:),2))
        nanmean(nanmean(ripples_trial_channel.latency_first(sel_ret_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.latency_first(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.latency_all        =  [...
        nanmean(nanmean(ripples_trial_channel.latency_all(sel_enc_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.latency_all(sel_enc_IS,:),2))
        nanmean(nanmean(ripples_trial_channel.latency_all(sel_ret_IO,:),2))
        nanmean(nanmean(ripples_trial_channel.latency_all(sel_ret_IS,:),2))]';

    ripples_condition{isubject}.channeldensity        =  [...
        nanmean(ripples_trial_channel.channeldensity(sel_enc_IO,:))
        nanmean(ripples_trial_channel.channeldensity(sel_enc_IS,:))
        nanmean(ripples_trial_channel.channeldensity(sel_ret_IO,:))
        nanmean(ripples_trial_channel.channeldensity(sel_ret_IS,:))]';

%     save([settings.save_dir,'eeg_session01_all_chan_nobadchan_cmntrim_ripple_',num2str(settings.subjects(isubject,:))],'alldat','ripples_condition', 'ripples_trial_channel', 'onsets','onsets_session','-v7.3')
%     clearvars -except settings isubject
end
toc
return

%% Gather data and plot

vars = {'number' 'density' 'duration' 'frequency' 'amplitude' 'activity' 'latency_first' 'latency_all' 'channeldensity'};


all_number          = [];
all_density         = [];
all_duration        = [];
all_frequency       = [];
all_amplitude       = [];
all_activity        = [];
all_latency_first   = [];
all_latency_all     = [];
all_channeldensity  = [];

for isubject=1:size(settings.subjects,1)

    all_number(isubject,:)          = ripples_condition{isubject}.number;
    all_density(isubject,:)         = ripples_condition{isubject}.density;
    all_duration(isubject,:)        = ripples_condition{isubject}.duration;
    all_frequency(isubject,:)       = ripples_condition{isubject}.frequency;
    all_amplitude(isubject,:)       = ripples_condition{isubject}.amplitude;
    all_activity(isubject,:)        = ripples_condition{isubject}.activity;
    all_latency_first(isubject,:)   = ripples_condition{isubject}.latency_first;
    all_latency_all(isubject,:)     = ripples_condition{isubject}.latency_all;
    all_channeldensity(isubject,:)  = ripples_condition{isubject}.channeldensity;

end


figure;

my_gray = settings.colour_scheme(8,:);

for ivar=1:numel(vars)

    d = eval(['all_' vars{ivar}]);
    m = nanmean(d);
    s = nanstd(d)./sqrt(size(d,1));

    bar(1:numel(m),m,'facecolor',my_gray)
    hold on
    errorbar(1:numel(m),m,s,'color',my_gray,'linestyle','none','linewidth',2)
    [~,pval] = ttest(d(:,1),d(:,2));
    %     [pval]   = signrank(d(:,1),d(:,2));
    text(1.5,mean(m(1:2)),sprintf('%0.3f',pval))
    [~,pval] = ttest(d(:,3),d(:,4));
    %     [pval]   = signrank(d(:,3),d(:,4));
    text(3.5,mean(m(3:4)),sprintf('%0.3f',pval))
    set(gca,'xtick',1:4,'xticklabel',strrep(settings.conditions,'_',' '))
    xlim([0 5])
    title(vars{ivar},'interpreter','none')

    pause
    clf
end
close all



[group1, group2] = deal(all_density(:,4), all_density(:,3));

[~, p_val, ~, stats]    = ttest(group1, group2);
t_stat                  = stats.tstat;


n1 = length(group1);
n2 = length(group2);


s1 = std(group1);
s2 = std(group2);


df = ((s1^2 / n1 + s2^2 / n2)^2) / ((s1^2 / n1)^2 / (n1 - 1) + (s2^2 / n2)^2 / (n2 - 1));

% Calculate Cohen's d
cohens_d = 2 * t_stat / sqrt(df);

fprintf('T-statistic: %.4f\n', t_stat);
fprintf('P-value: %.4f\n', p_val);
fprintf('Cohen''s d: %.4f\n', cohens_d);


%%