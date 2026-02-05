clear
%% automatic artifact detection (artifacts will be removed when realigning the trials in step 02b and 02c)
%                  Casper Kerren      [kerren@cbs.mpg.de]

settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subject            = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW','09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
settings.data_dir           = [settings.base_path_castle,'preprocessing/channel_removal_all_channels/common_trimmed_average/'];
settings.save_dir           = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];

% artifact detection values

settings.amp_thres           = 4.0;
settings.super_amp_thres     = 6.0;

settings.jump_thres          = 4.0;
settings.super_jump_thres    = 6.0;

settings.hifreq_thres        = 4.0;
settings.super_hifreq_thres  = 6.0;

settings.hpfreq              = 250;

% plotting
settings.plot_artifactsegments = 0; % flag, determining whether you want to plot all artifacts with a time window around them
settings.plot_allsegments      = 0; % flag, determining whether you want to plot raw data divided in segments with artifacts highlighted

%%
for nuPPP = 1:size(settings.subject,1)


    fprintf('processing subject %01d/%02d\n',nuPPP,size(settings.subject,1));


    data_in         = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_',num2str(settings.subject(nuPPP,:))],'data','onsets_session');
    data            = data_in.data;
    onsets_session  = data_in.onsets_session;
    around          = [5 10]; % how much time to include before first onset and after last onset


    %% restrict to relevant data (note that if you want the actual onset in ms after having cut around the relevant data onsets will have to be adjusted henceforth, eg onsets  = onsets_session - data.time{1}(1)*data.fsample;

    cfg         = [];
    cfg.latency = [onsets_session(1)/data.fsample-around(1) onsets_session(end)/data.fsample+around(2)];
    data        = ft_selectdata(cfg,data);

    data.sampleinfo = [data.time{1}(1)*data.fsample,data.time{1}(end)*data.fsample];
    %% read in all data for 1 channel
    for ichannel = 1:numel(data.label)

        cfg                     = [];
        cfg.channel             = data.label{ichannel};
        this_data               = ft_preprocessing(cfg,data);

        % If NaNs in recording, interpolate these. Later on,
        % I add these values as artifacts
        cfg             = [];
        cfg.prewindow   = 3;
        cfg.postwindow  = 3;
        this_data       = ft_interpolatenan(cfg,this_data);


        %% create high-pass filtered data
        data_hp     = ft_preproc_highpassfilter(this_data.trial{1},this_data.fsample,settings.hpfreq);

        %% ARTIFACT DETECTION

        %% -- amplitude
        trlmat                              = ft_preproc_detrend(this_data.trial{1});
        zscoremat                           = zscore(trlmat,[],2);
        trial_amplitude_artifact            = logical(abs(zscoremat) > settings.amp_thres);
        trial_amplitude_artifact_extreme    = logical(abs(zscoremat) > settings.super_amp_thres);
        zscoremat                           = [];

        %% -- gradient.
        jumps                       = [0 diff(this_data.trial{1},1,2)];
        zscoremat                   = zscore(jumps,[],2);
        trial_jump_artifact         = logical(abs(zscoremat) > settings.jump_thres);
        trial_jump_artifact_extreme = logical(abs(zscoremat) > settings.super_jump_thres);
        zscoremat                   = [];

        %% high frequency bursts
        zscoremat                       = zscore(data_hp,[],2);
        trial_hifreq_artifact           = logical(abs(zscoremat) > settings.hifreq_thres);
        trial_hifreq_artifact_extreme   = logical(abs(zscoremat) > settings.super_hifreq_thres);
        zscoremat                       = [];


        %% NaNs
        cfg                             = [];
        cfg.channel                     = data.label{ichannel};
        data_nan                        = ft_preprocessing(cfg,data);
        trial_NaNs_artifact             = isnan(data_nan.trial{1,1});

        %% combine all types of artifacts

        anyartifact = logical(...
            trial_amplitude_artifact & (trial_jump_artifact | trial_hifreq_artifact)    | ...
            trial_amplitude_artifact_extreme                                            | ...
            trial_jump_artifact_extreme                                                 | ...
            trial_hifreq_artifact_extreme                                               | ...
            trial_NaNs_artifact);

        %% mark artifact samples
        arttype  = anyartifact; % anyartifact flatline trial_amplitude_artifact_extreme trial_jump_artifact_extreme;

        artbeg   = find(diff([0 arttype]) ==  1);
        artend   = find(diff([arttype 0]) == -1);
        artifact = [artbeg(:) artend(:)];

        fprintf('total samples: %d, clean samples (unpadded): %d (%3.0f%%)\n',numel(arttype),sum(arttype==0), 100*sum(arttype==0)/numel(arttype));

        this_data.artifact = artifact;

        %% plot artifacts artifactwise
        if settings.plot_artifactsegments

            startingpoint   = 1;
            around          = 10000;
            nsegs           = round(numel(data.trial{1})/around);

            figure;
            for iart = startingpoint:1:size(artifact,1)

                dataplotidx = (artifact(iart,1)-around/2:artifact(iart,2)+around/2);
                dataplotidx = dataplotidx(dataplotidx > 0);

                dataline      = this_data.trial{1}(dataplotidx);
                dataline_hp   = data_hp(dataplotidx);

                dataplottime  = this_data.time{1}(dataplotidx);

                nanline    = nan(size(dataline));
                nanline_hp = nan(size(dataline));
                artidx     = artifact(iart,1):artifact(iart,2);
                artplotidx = find(dataplotidx >= artidx(1) & dataplotidx <= artidx(end));

                nanline(artplotidx)     = this_data.trial{1}(artidx);
                nanline_hp(artplotidx)  = data_hp(artidx);

                subplot(211)
                plot(dataplottime,dataline,'k');
                hold on
                plot(dataplottime,nanline,'r.-','linewidth',2)
                title(sprintf('seg %d/%d, art %d/%d',ceil(artifact(iart,1)/around),nsegs,iart,size(artifact,1)))
                axis tight

                subplot(212)
                plot(dataplottime,dataline_hp,'k');
                hold on
                plot(dataplottime,nanline_hp,'r.-','linewidth',2)
                title(sprintf('seg %d/%d, art %d/%d, HP-filtered',ceil(artifact(iart,1)/around),nsegs,iart,size(artifact,1)))
                axis tight

                pause()
                clf

            end
        end

        %% pad artifacts for visualization

        if settings.plot_allsegments

            refract        = 2000;
            artifactpadded = [artifact(:,1)-refract/2 artifact(:,2)+refract/2];

            thisline    = this_data.trial{1};
            thisline_hp = data_hp;

            nanline     = nan(size(thisline));
            nanline_hp  = nan(size(thisline));

            for iart=1:size(artifactpadded,1)
                fprintf('%d/%d\n',iart,size(artifactpadded,1))
                idx = artifactpadded(iart,1):artifactpadded(iart,2);
                idx = idx(idx > 0 & idx < numel(thisline));
                nanline(idx)    = thisline(idx);
                nanline_hp(idx) = thisline_hp(idx);
            end

        end
        %% plot padded artifacts segmentwise

        if settings.plot_allsegments
            segsize = 20000;
            nseg = numel(thisline)/segsize;

            figure
            for iseg=1:nseg

                idx = (iseg-1)*segsize+1 : iseg*segsize;

                subplot(211)
                plot(thisline(idx),'g');
                hold on
                plot(nanline(idx),'r-')
                title(sprintf('seg%d/%d',iseg,round(nseg)))
                set(gca,'tickdir','out')

                subplot(212)
                plot(thisline_hp(idx),'k');
                hold on
                plot(nanline_hp(idx),'r-')
                title(sprintf('seg%d/%d, HP-filtered',iseg,round(nseg)))
                set(gca,'tickdir','out')

                pause
                clf
            end
            close all
        end

        %% save data with artifact info

        nanline                             = [];
        thisline                            = [];
        anyartifact                         = [];
        trial_amplitude_artifact            = [];
        trial_amplitude_artifact_extreme    = [];
        trial_jump_artifact                 = [];
        trial_jump_artifact_extreme         = [];
        artifact                            = [];
        artifactpadded                      = [];
        artbeg                              = [];
        artend                              = [];
        data_hp                             = [];


        data.artifact{ichannel} = this_data.artifact;


    end


    data.artifactdef.gradient.jumppeak  = settings.jump_thres;
    data.artifactdef.gradient.extreme   = settings.super_jump_thres;
    data.artifactdef.hifreq.hipassfreq  = settings.hpfreq;
    data.artifactdef.hifreq.jumppeak    = settings.hifreq_thres;
    data.artifactdef.hifreq.extreme     = settings.super_hifreq_thres;
    data.artifactdef.amplitude.jumppeak = settings.amp_thres;
    data.artifactdef.amplitude.extreme  = settings.super_amp_thres;

    save([settings.save_dir,'eeg_session01_all_chan_nobadchan_cmntrim_artdet_',num2str(settings.subject(nuPPP,:))],'data','onsets_session','-v7.3')

    clearvars -except settings nuPPP
end
%%