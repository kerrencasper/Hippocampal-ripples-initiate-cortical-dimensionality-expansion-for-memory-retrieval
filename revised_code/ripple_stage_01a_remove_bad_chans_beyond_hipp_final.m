%% remove physically defect channels
% (important for later medianfilter, classification etc)
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear

addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses/general_scripts_matlab/fieldtrip-20230422')
ft_defaults
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subject            = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('s01_CF', 's02_JM', 's03_SO', 's06_AH','s07_FC', 's08_HW', 's09_AM', 's10_MH','s11_FS', 's12_AS', 's13_CB', 's14_KK');
settings.data_dir           = [settings.base_path_castle,'subjects/'];
settings.save_dir           = [settings.base_path_castle,'preprocessing/channel_removal_all_channels/'];

% filter settings
settings.lpfreq         = 250;
settings.hpfreq         = .1;
settings.hpfiltord      = 4;
settings.bsfreq         = [49 51;99 101;149 151;199 201];
settings.filttype       = 'but';


for nuPPP = 1:size(settings.subject,1)
    
    
    fprintf('processing subject %01d/%02d\n',nuPPP,size(settings.subject,1));
    load([settings.data_dir,num2str(settings.SubjectIDs(nuPPP,:)),'/eeg.mat'])
    
    
    %%
    cfg             = [];
    data_all        = ft_preprocessing(cfg,eeg);
    
    %% inspect session-wise
    
    this_session = 1;
    
    %% select session
    
    cfg             = [];
    cfg.trials      = this_session;
    data_in         = ft_selectdata(cfg,data_all);
    
    clear data_all
    

    %% browse data (High, low pass and band stop filter the data)
    
    cfg                         = [];
    cfg.viewmode                = 'vertical';
    cfg.blocksize               = 30;
    cfg.preproc.bsfilter        = 'yes';
    cfg.preproc.bsfreq          = settings.bsfreq;
    cfg.preproc.lpfilter        = 'yes';
    cfg.preproc.lpfreq          = settings.lpfreq;
    cfg.preproc.lpfilttype      = settings.filttype;
    cfg.preproc.hpfilter        = 'yes';
    cfg.preproc.hpfreq          = settings.hpfreq;
    cfg.preproc.hpfilttype      = settings.filttype;
    cfg.preproc.hpfiltord       = settings.hpfiltord;
    cfg.preproc.instabilityfix  = 'split';
    ft_databrowser(cfg,data_in)
    
    pause
    
    %% subject MH has a lot of noise from 2700000:3000000 - NaN these entries
    if nuPPP == 8
        datain.trial{1}(:,2700000:3000000) = NaN;
    end
    
    %% do the FT artifactrejection - OBS! only for channels! (High, low pass and band stop filter the data)
    
    display_max   = [];
    artrej_data   = data_in;
    rej_meth      = menu('select way of artifact rejection:','none','summary_hpfilt','summary','channel','trial');
    
    while rej_meth > 1
        
        cfg             = [];
        cfg.channel     = 'all';
        if numel(eeg.trial) > 1
            cfg.latency     = [onsets{this_session}(1)/1000-5 onsets{this_session}(end)/1000+10];
        else
            cfg.latency     = [onsets(1)/1000-5 onsets(end)/1000+10];
        end
        cfg.alim                    = display_max;
        cfg.keeptrial               = 'yes';
        cfg.preproc.lnfilter        = 'no';
        cfg.preproc.demean          = 'yes';
        cfg.preproc.detrend         = 'no';
        cfg.preproc.bsfilter        = 'yes';
        cfg.preproc.bsfreq          = settings.bsfreq;
        cfg.preproc.lpfilter        = 'yes';
        cfg.preproc.lpfreq          = settings.lpfreq;
        cfg.preproc.lpfilttype      = settings.filttype;
        cfg.preproc.hpfilter        = 'yes';
        cfg.preproc.hpfreq          = settings.hpfreq;
        cfg.preproc.hpfilttype      = settings.filttype;
        cfg.preproc.hpfiltord       = settings.hpfiltord;
        cfg.preproc.instabilityfix  = 'split';
        switch rej_meth
            case 2
                cfg.method           = 'summary';
                cfg.metric           = 'maxabs';
                cfg.preproc.hpfilter = 'yes';
                cfg.preproc.hpfreq   = 150;
            case 3
                cfg.method = 'summary';
                cfg.metric = 'maxabs';
            case 4
                cfg.method = 'channel';
            case 5
                cfg.method = 'trial';
        end
        cfg.keepchannel = 'no';
        
        artrej_data   = ft_rejectvisual(cfg, artrej_data);
        
        rej_meth      = menu('select way of artifact rejection:','none','summary_hpfilt','summary','channel','trial');
        
    end
    
    %% find channels in non-filtered data.
    
    artrej_data.rej_chan = setdiff(artrej_data.cfg.channel,artrej_data.label);

    %% save
    
    data = artrej_data;
    if numel(eeg.trial) > 1
        onsets_session  = onsets{this_session};
    else
        onsets_session  = onsets;
    end
    
    save(sprintf([settings.save_dir,'eeg_session%01d_all_chan_nobadchan_',num2str(settings.subject(nuPPP,:))],this_session),'data','onsets_session','-v7.3');
    
    clearvars -except nuPPP settings
end
%%