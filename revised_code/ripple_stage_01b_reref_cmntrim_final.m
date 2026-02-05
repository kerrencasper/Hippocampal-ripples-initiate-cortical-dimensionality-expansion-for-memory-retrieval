%% re-reference based on common mean
%                  Casper Kerren      [kerren@cbs.mpg.de]
clear
%% stage 1: focus on HIPP contacts and do common (trimmed) average rereferencing

addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects            = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW','09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
settings.data_dir           = [settings.base_path_castle,'preprocessing/channel_removal_all_channels/'];
settings.trimaway           = 20;
settings.save_dir           = [settings.base_path_castle,'preprocessing/channel_removal_all_channels/common_trimmed_average/'];
settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};

for nuPPP = 1:size(settings.subject,1)

    fprintf('processing subject %01d/%02d\n',nuPPP,size(settings.subject,1));
    load([settings.data_dir,'eeg_session01_all_chan_nobadchan_',num2str(settings.subject(nuPPP,:))])

    %% create common (trimmed) average

    cfg         = [];
    cfg.channel = find(~ismember(data.label,settings.scalp_channels));
    data        = ft_selectdata(cfg,data);


    trldat          = data.trial{1};
    ref             = trimmean(trldat,settings.trimaway);
    ref             = repmat(ref,[numel(data.label) 1]);
    data.trial{1}   = trldat - ref;

    % browse data

    cfg             = [];
    cfg.viewmode    = 'vertical';
    cfg.blocksize   = 30;
    ft_databrowser(cfg,data);

    pause
    %% save

    save([settings.save_dir,'eeg_session01_all_chan_nobadchan_cmntrim_',num2str(settings.subject(nuPPP,:))],'data','onsets_session','-v7.3')


    clearvars -except settings nuPPP
end
%%

