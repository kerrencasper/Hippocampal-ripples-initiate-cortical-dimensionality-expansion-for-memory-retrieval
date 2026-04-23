clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults
%% automatic artifact detection (artifacts will be removed when realigning the trials in step 02b and 02c)
%                  Bernhard Staresina [bernhard.staresina@psy.ox.ac.uk]
%                  Casper Kerren      [kerren@cbs.mpg.de]

paths = config_paths();

settings = [];
settings.base_path_castle = paths.base_path;

settings.save_dir           = paths.save_dir;
settings.data_dir_channels  = paths.channels_dir;
settings.anatomy_dir        = paths.anatomy_dir;
settings.AAL_dir            = paths.AAL_dir;
settings.SPM_dir            = paths.SPM_dir;

settings.data_dir           = [settings.base_path_castle,'preprocessing/channel_removal_all_channels/common_trimmed_average/'];
settings.data_dir_art       = [settings.base_path_castle,'preprocessing/artifact_rejected_data/'];


load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;


settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};
settings.subject            = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW','09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
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
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur
settings.do_surrogates              = 0; % switch time of ripples between trials, 1 == for all trials, 2 == for correct trials only

%% decoding settings

% data preprocessing
settings.ori_fsample            = 1000; % original sample frequency
settings.do_resample            = 100; % [] or sample frequency

% baseline and zscoring settings
settings.zscore_data4class      = 1;
settings.bs_correct             = 1;
settings.bs_period              = [-.2 0]; % [-.5 -.1]
settings.bs_trim                = 0; % can be 0. amount of % to trim away when calculating the baseling


% smoothing options
settings.do_smoothdata          = 1; % use matlabs smoothdata function for running average
settings.smooth_win             = .200; % .100, [] running average time window in seconds;

% time of interest
settings.TOI_train              = [-.5 3];
settings.timesteps_train        = ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample); % in s. if you want it to take less sample points multiply [e.g., ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample)*2

settings.TOI_test               = [-1.2 1.2]; % time around ripple (I take this time window to get a proper estimate around the edges too. Only look at -1 to 1 later.
settings.timesteps_test         = ((settings.ori_fsample/settings.do_resample)/settings.ori_fsample); % in s

settings.classifier             = 'lda';
settings.metric                 = 'auc';

%% channel settings

settings.channel_sel            = 1; % 1 exclude hippo, 2 only hippo, 3 all channels

%% settings pca

settings.smooth_before_dim      = 1; % smooth Nans before doing dim reduction
settings.nu_time_points         = 60; % time of sliding window in ms
settings.prc_overlap            = .9; % percentage overlap sliding window
settings.decode_components      = 0; % Decode the PCA components I picked.
%% start for loop

timeaxis = settings.TOI_test(1):settings.timesteps_test:settings.TOI_test(2);
freqaxis = settings.TOI_train(1):settings.timesteps_train:settings.TOI_train(2);

perf    = cell(1,numel(subjects));
channs  = cell(1,numel(subjects));

RT_all_subj_correct     = cell(1,numel(subjects));
RT_all_subj_incorrect   = cell(1,numel(subjects));
ripple_time_correct     = cell(1,numel(subjects));
ripple_time_incorrect   = cell(1,numel(subjects));
% artifact detection values


% plotting
settings.plot_STA      = 0; % If plot spike-triggered average

numWorkers = 8; %

parpool('local', numWorkers);


%%
parfor isubject = 1:size(settings.subject,1)

    %% find spikes
 
    data_in         = load([settings.data_dir,'eeg_session01_all_chan_nobadchan_cmntrim_',num2str(settings.subject(isubject,:))],'data','onsets_session');
    data            = data_in.data;
    onsets_session  = data_in.onsets_session;
    around          = [5 10]; % how much time to include before first onset and after last onset


    %% restrict to relevant data (note that if you want the actual onset in ms after having cut around the relevant data onsets will have to be adjusted henceforth, eg onsets  = onsets_session - data.time{1}(1)*data.fsample;

    cfg         = [];
    cfg.latency = [onsets_session(1)/data.fsample-around(1) onsets_session(end)/data.fsample+around(2)];
    data        = ft_selectdata(cfg,data);

    sampleinfo_tmp = data.sampleinfo(1); 

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

    if isubject == 8
    cfg             = [];
    cfg.prewindow   = 3;
    cfg.postwindow  = 3;
    data            = ft_interpolatenan(cfg,data);
    end

    X = data.trial{1}';      % Nsamples x Nchannels
    fs = data.fsample;
    exclusion_range = round(0.01 * fs);

   [xnew, spkindx] = despike_dilate_peak_spike_old(X,2, exclusion_range, false, 0);

%    [xnew, spkindx, spkinfo] = despike_dilate_peak_spike(X, 2, exclusion_range, false, 0, [], fs);

%    figure; plot(spkinfo.sta_time, spkinfo.sta_raw); title('STA raw');
%     figure; plot(spkinfo.sta_time, spkinfo.sta_bc);  title('STA baseline-corrected');

   fs                 = data.fsample;
   [Nsamp, Nchan]     = size(X);

   [spk_samp, spk_chan] = ind2sub(size(X), spkindx);

  
counts = accumarray(spk_chan, 1, [size(X,2) 1]);
disp(counts')

   %% plotting

   fs             = data.fsample;
   [Nsamp, Nchan] = size(X);
   [spk_samp, spk_chan] = ind2sub(size(X), spkindx);


    refrac_ms = 100;                 % e.g., 50 ms
refrac_samp = round(refrac_ms/1000 * fs);

keep = true(size(spk_samp));

for c = 1:Nchan
    idx = find(spk_chan == c);
    [~,ord] = sort(spk_samp(idx));
    idx = idx(ord);

    last = -inf;
    for j = 1:numel(idx)
        s = spk_samp(idx(j));
        if s - last < refrac_samp
            keep(idx(j)) = false;
        else
            last = s;
        end
    end
end

spk_samp = spk_samp(keep);
spk_chan = spk_chan(keep);
spkindx  = spkindx(keep);

counts = accumarray(spk_chan, 1, [size(X,2) 1]);
disp(counts')


   win_sec   = [-0.1 0.1];
   win_samp  = round(win_sec * fs);
   lags      = win_samp(1):win_samp(2);
   nLags     = numel(lags);
   t         = lags / fs;

   % baseline window (pre-spike)
   baseline_sec  = [-0.1 -0.02];
   baseline_samp = round(baseline_sec * fs);
   base_idx      = baseline_samp(1):baseline_samp(2);
   base_mask     = ismember(lags, base_idx);

   STA_chan = zeros(nLags, Nchan);
   count    = zeros(Nchan, 1);
    c = [];
    s = [];
   for k = 1:numel(spk_samp)
       s = spk_samp(k);
       c = spk_chan(k);

       if s + lags(1) < 1 || s + lags(end) > Nsamp
           continue
       end

       seg = X(s + lags, c);

       [b,a] = butter(2, 5/(fs/2), 'high');   % 5 Hz HP
        seg = filtfilt(b,a, seg);

       % epoch-wise baseline correction
       b   = mean(seg(base_mask), 'omitnan');
       seg = seg - b;

%        seg = detrend(seg, 'linear');

       STA_chan(:, c) = STA_chan(:, c) + seg;
       count(c)       = count(c) + 1;
   end

   for c_idx = 1:Nchan
       if count(c_idx) > 0
           STA_chan(:, c_idx) = abs(STA_chan(:, c_idx) ./ count(c_idx));
       else
           STA_chan(:, c_idx) = NaN;
       end
   end

   valid_ch = count > 0;
   STA_mean = nanmean(STA_chan(:, valid_ch), 2);

    perf{isubject}.STA              = STA_chan(:, valid_ch);
    perf{isubject}.STA_mean         = STA_mean;
    perf{isubject}.STA_mean_time    = t;

   figure;
   plot(t, STA_mean, 'k', 'LineWidth', 1.5);
   xlabel('Time from spike (s)');
   ylabel('Voltage (\muV)');
   title('Spike-triggered average (baseline-corrected)');
   xline(0, '--'); grid on;


   figure;
   for ichan = 1:Nchan
       subplot(3,3,ichan)

       plot(t, abs(STA_chan(:,ichan)), 'k', 'LineWidth', 1.5);
       xlabel('Time from spike (s)');
       ylabel('Voltage (\muV)');
       title('Spike-triggered average (all channels)');
       xline(0, '--');
       grid on;
   end
  

   figure;
    plot(t, abs(STA_chan)./max(abs(STA_chan)), 'k', 'LineWidth', 1.5);
    hold on
    plot(t, abs(STA_mean)./max(abs(STA_mean)), 'k', 'LineWidth', 3);

end