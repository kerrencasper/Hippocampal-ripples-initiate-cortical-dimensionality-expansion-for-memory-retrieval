%% plot supplemental material

%                  Casper Kerren      [kerren@cbs.mpg.de]

%%


clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults
[~,ftpath]=ft_version;

paths = config_paths();

settings = [];
settings.base_path_castle   = paths.base_path;
settings.data_dir           = paths.data_dir;
settings.data_dir_channels  = paths.channels_dir;
settings.anatomy_dir        = paths.anatomy_dir;
settings.AAL_dir            = paths.AAL_dir;
settings.SPM_dir            = paths.SPM_dir;

settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

addpath(genpath([paths.MVPA_Light_master]))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([paths.subfunctions])
addpath(genpath([paths.help_functions]))
addpath(genpath([paths.plotting]))

settings.data_dir_channels  = [settings.base_path_castle,'ripple_project_publication_for_replication/templates'];
settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};

settings.colour_scheme_1 = brewermap(30,'RdBu');
settings.colour_scheme_1 = settings.colour_scheme_1;

settings.nu_perm = 4096;


%%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------- SUPPLEMENTAL FIGURE 1 ------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------
%  - plot all channels per participant


addpath('/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/ripple_project_publication_for_replication/additional_analyses/visualisation')

load all_channels

load mean_ripple_per_participant
load std_ripple_per_participant
load nu_ripples_per_participant

load perf_time_resolved_power_spectra

perf = perf_time_resolved_power_spectra;

freqaxis = perf{1}.freq;

% restrict to gamma search range, e.g. 40–140 Hz
gamma_idx = freqaxis >= 55 & freqaxis <= 140;

FOI = nearest(freqaxis,55):nearest(freqaxis,140);
f_gamma = freqaxis(FOI);

powspctrm = [];
for isubject = 1:size(settings.subjects,1)
   
    powspctrm(isubject,:) = perf{isubject}.powspctrm(FOI)/ max(perf{isubject}.powspctrm(FOI));
end



load TF_ripple_aligned

TF = TF_ripple_aligned;

timeaxis    = TF{1,1}.time;
freqaxis     = TF{1,1}.freq;

toi_idx = -1:1/1000:1;

channelsize = 20;
transp      = .05;
col         = 'r';

tmp                             = load([settings.data_dir_channels,'/channels_to_exclude_all_hipp_both_hem.mat']);
channels_to_exclude_all_hipp    = tmp.channels_to_exclude_all_hipp;
tmp                             = load([settings.data_dir_channels,'/channels_hipp_ripples.mat']);
channels_hipp_ripples           = tmp.channels_hipp_ripples;


channelsize = 15;
transp      = .05;
col         = 'r';

this_view   = 'horizontal'; % 'horizontal' 'sagittal' 'coronal'

this_brain  = 'surface_inflated_both';

ld          = load(fullfile(ftpath,'template/anatomy/', [this_brain '.mat']));
brain       = ld.mesh;


epos_all        = [];
subject_colors  = [];  % To store the color for each contact
count_hipp_chan = [];

% Generate a colormap with distinct colors for each subject
num_subjects = size(settings.subjects,1);
cmap = lines(num_subjects);  % Use 'lines' colormap for distinct colors

for isubject = 1:num_subjects
    cfg.channel     = [];
    chan_to_keep    = [];
    idx             = [];
    coordinates     = [];

    cfg.channel                 = {all_channels{isubject}.names};
    chan_to_keep                =  cfg.channel;
    idx                         = ismember(cfg.channel,chan_to_keep);
    count_hipp_chan(isubject,:) = sum(idx);
    coordinates                 = {all_channels{isubject}.coords};
    coordinates                 = coordinates(idx);
    
    epos_all                    = [epos_all; cell2mat(coordinates')];
    
    % Assign the subject's color to the corresponding contacts
    subject_colors = [subject_colors; repmat(cmap(isubject, :), sum(idx), 1)];
end


figure('units','normalized','outerposition',[0 0 1 1]);
views = [180 -90];


for i = 1:12

    % Initialize variables for each iteration
    subject_colors  = [];
    count_hipp_chan = [];
    epos_all        = [];
    cfg.channel     = [];
    chan_to_keep    = [];
    idx             = [];
    coordinates     = [];

    % Select relevant channels
    cfg.channel                 = {all_channels{i}.names};
    chan_to_keep                = intersect(cellstr(setdiff([cfg.channel], settings.scalp_channels)), char(channels_hipp_ripples{i, :}));
%     chan_to_keep                = cfg.channel;
    idx                         = ismember(cfg.channel, chan_to_keep);
    count_hipp_chan(i, :)       = sum(idx);
    coordinates                 = {all_channels{i}.coords};
    coordinates                 = coordinates(idx);
    
    epos_all                    = [epos_all; cell2mat(coordinates')];
    subject_colors = [subject_colors; repmat(cmap(i, :), sum(idx), 1)];

    % Plotting the brain mesh at the top
    if i < 5
        h1 = subplot(6, 4, i);
        pos_h    = get(h1, 'Position');
        pos_h(1) = pos_h(1)-.045;
        pos_h(2) = pos_h(2)-.01;
        pos_h(3) = pos_h(3)+.09;
        pos_h(4) = pos_h(4)+.09;
        set(h1, 'Position', pos_h);


    elseif i >= 5 && i < 9
        h1 = subplot(6, 4, i + 4);
        pos_h    = get(h1, 'Position');
        pos_h(1) = pos_h(1)-.045;
        pos_h(2) = pos_h(2)-.02;
        pos_h(3) = pos_h(3)+.09;
        pos_h(4) = pos_h(4)+.09;
        set(h1, 'Position', pos_h);

    else
        h1 = subplot(6, 4, i + 8);
        pos_h    = get(h1, 'Position');
        pos_h(1) = pos_h(1)-.045;
        pos_h(2) = pos_h(2)-.02;
        pos_h(3) = pos_h(3)+.09;
        pos_h(4) = pos_h(4)+.09;
        set(h1, 'Position', pos_h);

    end

    ft_plot_mesh(brain, 'facealpha', transp, 'facecolor', 'cortex')
    
    for ii = 1:size(epos_all, 1)
        this_MNI = epos_all(ii, :);
        ft_plot_mesh(this_MNI, 'vertexcolor', subject_colors(ii, :), 'vertexsize', channelsize);
    end

    camlight headlight;
    set(gcf, 'color', 'w', 'alpha', 0);
    view(views)

    % Determine subplot position for mean_ripple_per_participant
    if i < 5
        h = subplot(6, 4, i + 4);
        pos = get(h, 'Position'); 
        pos(2) = pos(2)+.01;
        set(h, 'Position', pos);
    elseif i >= 5 && i < 9
        h = subplot(6, 4, i + 8);
        pos = get(h, 'Position'); 
         pos(2) = pos(2)-.03;
        set(h, 'Position', pos);

        pos_h(2) = pos_h(2)-.04;
        set(h1, 'Position', pos_h);

    else
        h = subplot(6, 4, i + 12);
        pos = get(h, 'Position'); 
        pos(2) = pos(2)-.06;
        set(h, 'Position', pos);

        pos_h(2) = pos_h(2)-.08;
        set(h1, 'Position', pos_h);

    end

    % Plot mean_ripple_per_participant with bounded line
    boundedline(toi_idx, mean_ripple_per_participant(i, :), std_ripple_per_participant(i, :), 'k', 'alpha')
    xlabel('Time (sec)')
    ylabel('\muV')
    title(sprintf('patient %d, %d events', i, nu_ripples_per_participant(i, :)), 'interpreter', 'none')
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', 16)
    axis tight


    h2 = subplot('Position',[.01, .01, .01, .01]);
    pos2    = get(h2, 'Position');
    pos2(1) = pos(1)+.12;
    pos2(2) = pos(2)+.06;
    pos2(3) = .04;
    pos2(4) = .04;
    set(h2, 'Position', pos2);
    to_plot_ripple = squeeze(mean(TF{i}.all));


    to_plot_ripple = [];
    to_plot_ripple = squeeze(mean(TF{i}.all));

    imagesc(...
        timeaxis,...
        freqaxis,...
        to_plot_ripple);
    colormap(hot)
    caxis([0 .7])
    hline(90,'w')
    axis off


    axis xy
    set(gca,'TickDir','out')

    h3 = subplot('Position',[.01, .01, .01, .01]);
    pos3    = get(h3, 'Position');
    pos3(1) = pos(1)+.01;
    pos3(2) = pos(2)+.06;
    pos3(3) = .04;
    pos3(4) = .04;
    set(h3, 'Position', pos3);

    plot(f_gamma,powspctrm(i,:),'k')
   
end

clearvars -except settings

%% dimensionality (100ms)

% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------- SUPPLEMENTAL FIGURE 3 ------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------


figure('units','normalized','outerposition',[0 0 1 1]);

% dimensionality 100ms sliding window
load perf_dimensionality_100ms

subplot(3,3,1)

perf = perf_dimensionality_100ms;


correct_incorrect = {};

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(isubject,1,:)   = perf{isubject}.correct.accuracy;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = perf{isubject}.incorrect.accuracy;
  
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);


xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

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

cfg.numrandomization    = 1000;%'all';

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

d = squeeze(correct_incorrect{1}.individual-correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(perf{1,1}.time_train,m,s,'k','cmap', flipud(settings.colour_scheme_1(23,:)));
plot(perf{1,1}.time_train,m,'k','linewidth',3);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',20)
xlabel('ripple time (sec)')
ylabel('dim. difference')
set(gca,'TickDir','out')
xticks([-.9 0 .9])
xticklabels({'-1', '0', '1'})

axis tight
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])

%% dimensionality 200ms sliding window
load perf_dimensionality_200ms

subplot(3,3,2)

perf = perf_dimensionality_200ms;



correct_incorrect = {};

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(isubject,1,:)   = perf{isubject}.correct.accuracy;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = perf{isubject}.incorrect.accuracy;
  
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);


xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

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

cfg.numrandomization    = 1000;%'all';

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

d = squeeze(correct_incorrect{1}.individual-correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(perf{1,1}.time_train,m,s,'k','cmap', flipud(settings.colour_scheme_1(23,:)));
plot(perf{1,1}.time_train,m,'k','linewidth',3);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',20)
xlabel('ripple time (sec)')
ylabel('dim. difference')
set(gca,'TickDir','out')
xticks([-.9 0 .9])
xticklabels({'-1', '0', '1'})

axis tight
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])




clearvars -except settings




%% reconstructing using components



load perf_PCA_LDA

perf = perf_PCA_LDA;

freqaxis = perf{1, 1}.time_train;
timeaxis = perf{1, 1}.time_test;

    correct                 = cell(1,size(settings.subjects,1));
    incorrect               = cell(1,size(settings.subjects,1));
    correct_dec     = [];
    incorrect_dec   = [];
sigma = .5; % STD of 2D gaussian smoothing
    for isubject = 1:size(settings.subjects,1)

        % correct
        correct{1,isubject}                                 = struct;
        correct{1,isubject}.label                           = {'chan'};
        correct{1,isubject}.dimord                          = 'chan_freq_time';
        correct{1,isubject}.freq                            = perf{1, 1}.time_train;
        correct{1,isubject}.time                            = perf{1, 1}.time_test;
        correct{1,isubject}.powspctrm(1,:,:)                = perf{isubject}.correct.accuracy;


        % incorrect
        incorrect{1,isubject}                               = correct{1,isubject};
        incorrect{1,isubject}.powspctrm(1,:,:)              = perf{isubject}.incorrect.accuracy;

        baseline_all{1,isubject}                            = correct{1,isubject};
        baseline_all{1,isubject}.powspctrm(1,:,:)           = .5*ones(size(perf{isubject}.correct.accuracy));

        correct_dec(isubject,:,:)                           = imgaussfilt(perf{isubject}.correct.accuracy,sigma);
        incorrect_dec(isubject,:,:)                         = imgaussfilt(perf{isubject}.incorrect.accuracy,sigma);
    end


%% FT stats (decoding)

cfg                     = [];
cfg.latency             = [-1 1];
cfg.frequency           = [-.2 3]; % encoding time
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

cfg.numrandomization    = 1000;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = size(settings.subjects,1);
% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

[Fieldtripstats] = ft_freqstatistics(cfg, correct{:}, incorrect{:});

length(find(Fieldtripstats.mask==1))


stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);


%% Correlate with original data (still debating whether we want to include this or not)

% get size of the different decoding analyses

xlimits = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);

% limit decoding from reconstructed data
correct_dec_rec     = correct_dec(:, :, xlimits);
incorrect_dec_rec   = incorrect_dec(:, :, xlimits);
% load original data
load('tvals_coarse_enc_ripple.mat')
load('mask_t_vals_coarse_enc_ripple.mat')

tvals_single_ripple_coarse  = tvals_coarse_enc_ripple;
mask_t_vals_coarse          = squeeze(mask_t_vals_coarse_enc_ripple);

% Load data from decoding analysis
load('correct_dec_all.mat')
load('incorrect_dec_all.mat')

correct_dec_all = correct_dec_all;
incorrect_dec_all = incorrect_dec_all;

[rows_large, cols_large] = size(tvals_single_ripple_coarse);

% Get the size of the smaller matrix
[rows_small, cols_small] = size(tvals);

% Calculate the number of columns to add
cols_to_add = cols_large - cols_small;

% Pad the smaller matrix with zeros to match the size of the larger one
% (only a few data points)
tvals_reconstructed_resized                 = [tvals, zeros(rows_small, cols_to_add)];
correct_dec_rec(:,:,end:end+cols_to_add)    = 0;
incorrect_dec_rec(:,:,end:end+cols_to_add)  = 0;
TOI_to_plot                                 = linspace(-1,1, size(tvals_single_ripple_coarse,2));
FOI_to_plot                                 = correct{1, 1}.freq;


[r, p] = corr(tvals_single_ripple_coarse(logical(mask_t_vals_coarse)),tvals_reconstructed_resized(logical(mask_t_vals_coarse)))



for isubject = 1:size(settings.subjects,1)
    tmp_sub = [];
    tmp_sub = squeeze(correct_dec_rec(isubject, :,:));
%     tmp_sub(~mask_t_vals_coarse) = NaN;
    correct_dec_rec(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(correct_dec_all(isubject, :,:));
%     tmp_sub(~mask_t_vals_coarse) = NaN;
    correct_dec_all(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(incorrect_dec_rec(isubject, :,:));
%     tmp_sub(~mask_t_vals_coarse) = NaN;
    incorrect_dec_rec(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(incorrect_dec_all(isubject, :,:));
%     tmp_sub(~mask_t_vals_coarse) = NaN;
    incorrect_dec_all(isubject, :,:) = tmp_sub;
end

correct_dec_mean    = squeeze(mean(correct_dec_all,2));
incorrect_dec_mean  = squeeze(mean(incorrect_dec_all,2));

correct_dec_rec_mean    = squeeze(mean(correct_dec_rec,2));
incorrect_dec_rec_mean  = squeeze(mean(incorrect_dec_rec,2));

[rho_ppp, p_ppp] = corr(correct_dec_mean-incorrect_dec_mean,correct_dec_rec_mean-incorrect_dec_rec_mean,'type', 'Spearman');
subplot(3,3,3)
 imagesc(...
    timeaxis(stats_time),...
    timeaxis(stats_time),...
    rho_ppp(1:end-5,1:end-5))
colormap(flipud(settings.colour_scheme_1))
caxis([0 1])
vline(0)


axis xy
set(gca,'FontSize',18)
ylabel('ripple time (sec) original')
xlabel('ripple time (sec) reconstructed')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');

ylabel(hcb, 'Correlation', 'FontSize', 20); 

correct_dec_mean    = mean(mean(correct_dec_all,3),2);
incorrect_dec_mean  = mean(mean(incorrect_dec_all,3),2);

correct_dec_rec_mean    = mean(mean(correct_dec_rec,3),2);
incorrect_dec_rec_mean  = mean(mean(incorrect_dec_rec,3),2);

[rho_ppp, p_ppp] = corr(correct_dec_mean-incorrect_dec_mean,correct_dec_rec_mean-incorrect_dec_rec_mean,'type', 'Spearman');




%% power law

hold on
subplot(3,3,4)


load perf_dimensionality_fit_power_law

perf = perf_dim_power_law;


correct_incorrect = {};

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(isubject,1,:)   = perf{isubject}.correct.accuracy;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = perf{isubject}.incorrect.accuracy;
  
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);

%%
xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

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

xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);

d = squeeze(correct_incorrect{1}.individual(:,xlimits))-squeeze(correct_incorrect{2}.individual(:,xlimits));

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(correct_incorrect{1, 1}.time(xlimits)  ,m,s, 'cmap', settings.colour_scheme_1(23,:));
plot(correct_incorrect{1, 1}.time(xlimits)  ,m,'k','linewidth',3);
hold on


stats_time = nearest(correct_incorrect{1}.time(xlimits),cfg.latency(1)):nearest(correct_incorrect{1}.time(xlimits),cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time(xlimits)));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time(xlimits),sigline,'k','linewidth',4);

set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Distance from power law')
set(gca,'TickDir','out')
xticks([-.95 0 .95])
xticklabels({'-1', '0', '1'})
vline(0)
hline(0)



%% dimensionality LME effective dimensionality


load perf_dimensionality_effective_dim

perf = perf_dimensionality_effective_dim;

load RT_all_subj_correct
load RT_all_subj_incorrect

RT_correct = cellfun(@mean, RT_all_subj_correct);
RT_incorrect = cellfun(@mean, RT_all_subj_incorrect);

TOI = nearest(perf{1, 1}.time_test,-1):nearest(perf{1, 1}.time_test,1);

subjectIDs  = [];
blockIDs    = [];
conditions  = [];
halves      = [];
accuracies  = [];

for isubject = 1:numel(perf)
    perf_correct     = perf{isubject}.correct.accuracy(:,TOI);
    perf_incorrect   = perf{isubject}.incorrect.accuracy(:,TOI);

    idx = (perf_correct(:,1) ~= 0 & perf_incorrect(:,1) ~= 0);

    perf_correct    = perf_correct(idx,:);
    perf_incorrect  = perf_incorrect(idx,:);

    % Correct trials
    num_blocks_correct = size(perf_correct, 1);
    for iblock = 1:num_blocks_correct
        num_samples = size(perf_correct, 2);
        half_point  = floor(num_samples / 2);
        % First half
        subjectIDs  = [subjectIDs; repmat(isubject, half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, half_point, 1)];
        conditions  = [conditions; repmat(1, half_point, 1)];  % 1 for correct
        halves      = [halves; repmat(1, half_point, 1)];  % 1 for first half
        accuracies  = [accuracies; perf_correct(iblock, 1:half_point)'];
        % Second half
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples - half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples - half_point, 1)];
        conditions  = [conditions; repmat(1, num_samples - half_point, 1)];  % 1 for correct
        halves      = [halves; repmat(2, num_samples - half_point, 1)];  % 2 for second half
        accuracies  = [accuracies; perf_correct(iblock, half_point+1:end)'];
    end

    % Incorrect trials
    num_blocks_incorrect = size(perf_incorrect, 1);
    for iblock = 1:num_blocks_incorrect
        num_samples = size(perf_incorrect, 2);
        half_point  = floor(num_samples / 2);
        % First half
        subjectIDs  = [subjectIDs; repmat(isubject, half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, half_point, 1)];
        conditions  = [conditions; repmat(2, half_point, 1)];  % 2 for incorrect
        halves      = [halves; repmat(1, half_point, 1)];  % 1 for first half
        accuracies  = [accuracies; perf_incorrect(iblock, 1:half_point)'];
        % Second half
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples - half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples - half_point, 1)];
        conditions  = [conditions; repmat(2, num_samples - half_point, 1)];  % 2 for incorrect
        halves      = [halves; repmat(2, num_samples - half_point, 1)];  % 2 for second half
        accuracies  = [accuracies; perf_incorrect(iblock, half_point+1:end)'];
    end
end

% Create a table
dataTable = table(subjectIDs, blockIDs, conditions, halves, accuracies, ...
    'VariableNames', {'SubjectID', 'BlockID', 'Condition', 'Half', 'Accuracy'});


% Include interaction between Condition and Half
lme = fitlme(dataTable, 'Accuracy ~ Condition * Half + (1|SubjectID) + (1|BlockID)');

% Display the results
disp(lme)

%% PLOT
% Extract unique subject IDs
uniqueSubjects = unique(dataTable.SubjectID);

% Initialize arrays to store the means and SEMs
means = zeros(numel(uniqueSubjects), 2, 2);  % Dimensions: subject x condition x half
sems = zeros(numel(uniqueSubjects), 2, 2);

% Loop through each subject, condition, and half to calculate means and SEMs
for isubject = 1:numel(uniqueSubjects)
    for condition = 1:2
        for half = 1:2
            % Filter data for the current subject, condition, and half
            subset = dataTable(dataTable.SubjectID == uniqueSubjects(isubject) & ...
                dataTable.Condition == condition & ...
                dataTable.Half == half, :);
            % Calculate the mean and SEM
            means(isubject, condition, half) = mean(subset.Accuracy);
            sems(isubject, condition, half) = std(subset.Accuracy) / sqrt(height(subset));
        end
    end
end


% Create labels for the plot
conditionLabels = {'Correct', 'Incorrect'};
halfLabels = {'First Half', 'Second Half'};

% Plot the interaction effect for each participant

mean_correct    = squeeze(means(:,1,:));
mean_incorrect  = squeeze(means(:,2,:));

data_lme = {};
data_lme{1, 1} = mean_correct(:,2)-mean_correct(:,1);
data_lme{2, 1} = mean_incorrect(:,2)-mean_incorrect(:,1);

subplot(3,3,5)

N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector

scattercorrect = scatter(.9, mean_correct(:,2)-mean_correct(:,1), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;


scattercorrect = scatter(1.1, mean_incorrect(:,2)-mean_incorrect(:,1), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;





hold on;
[f, xi] = ksdensity(mean_correct(:,2)-mean_correct(:,1),'Bandwidth',.1); 
fill(.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2)

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(mean_incorrect(:,2)-mean_incorrect(:,1),'Bandwidth',.1); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;

    line([.9 1.1], [mean_correct(:,2)-mean_correct(:,1), mean_incorrect(:,2)-mean_incorrect(:,1)], 'Color', [0.5 0.5 0.5]);

    
end


hold on
set(gca, 'TickDir', 'out');
ylabel(sprintf('Dimensionality change\n(post-pre ripple)')); 
xticks([.8 1.2])
xticklabels({'AM+', 'AM-'}); 
xlabel('Trial type');
set(gca, 'FontSize', 20);


hline(0)
box off



[a, b, c, d] = ttest(data_lme{1, 1},0)


[a, b, c, d] = ttest(data_lme{2, 1},0)


[a, b, c, d] = ttest(data_lme{1, 1},data_lme{2, 1})



%% LME MP dim

load perf_dim_LME_MPdim

perf = perf_dim_LME_MPdim;

load RT_all_subj_correct
load RT_all_subj_incorrect

RT_correct = cellfun(@mean, RT_all_subj_correct);
RT_incorrect = cellfun(@mean, RT_all_subj_incorrect);

TOI = nearest(perf{1, 1}.time_test,-1):nearest(perf{1, 1}.time_test,1);

subjectIDs  = [];
blockIDs    = [];
conditions  = [];
halves      = [];
accuracies  = [];

for isubject = 1:numel(perf)
    perf_correct     = perf{isubject}.correct.accuracy(:,TOI);
    perf_incorrect   = perf{isubject}.incorrect.accuracy(:,TOI);

    idx = (perf_correct(:,1) ~= 0 & perf_incorrect(:,1) ~= 0);

    perf_correct    = perf_correct(idx,:);
    perf_incorrect  = perf_incorrect(idx,:);

    % Correct trials
    num_blocks_correct = size(perf_correct, 1);
    for iblock = 1:num_blocks_correct
        num_samples = size(perf_correct, 2);
        half_point  = floor(num_samples / 2);
        % First half
        subjectIDs  = [subjectIDs; repmat(isubject, half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, half_point, 1)];
        conditions  = [conditions; repmat(1, half_point, 1)];  % 1 for correct
        halves      = [halves; repmat(1, half_point, 1)];  % 1 for first half
        accuracies  = [accuracies; perf_correct(iblock, 1:half_point)'];
        % Second half
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples - half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples - half_point, 1)];
        conditions  = [conditions; repmat(1, num_samples - half_point, 1)];  % 1 for correct
        halves      = [halves; repmat(2, num_samples - half_point, 1)];  % 2 for second half
        accuracies  = [accuracies; perf_correct(iblock, half_point+1:end)'];
    end

    % Incorrect trials
    num_blocks_incorrect = size(perf_incorrect, 1);
    for iblock = 1:num_blocks_incorrect
        num_samples = size(perf_incorrect, 2);
        half_point  = floor(num_samples / 2);
        % First half
        subjectIDs  = [subjectIDs; repmat(isubject, half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, half_point, 1)];
        conditions  = [conditions; repmat(2, half_point, 1)];  % 2 for incorrect
        halves      = [halves; repmat(1, half_point, 1)];  % 1 for first half
        accuracies  = [accuracies; perf_incorrect(iblock, 1:half_point)'];
        % Second half
        subjectIDs  = [subjectIDs; repmat(isubject, num_samples - half_point, 1)];
        blockIDs    = [blockIDs; repmat(iblock, num_samples - half_point, 1)];
        conditions  = [conditions; repmat(2, num_samples - half_point, 1)];  % 2 for incorrect
        halves      = [halves; repmat(2, num_samples - half_point, 1)];  % 2 for second half
        accuracies  = [accuracies; perf_incorrect(iblock, half_point+1:end)'];
    end
end

% Create a table
dataTable = table(subjectIDs, blockIDs, conditions, halves, accuracies, ...
    'VariableNames', {'SubjectID', 'BlockID', 'Condition', 'Half', 'Accuracy'});


% Include interaction between Condition and Half
lme = fitlme(dataTable, 'Accuracy ~ Condition * Half + (1|SubjectID) + (1|BlockID)');

% Display the results
disp(lme)

%% PLOT
% Extract unique subject IDs
uniqueSubjects = unique(dataTable.SubjectID);

% Initialize arrays to store the means and SEMs
means = zeros(numel(uniqueSubjects), 2, 2);  % Dimensions: subject x condition x half
sems = zeros(numel(uniqueSubjects), 2, 2);

% Loop through each subject, condition, and half to calculate means and SEMs
for isubject = 1:numel(uniqueSubjects)
    for condition = 1:2
        for half = 1:2
            % Filter data for the current subject, condition, and half
            subset = dataTable(dataTable.SubjectID == uniqueSubjects(isubject) & ...
                dataTable.Condition == condition & ...
                dataTable.Half == half, :);
            % Calculate the mean and SEM
            means(isubject, condition, half) = mean(subset.Accuracy);
            sems(isubject, condition, half) = std(subset.Accuracy) / sqrt(height(subset));
        end
    end
end


% Create labels for the plot
conditionLabels = {'Correct', 'Incorrect'};
halfLabels = {'First Half', 'Second Half'};

% Plot the interaction effect for each participant

mean_correct    = squeeze(means(:,1,:));
mean_incorrect  = squeeze(means(:,2,:));

data_lme = {};
data_lme{1, 1} = mean_correct(:,2)-mean_correct(:,1);
data_lme{2, 1} = mean_incorrect(:,2)-mean_incorrect(:,1);

subplot(3,3,6)

N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector

scattercorrect = scatter(.9, mean_correct(:,2)-mean_correct(:,1), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;


scattercorrect = scatter(1.1, mean_incorrect(:,2)-mean_incorrect(:,1), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;





hold on;
[f, xi] = ksdensity(mean_correct(:,2)-mean_correct(:,1),'Bandwidth',.1); 
fill(.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2)

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(mean_incorrect(:,2)-mean_incorrect(:,1),'Bandwidth',.1); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;

    line([.9 1.1], [mean_correct(:,2)-mean_correct(:,1), mean_incorrect(:,2)-mean_incorrect(:,1)], 'Color', [0.5 0.5 0.5]);

    
end


hold on
set(gca, 'TickDir', 'out');
ylabel(sprintf('Dimensionality change\n(post-pre ripple)')); 
xticks([.8 1.2])
xticklabels({'AM+', 'AM-'}); 
xlabel('Trial type');
set(gca, 'FontSize', 20);


hline(0)
box off



[a, b, c, d] = ttest(data_lme{1, 1},0)


[a, b, c, d] = ttest(data_lme{2, 1},0)


[a, b, c, d] = ttest(data_lme{1, 1},data_lme{2, 1})




%% Supplemental spike-triggered analyses

figure('units','normalized','outerposition',[0 0 1 1]);


load perf_dim_spike_trig_avg

perf = perf_dim_spike_trig_avg;

correct_incorrect = {};
num_trial_correct = [];
num_trial_incorrect = [];
to_plot = [];
count_trls = 1;
for isubject = 1:numel(perf)

    if ~isfield(perf{isubject}, 'correct') | ~isfield(perf{isubject}, 'incorrect') |  ~isfield(perf{isubject}.correct,'accuracy') | ~isfield(perf{isubject}.incorrect,'accuracy')
        continue;
    end

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(count_trls,1,:)   = perf{isubject}.correct.accuracy;
%     correct_incorrect{1}.individual(isubject,1,:)   = explained_var_corr(isubject,:);
    correct_incorrect{1}.dimord              = 'subj_chan_time';
    num_trial_correct(count_trls,:) = perf{isubject}.correct.trl_num_test;
    num_trial_incorrect(count_trls,:) = perf{isubject}.incorrect.trl_num_test;

    to_plot(count_trls,:) = abs(perf{isubject}.STA_mean)./max(abs(perf{isubject}.STA_mean));

    count_trls = count_trls+1;

end

to_plot_time = perf{1}.STA_mean_time;

correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);


correct_incorrect{2}                = correct_incorrect{1};

count_trls = 1;
for isubject = 1:numel(perf)

     
    if ~isfield(perf{isubject}, 'correct') | ~isfield(perf{isubject}, 'incorrect') |  ~isfield(perf{isubject}.correct,'accuracy') | ~isfield(perf{isubject}.incorrect,'accuracy')
        continue;
    end

    correct_incorrect{2}.individual(count_trls,1,:)           = perf{isubject}.incorrect.accuracy;
%     correct_incorrect{2}.individual(isubject,1,:)           = explained_var_incorr(isubject,:);
 count_trls = count_trls+1;
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);

%% plot spike-triggered average


nu_spikes_per_participant = round(mean(num_trial_correct,2)+mean(num_trial_incorrect,2));

for isubject = 1:numel(perf)

    
    subplot(4,3,isubject)

    tmp = (abs(perf{isubject}.STA));
    
    m = abs(perf{isubject}.STA_mean)./max(abs(perf{isubject}.STA_mean));
    s = std((abs(perf{isubject}.STA)')./max(abs(perf{isubject}.STA)'))./sqrt(size(perf{isubject}.STA,2));

%     plot(to_plot_time,abs(perf{isubject}.STA_mean),'MarkerFaceColor',[.5 .5 .5])    

%     hold on
    shadedErrorBar(to_plot_time,m,s);
%     plot(to_plot_time,to_plot(isubject,:),'k', 'LineWidth',3)
    title(sprintf('patient %d, %d analysed spikes', isubject, nu_spikes_per_participant(isubject, :)), 'interpreter', 'none')
    set(gca, 'FontSize', 20);
end


[a,b,c,d] = ttest(mean(num_trial_correct,2),mean(num_trial_incorrect,2));

%% STATS

figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,3,1)
xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

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

cfg.numrandomization    = 'all';%'all';

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


boundedline(perf{1,1}.time_test,m,s,'b');
plot(perf{1,1}.time_test,m,'k','linewidth',2);
hold on

d = squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));


boundedline(perf{1,1}.time_train,m,s,'r');
plot(perf{1,1}.time_train,m,'k','linewidth',2);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 2;

plot(correct_incorrect{1}.time,sigline,'r','linewidth',4);

set(gca,'FontSize',16,'FontName','Arial')
xlabel('ripple time (s)')
ylabel('dimensionality difference')
set(gca,'TickDir','out')

axis tight
vline(0)


xlim([cfg.latency(1), cfg.latency(end)])


%% decoding

load perf_dec_spike_trig_avg

perf = perf_dec_spike_trig_avg;





correct_dec_all     = [];
incorrect_dec_all   = [];

sigma = .5; % STD of 2D gaussian smoothing

correct                 = cell(1);
incorrect              = cell(1);
correct_colour         = cell(1);
incorrect_colour        = cell(1);
correct_scene            = cell(1);
incorrect_scene         = cell(1);
baseline_all            = cell(1);
correct_scene_colour     = cell(1);
incorrect_scene_colour   = cell(1);


sigma = .5; % STD of 2D gaussian smoothing

count_trls = 1;
for isubject = 1:numel(subjects)

    if ~isfield(perf{isubject}, 'correct') | ~isfield(perf{isubject}, 'incorrect') |...
            ~isfield(perf{isubject}.correct,'accuracy') | ~isfield(perf{isubject}.incorrect,'accuracy') |...
            ~isfield(perf{isubject}.colour.correct,'accuracy') | ~isfield(perf{isubject}.colour.incorrect,'accuracy') |...
            ~isfield(perf{isubject}.scene.correct,'accuracy') | ~isfield(perf{isubject}.scene.incorrect,'accuracy')
        continue;
    end

    nu_trials(count_trls,:) = sum(perf{isubject}.correct.trl_num_test);

    % correct
    correct{1,count_trls}                                 = struct;
    correct{1,count_trls}.label                           = {'chan'};
    correct{1,count_trls}.dimord                          = 'chan_freq_time';
    correct{1,count_trls}.freq                            = perf{1,1}.time_train;
    correct{1,count_trls}.time                            = perf{1,1}.time_test;
    correct{1,count_trls}.powspctrm(1,:,:)                = perf{isubject}.correct.accuracy;
   

    correct_dec_all(count_trls,:,:) = imgaussfilt(perf{isubject}.correct.accuracy, sigma);

    % incorrect
    incorrect{1,count_trls}                               = correct{1,count_trls};
    incorrect{1,count_trls}.powspctrm(1,:,:)          = perf{isubject}.incorrect.accuracy;
   
    incorrect_dec_all(count_trls,:,:) = imgaussfilt(perf{isubject}.incorrect.accuracy, sigma);

    % correct colour
    correct_colour{1,count_trls}                          = correct{1,count_trls};
    correct_colour{1,count_trls}.powspctrm(1,:,:)     = perf{isubject}.colour.correct.accuracy;
  

    % incorrect colour
    incorrect_colour{1,count_trls}                        = correct{1,count_trls};
    incorrect_colour{1,count_trls}.powspctrm(1,:,:)   = perf{isubject}.colour.incorrect.accuracy;
   
    % correct scene
    correct_scene{1,count_trls}                           = correct{1,count_trls};
    correct_scene{1,count_trls}.powspctrm(1,:,:)     = perf{isubject}.scene.correct.accuracy;
   
    % incorrect scene
    incorrect_scene{1,count_trls}                        = correct{1,count_trls};
    incorrect_scene{1,count_trls}.powspctrm(1,:,:)   = perf{isubject}.scene.incorrect.accuracy;
   
    % correct scene and colour collapsed
    correct_scene_colour{1,count_trls}                    = correct{1,count_trls};
    correct_scene_colour{1,count_trls}.powspctrm(1,:,:) = mean(cat(1,correct_colour{1,count_trls}.powspctrm,correct_scene{1,count_trls}.powspctrm));
   
    correct_dec_fine_all(count_trls,:,:) = imgaussfilt(mean(cat(1,correct_colour{1,count_trls}.powspctrm,correct_scene{1,count_trls}.powspctrm)), sigma);

    % incorrect scene and colour collapsed
    incorrect_scene_colour{1,count_trls}                  = correct{1,count_trls};
    incorrect_scene_colour{1,count_trls}.powspctrm(1,:,:) = mean(cat(1,incorrect_colour{1,count_trls}.powspctrm,incorrect_scene{1,count_trls}.powspctrm));
  
    incorrect_dec_fine_all(count_trls,:,:) = imgaussfilt(mean(cat(1,incorrect_colour{1,count_trls}.powspctrm,incorrect_scene{1,count_trls}.powspctrm)), sigma);


    % baseline
    baseline_all{1,count_trls}                            = correct{1,count_trls};
    baseline_all{1,count_trls}.powspctrm(1,:,:)           = .5*ones(size(perf{isubject}.correct.accuracy));

    count_trls = count_trls+1;
end


%% Run stats scene and colour collapsed (fine-grained)

timeaxis = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);
timeaxis = correct{1, 1}.time(timeaxis);



contrast = 'incorrect'; % baseline incorrect

cfg                     = [];
cfg.latency             = [timeaxis(1) timeaxis(end)]; % ripple time
cfg.frequency           = [-.2 freqaxis(end)]; % encoding time  [freqaxis(1) freqaxis(end)]
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

cfg.numrandomization    = 500;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = numel(correct);
% set up design matrix

% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];
cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
if strcmp(contrast,'baseline')
    [Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, baseline_all{:});
elseif strcmp(contrast,'incorrect')
    [Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, incorrect_scene_colour{:});
end
length(find(Fieldtripstats.mask==1))

%% plot (significant) t vals

stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);

subplot(2,3,2)
imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    tvals);
colormap(jet)
caxis([-5 5])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',16)
xlabel('ripple time (sec)')
ylabel('encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');

if any(Fieldtripstats.mask(:))
    plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats.mask)))
end




%% Gamma symmetry

clear
close all


% path settings

paths = config_paths();

settings = [];
settings.base_path_castle = paths.base_path;
settings.data_dir           = paths.data_dir;
settings.data_dir_channels  = paths.channels_dir;
settings.anatomy_dir        = paths.anatomy_dir;
settings.AAL_dir            = paths.AAL_dir;
settings.SPM_dir            = paths.SPM_dir;

load("colour_scheme.mat")
settings.colour_scheme = colour_scheme;

settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');
subjects                    = {'CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK'};
SubjectIDs                  = {'01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK'};

addpath(genpath([paths.MVPA_Light_master]))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([paths.subfunctions])
addpath(genpath([paths.help_functions]))
addpath(genpath([paths.plotting]))

% ripple settings

settings.remove_falsepositives      = 1; % decide whether or not to exclude ripples deemed false positives based on spectral peak detection
settings.full_enc_trial             = 1; % set to 0 if you want encoding trial to end with RT and to 1 if it should end at 3 sec
settings.remove_ripple_duplicates   = 1; % remove co-occuring ripples
settings.time_to_excl_RT            = .25; % exclude last 250 ms of trials, to make sure ripple event was in trial
settings.rippleselection            = 0; % use all (0), short (1) or long (2) ripples
settings.ripples_most_chan          = 0; % use only ripples from the channel with greatest number of ripples
settings.solo_ripple                = 1; % pick one ripple per trial if multiple ripple events are found
settings.ripple_latency             = [.25 5]; % define time window at retrieval in which the ripple events need to occur, eg [.5 1.5]

% TFR settings

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

figure('units','normalized','outerposition',[0 0 1 1]);

load PAC_gamma_symmetry

PAC = PAC_gamma_symmetry;


phase_to_plot       = [];
pow_to_plot         = [];
phase_to_plot_avg       = [];
pow_to_plot_avg         = [];
% figure;
for participant = 1:numel(subjects)
%     hold on
%     subplot(4,3,participant)

    phase_to_plot{participant}       = PAC{participant}.max_theta;
    pow_to_plot{participant}         = PAC{participant}.max_gamma;

    phase_to_plot_avg(participant,:)       = PAC{participant}.max_theta_average;
    pow_to_plot_avg(participant,:)         = PAC{participant}.max_gamma_average;

%     hist(PAC{participant}.max_theta)

end


theta_dispersion   = [];
gamma_dispersion   = [];

for participant = 1:numel(subjects)

    theta_vals = PAC{participant}.max_theta;
    theta_avg  = PAC{participant}.max_theta_average;

    gamma_vals = PAC{participant}.max_gamma;
    gamma_avg  = PAC{participant}.max_gamma_average;

    % --- Theta dispersion ---
    theta_dispersion.std(participant)  = std(theta_vals);
    theta_dispersion.mad(participant)  = mad(theta_vals,1);   % robust measure
    theta_dispersion.range(participant,:) = [min(theta_vals), max(theta_vals)];
    theta_dispersion.within2Hz(participant) = mean(abs(theta_vals - theta_avg) <= 2);

    % --- Gamma dispersion ---
    gamma_dispersion.std(participant)  = std(gamma_vals);
    gamma_dispersion.mad(participant)  = mad(gamma_vals,1);
    gamma_dispersion.range(participant,:) = [min(gamma_vals), max(gamma_vals)];
    gamma_dispersion.within10Hz(participant) = mean(abs(gamma_vals - gamma_avg) <= 10);

end


mean(theta_dispersion.std)
mean(theta_dispersion.mad)
mean(theta_dispersion.within2Hz)



mean(gamma_dispersion.std)
mean(gamma_dispersion.mad)
mean(gamma_dispersion.within10Hz)




%% waveform metrics

all_peaktrough = [];
all_risedecay = [];
all_skew = [];
all_kurt = [];
all_gammafreq = [];

for s = 1:numel(PAC)
    all_peaktrough = [all_peaktrough; PAC{s}.peak_trough_ratio(:)];
    all_risedecay = [all_risedecay; PAC{s}.rise_decay_ratio(:)];
    all_skew = [all_skew; PAC{s}.skew(:)];
    all_kurt = [all_kurt; PAC{s}.kurt(:)];
    all_gammafreq = [all_gammafreq; PAC{s}.max_gamma(:)];
end

all_peaktrough_ppp = [];
all_risedecay_ppp = [];
all_skew_ppp = [];
all_kurt_ppp = [];
all_gammafreq_ppp = [];
gamma_to_plot = [];

for s = 1:numel(PAC)
    all_peaktrough_ppp(:,s) = mean(PAC{s}.peak_trough_ratio(:));
    all_risedecay_ppp(:,s) = mean(PAC{s}.rise_decay_ratio(:));
    all_skew_ppp(:,s) = mean(PAC{s}.skew(:));
    all_kurt_ppp(:,s) = mean(PAC{s}.kurt(:));
    all_gammafreq_ppp (:,s)= mean(PAC{s}.max_gamma(:));
    gamma_to_plot(s,:)         = PAC{s}.gamma_wave_grand_mean;

end


m = mean(gamma_to_plot);
s = std(gamma_to_plot)/sqrt(12);
subplot(3,2,1)
plot(x, gamma_to_plot', 'Color', [0.4 0.4 0.4 0.3]);
hold on
plot(x, m, 'k', 'LineWidth', 2);
xlabel(x_label);
ylabel('Gamma amplitude');
title('Cortical gamma waveform');
set(gca,'FontSize',20)
set(gca,'TickDir','out')


box off

mean_peaktr = mean(all_peaktrough,'omitnan')
std_peaktr  = std(all_peaktrough,'omitnan')

mean_risedecay = mean(all_risedecay,'omitnan')
std_risedecay  = std(all_risedecay,'omitnan')

mean_skew = mean(all_skew,'omitnan')
std_skew  = std(all_skew,'omitnan')

mean_kurt = mean(all_kurt,'omitnan')
std_kurt  = std(all_kurt,'omitnan')



[r_pt, p_pt] = corr(all_peaktrough, all_gammafreq, 'rows','complete');
[r_rd, p_rd] = corr(all_risedecay, all_gammafreq, 'rows','complete');
[r_sk, p_sk] = corr(all_skew, all_gammafreq, 'rows','complete');
[r_ku, p_ku] = corr(all_kurt, all_gammafreq, 'rows','complete');

fprintf('Peak–trough corr = %.3f, p=%.3g\n', r_pt, p_pt);
fprintf('Rise–decay corr = %.3f, p=%.3g\n', r_rd, p_rd);
fprintf('Skew corr = %.3f, p=%.3g\n', r_sk, p_sk);
fprintf('Kurtosis corr = %.3f, p=%.3g\n', r_ku, p_ku);



subplot(3,2,2); nhist(all_peaktrough,'proportion'); title('Peak–Trough Ratio');
set(gca,'FontSize',20)
set(gca,'TickDir','out')
subplot(3,2,3); nhist(all_risedecay,'proportion'); title('Rise-Decay Ratio');
set(gca,'FontSize',20)
set(gca,'TickDir','out')
subplot(3,2,4); nhist(all_skew,'proportion'); title('Skewness');
set(gca,'FontSize',20)
set(gca,'TickDir','out')
subplot(3,2,5); nhist(all_kurt,'proportion'); title('Kurtosis');
set(gca,'FontSize',20)
set(gca,'TickDir','out')



% relate to PAC effect

load PAC_centred_power_mean_vector_length


PAC = PAC_centred_power_mean_vector_length;

timeaxis = -2:2;
freqaxis = -10:10;


correct_PAC     = {};
baseline_PAC    = {};
baseline_all    = [];
correct_all     = [];
max_theta       = [];
max_gamma       = [];


for isubject = 1:numel(PAC)
    
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

    max_theta(isubject,:) = PAC{1,isubject}.max_theta;
    max_gamma(isubject,:) = PAC{1,isubject}.max_gamma;
    
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

nSub = numel(PAC);
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
subplot('Position',[.35 .6 .2 .3]);

stats_time = nearest(timeaxis,cfg.latency(1)):nearest(timeaxis,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

% tvalues
tvals = squeeze(Fieldtripstats.stat);
mask = find(squeeze(Fieldtripstats.mask));




for isubject = 1:numel(correct_PAC)
    tmp = squeeze(correct_all(isubject,:,:)-baseline_all(isubject,:,:));


    PAC_to_corr(isubject,:) = mean(tmp(mask));

end


[r_pt, p_pt] = corr(all_peaktrough_ppp', PAC_to_corr, 'rows','complete');
[r_rd, p_rd] = corr(all_risedecay_ppp', PAC_to_corr, 'rows','complete');
[r_sk, p_sk] = corr(all_skew_ppp', PAC_to_corr, 'rows','complete');
[r_ku, p_ku] = corr(all_kurt_ppp', PAC_to_corr, 'rows','complete');

fprintf('Peak–trough corr = %.3f, p=%.3g\n', r_pt, p_pt);
fprintf('Rise–decay corr = %.3f, p=%.3g\n', r_rd, p_rd);
fprintf('Skew corr = %.3f, p=%.3g\n', r_sk, p_sk);
fprintf('Kurtosis corr = %.3f, p=%.3g\n', r_ku, p_ku);



%% plot decoding cross-validated





%% decoding

load perf_dec_cross_validated

perf = perf_dec_cross_validated;





correct_dec_all     = [];
incorrect_dec_all   = [];

sigma = .5; % STD of 2D gaussian smoothing

correct                 = cell(1);
incorrect              = cell(1);
correct_colour         = cell(1);
incorrect_colour        = cell(1);
correct_scene            = cell(1);
incorrect_scene         = cell(1);
baseline_all            = cell(1);
correct_scene_colour     = cell(1);
incorrect_scene_colour   = cell(1);


sigma = .5; % STD of 2D gaussian smoothing

count_trls = 1;
for isubject = 1:size(settings.subjects,1)

    
    % correct
    correct{1,count_trls}                                 = struct;
    correct{1,count_trls}.label                           = {'chan'};
    correct{1,count_trls}.dimord                          = 'chan_freq_time';
    correct{1,count_trls}.freq                            = perf{1,1}.time_train;
    correct{1,count_trls}.time                            = perf{1,1}.time_test;
    correct{1,count_trls}.powspctrm(1,:,:)                = perf{isubject}.correct.accuracy;
   

    correct_dec_all(count_trls,:,:) = imgaussfilt(perf{isubject}.correct.accuracy, sigma);

    % incorrect
    incorrect{1,count_trls}                               = correct{1,count_trls};
    incorrect{1,count_trls}.powspctrm(1,:,:)          = perf{isubject}.incorrect.accuracy;
   
    incorrect_dec_all(count_trls,:,:) = imgaussfilt(perf{isubject}.incorrect.accuracy, sigma);

    % correct colour
    correct_colour{1,count_trls}                          = correct{1,count_trls};
    correct_colour{1,count_trls}.powspctrm(1,:,:)     = perf{isubject}.colour.correct.accuracy;
  

    % incorrect colour
    incorrect_colour{1,count_trls}                        = correct{1,count_trls};
    incorrect_colour{1,count_trls}.powspctrm(1,:,:)   = perf{isubject}.colour.incorrect.accuracy;
   
    % correct scene
    correct_scene{1,count_trls}                           = correct{1,count_trls};
    correct_scene{1,count_trls}.powspctrm(1,:,:)     = perf{isubject}.scene.correct.accuracy;
   
    % incorrect scene
    incorrect_scene{1,count_trls}                        = correct{1,count_trls};
    incorrect_scene{1,count_trls}.powspctrm(1,:,:)   = perf{isubject}.scene.incorrect.accuracy;
   
    % correct scene and colour collapsed
    correct_scene_colour{1,count_trls}                    = correct{1,count_trls};
    correct_scene_colour{1,count_trls}.powspctrm(1,:,:) = mean(cat(1,correct_colour{1,count_trls}.powspctrm,correct_scene{1,count_trls}.powspctrm));
   
    correct_dec_fine_all(count_trls,:,:) = imgaussfilt(mean(cat(1,correct_colour{1,count_trls}.powspctrm,correct_scene{1,count_trls}.powspctrm)), sigma);

    % incorrect scene and colour collapsed
    incorrect_scene_colour{1,count_trls}                  = correct{1,count_trls};
    incorrect_scene_colour{1,count_trls}.powspctrm(1,:,:) = mean(cat(1,incorrect_colour{1,count_trls}.powspctrm,incorrect_scene{1,count_trls}.powspctrm));
  
    incorrect_dec_fine_all(count_trls,:,:) = imgaussfilt(mean(cat(1,incorrect_colour{1,count_trls}.powspctrm,incorrect_scene{1,count_trls}.powspctrm)), sigma);


    % baseline
    baseline_all{1,count_trls}                            = correct{1,count_trls};
    baseline_all{1,count_trls}.powspctrm(1,:,:)           = .5*ones(size(perf{isubject}.correct.accuracy));

    count_trls = count_trls+1;
end


%% Run stats scene and colour collapsed (fine-grained)

timeaxis = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);
timeaxis = correct{1, 1}.time(timeaxis);


freqaxis = correct{1, 1}.freq;

contrast = 'incorrect'; % baseline incorrect

cfg                     = [];
cfg.latency             = [correct{1, 1}.time(1) correct{1, 1}.time(end)]; % ripple time
cfg.frequency           = [-.2   correct{1, 1}.freq(end)]; % encoding time  [freqaxis(1) freqaxis(end)]
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

cfg.numrandomization    = 500;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = numel(correct);
% set up design matrix

% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];
cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
if strcmp(contrast,'baseline')
    [Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, baseline_all{:});
elseif strcmp(contrast,'incorrect')
    [Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, incorrect_scene_colour{:});
end
length(find(Fieldtripstats.mask==1))

%% plot (significant) t vals

stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);

subplot(2,3,2)
imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    tvals);
colormap(flipud(settings.colour_scheme_1))
caxis([-5 5])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',16)
xlabel('ripple time (sec)')
ylabel('encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');

if any(Fieldtripstats.mask(:))
    plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats.mask)))
end



