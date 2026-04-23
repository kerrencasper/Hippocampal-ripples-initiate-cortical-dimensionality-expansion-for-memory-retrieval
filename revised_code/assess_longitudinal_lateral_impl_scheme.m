%% plot all figures for paper

%                  Casper Kerren      [kerren@cbs.mpg.de]


% Figure 1 - Paradigm, implantation scheme, hypothesis
% Figure 2 - Hippocampal ripple density increases during successful memory retrieval
% Figure 3 - Target memory decoding and dimensionality transformation are locked to ripple events.
% Figure 4 - Phase-amplitude coupling following ripple events are related to dimensionality expansion


clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

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
addpath('/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/scripts_to_publish')


settings.colour_scheme_1 = brewermap(30,'RdBu');
settings.colour_scheme_1 = settings.colour_scheme_1;

settings.nu_perm = 4096;

%% DECODING ripple-locked - fine

load '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/scripts_to_publish/to_plot_for_paper/perf_decoding_ripple_aligned'


perf = perf_decoding_ripple_aligned;

timeaxis = linspace(-1,1, size(perf{1, 1}.colour.correct.accuracy,2));
freqaxis = linspace(-.5,3, size(perf{1, 1}.colour.correct.accuracy,1));


subjects = 1:12;

correct                 = cell(1,numel(subjects));
incorrect               = cell(1,numel(subjects));
correct_colour          = cell(1,numel(subjects));
incorrect_colour        = cell(1,numel(subjects));
correct_scene           = cell(1,numel(subjects));
incorrect_scene         = cell(1,numel(subjects));
baseline_all            = cell(1,numel(subjects));
correct_scene_colour    = cell(1,numel(subjects));
incorrect_scene_colour  = cell(1,numel(subjects));

correct_dec_all     = [];
incorrect_dec_all   = [];

sigma = .5; % STD of 2D gaussian smoothing


for isubject = 1:12

    % correct
    correct{1,isubject}                                 = struct;
    correct{1,isubject}.label                           = {'chan'};
    correct{1,isubject}.dimord                          = 'chan_freq_time';
    correct{1,isubject}.freq                            = perf{1,1}.time_train;
    correct{1,isubject}.time                            = perf{1,1}.time_test;
    correct{1,isubject}.powspctrm(1,:,:)                = perf{isubject}.correct.accuracy;
    if size(perf{isubject}.correct.accuracy,1) == numel(freqaxis) && size(perf{isubject}.correct.accuracy,2) == numel(timeaxis)
        correct{1,isubject}.powspctrm(1,:,:)            = perf{isubject}.correct.accuracy;
    else
        correct{1,isubject}.powspctrm(1,:,:)            = nan(numel(freqaxis),numel(timeaxis));
    end

    correct_dec_all(isubject,:,:) = imgaussfilt(perf{isubject}.correct.accuracy, sigma);

    % incorrect
    incorrect{1,isubject}                               = correct{1,isubject};
    if size(perf{isubject}.incorrect.accuracy,1) == numel(freqaxis) && size(perf{isubject}.incorrect.accuracy,2) == numel(timeaxis)
        incorrect{1,isubject}.powspctrm(1,:,:)          = perf{isubject}.incorrect.accuracy;
    else
        incorrect{1,isubject}.powspctrm(1,:,:)          = nan(numel(freqaxis),numel(timeaxis));
    end

    incorrect_dec_all(isubject,:,:) = imgaussfilt(perf{isubject}.incorrect.accuracy, sigma);

    % correct colour
    correct_colour{1,isubject}                          = correct{1,isubject};
    if size(perf{isubject}.colour.correct.accuracy,1) == numel(freqaxis) && size(perf{isubject}.colour.correct.accuracy,2) == numel(timeaxis)
        correct_colour{1,isubject}.powspctrm(1,:,:)     = perf{isubject}.colour.correct.accuracy;
    else
        correct_colour{1,isubject}.powspctrm(1,:,:)     = nan(numel(freqaxis),numel(timeaxis));
    end

    % incorrect colour
    incorrect_colour{1,isubject}                        = correct{1,isubject};
    if size(perf{isubject}.colour.incorrect.accuracy,1) == numel(freqaxis) && size(perf{isubject}.colour.incorrect.accuracy,2) == numel(timeaxis)
        incorrect_colour{1,isubject}.powspctrm(1,:,:)   = perf{isubject}.colour.incorrect.accuracy;
    else
        incorrect_colour{1,isubject}.powspctrm(1,:,:)   = nan(numel(freqaxis),numel(timeaxis));
    end

    % correct scene
    correct_scene{1,isubject}                           = correct{1,isubject};
    if size(perf{isubject}.scene.correct.accuracy,1) == numel(freqaxis) && size(perf{isubject}.scene.correct.accuracy,2) == numel(timeaxis)
        correct_scene{1,isubject}.powspctrm(1,:,:)     = perf{isubject}.scene.correct.accuracy;
    else
        correct_scene{1,isubject}.powspctrm(1,:,:)     = nan(numel(freqaxis),numel(timeaxis));
    end

    % incorrect scene
    incorrect_scene{1,isubject}                        = correct{1,isubject};
    if size(perf{isubject}.scene.incorrect.accuracy,1) == numel(freqaxis) && size(perf{isubject}.scene.incorrect.accuracy,2) == numel(timeaxis)
        incorrect_scene{1,isubject}.powspctrm(1,:,:)   = perf{isubject}.scene.incorrect.accuracy;
    else
        incorrect_scene{1,isubject}.powspctrm(1,:,:)   = nan(numel(freqaxis),numel(timeaxis));
    end

    % correct scene and colour collapsed
    correct_scene_colour{1,isubject}                    = correct{1,isubject};
    if ...
            size(perf{isubject}.colour.correct.accuracy,1) == numel(freqaxis) && size(perf{isubject}.colour.correct.accuracy,2) == numel(timeaxis) && ...
            size(perf{isubject}.scene.correct.accuracy,1) == numel(freqaxis) && size(perf{isubject}.scene.correct.accuracy,2) == numel(timeaxis)

        correct_scene_colour{1,isubject}.powspctrm(1,:,:) = mean(cat(1,correct_colour{1,isubject}.powspctrm,correct_scene{1,isubject}.powspctrm));
    else
        correct_scene_colour{1,isubject}.powspctrm(1,:,:) = nan(numel(freqaxis),numel(timeaxis));
    end

    correct_dec_fine_all(isubject,:,:) = imgaussfilt(mean(cat(1,correct_colour{1,isubject}.powspctrm,correct_scene{1,isubject}.powspctrm)), sigma);

    % incorrect scene and colour collapsed
    incorrect_scene_colour{1,isubject}                  = correct{1,isubject};
    if ...
            size(perf{isubject}.colour.incorrect.accuracy,1) == numel(freqaxis) && size(perf{isubject}.colour.incorrect.accuracy,2) == numel(timeaxis) && ...
            size(perf{isubject}.scene.incorrect.accuracy,1) == numel(freqaxis) && size(perf{isubject}.scene.incorrect.accuracy,2) == numel(timeaxis)

        incorrect_scene_colour{1,isubject}.powspctrm(1,:,:) = mean(cat(1,incorrect_colour{1,isubject}.powspctrm,incorrect_scene{1,isubject}.powspctrm));
    else
        incorrect_scene_colour{1,isubject}.powspctrm(1,:,:) = nan(numel(freqaxis),numel(timeaxis));
    end

    incorrect_dec_fine_all(isubject,:,:) = imgaussfilt(mean(cat(1,incorrect_colour{1,isubject}.powspctrm,incorrect_scene{1,isubject}.powspctrm)), sigma);


    % baseline
    baseline_all{1,isubject}                            = correct{1,isubject};
    baseline_all{1,isubject}.powspctrm(1,:,:)           = .5*ones(size(perf{isubject}.correct.accuracy));

end

%% Run stats scene and colour collapsed (fine-grained)
figure('units','normalized','outerposition',[0 0 1 1]);
gcf = subplot('Position',[.06 .55 .22 .3]);

timeaxis = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);
timeaxis = correct{1, 1}.time(timeaxis);

% longi

longi_scheme_correct = correct_scene_colour(:,1:8);
lat_scheme_correct = correct_scene_colour(:,9:12);

longi_scheme_incorrect = incorrect_scene_colour(:,1:8);
lat_scheme_incorrect = incorrect_scene_colour(:,9:12);


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

nSub = numel(longi_scheme_correct);
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
    [Fieldtripstats] = ft_freqstatistics(cfg, longi_scheme_correct{:}, longi_scheme_incorrect{:});
end
length(find(Fieldtripstats.mask==1))

% plot (significant) t vals

stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);

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
title('Longitudinal scheme')

if any(Fieldtripstats.mask(:))
    plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats.mask)))
end



%% lat

hold on
gcf = subplot('Position',[.35 .55 .22 .3]);

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

nSub = numel(lat_scheme_correct);
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
    [Fieldtripstats] = ft_freqstatistics(cfg, lat_scheme_correct{:}, lat_scheme_incorrect{:});
end
length(find(Fieldtripstats.mask==1))

% plot (significant) t vals

stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);

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
title('Lateral scheme')

if any(Fieldtripstats.mask(:))
    plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats.mask)))
end


%% Dimensionality - ripple-locked

hold on
gcf = subplot('Position',[.06 .15 .22 .3]);


load '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/scripts_dimensionality/scripts_to_publish/to_plot_for_paper/perf_dimensionality'

perf = perf_dimensionality;

RT_correct = cellfun(@mean, perf{1}.RT.correct);
RT_max_correct = cellfun(@max, perf{1}.RT.correct);
RT_min_correct = cellfun(@min, perf{1}.RT.correct);

RT_incorrect = cellfun(@mean, perf{1}.RT.incorrect);
RT_max_incorrect = cellfun(@max, perf{1}.RT.incorrect);
RT_min_incorrect = cellfun(@min, perf{1}.RT.incorrect);



%% longi
correct_incorrect = {};

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(isubject,1,:)   = perf{isubject}.correct.accuracy;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end

correct_incorrect{1}.individual = correct_incorrect{1, 1}.individual(1:8,:,:);
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = perf{isubject}.incorrect.accuracy;
  
end

correct_incorrect{2}.individual = correct_incorrect{2}.individual(1:8,:,:);

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

cfg.numrandomization    = settings.nu_perm;%'all';

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


% plot significant vals

d = squeeze(correct_incorrect{1}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(correct_incorrect{1, 1}.time  ,m,s, 'cmap', settings.colour_scheme_1(23,:));
plot(correct_incorrect{1, 1}.time  ,m,'k','linewidth',3);
hold on

d = squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(correct_incorrect{1, 1}.time  ,m,s, 'cmap', settings.colour_scheme_1(7,:));
plot(correct_incorrect{1, 1}.time  ,m,'k','linewidth',3);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 3.7;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',16)
xlabel('Ripple time (sec)')
ylabel('Dimensionality')
set(gca,'TickDir','out')
xticks([-.95 0 .95])
xticklabels({'-1', '0', '1'})
yticks([2 2.5 3 3.5 4])
yticklabels({'2', '2.5', '3', '3.5' '4'})
title('Longitudinal scheme')

ylim([2 4])
vline(0)


xlim([cfg.latency(1), cfg.latency(end)])


%% lat

hold on
gcf = subplot('Position',[.35 .15 .22 .3]);




correct_incorrect = {};

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(isubject,1,:)   = perf{isubject}.correct.accuracy;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end

correct_incorrect{1}.individual = correct_incorrect{1, 1}.individual(9:12,:,:);
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = perf{isubject}.incorrect.accuracy;
  
end

correct_incorrect{2}.individual = correct_incorrect{2}.individual(9:12,:,:);

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

cfg.numrandomization    = settings.nu_perm;%'all';

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

boundedline(correct_incorrect{1, 1}.time  ,m,s, 'cmap', settings.colour_scheme_1(23,:));
plot(correct_incorrect{1, 1}.time  ,m,'k','linewidth',3);
hold on

d = squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(correct_incorrect{1, 1}.time  ,m,s, 'cmap', settings.colour_scheme_1(7,:));
plot(correct_incorrect{1, 1}.time  ,m,'k','linewidth',3);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = 3.7;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',16)
xlabel('Ripple time (sec)')
ylabel('Dimensionality')
set(gca,'TickDir','out')
xticks([-.95 0 .95])
xticklabels({'-1', '0', '1'})
yticks([2 2.5 3 3.5 4])
yticklabels({'2', '2.5', '3', '3.5' '4'})
title('Lateral scheme')
ylim([2 4])
vline(0)


xlim([cfg.latency(1), cfg.latency(end)])
