%% plot all figures for paper

%                  Bernhard Staresina [bernhard.staresina@psy.ox.ac.uk]
%                  Casper Kerren      [kerren@cbs.mpg.de]


% Figure 1 - hypothesis, paradigm
% Figure 2 - electrode coverage and all relevant characteristics
% Figure 3 - decoding and dimensionality + LME results 
% Figure 4 - dPCA


clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults
settings                    = [];
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/MVPA-Light-master'))
addpath(genpath([settings.base_path_castle,'ripple_project_publication_for_replication/main_analyses/Slythm']))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/help_functions'))
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))


settings.colour_scheme_1 = brewermap(30,'RdBu');
settings.colour_scheme_1 = settings.colour_scheme_1;

settings.nu_perm = 4096;

%%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------FIGURE 1--------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------


figure('units','normalized','outerposition',[0 0 1 1]);
gcf = subplot(2,2,1);

%% GA Ripple

load GA_ripple
fs = 1000;

TOI = -6:1/fs:4;
TOI_idx = nearest(TOI,-1):nearest(TOI,1);


m = squeeze(nanmean(GA_ripple(:,TOI_idx)));
s = squeeze(std(GA_ripple(:,TOI_idx)))/sqrt(size(GA_ripple,1));

TOI = TOI(TOI_idx);
boundedline(TOI,m,s,'k','alpha')
plot(TOI,m,'k','LineWidth',2)
axis tight

xlabel('Ripple time (sec)')
ylabel('\muV')
set(gca,'FontSize',20)



%% TF ripple-aligned

load TF_ripple_aligned

hold on
h2 = subplot(2,2,2);
pos2    = get(h2, 'Position');
pos2(1) = .37;
pos2(2) = .82;
pos2(3) = .09;
pos2(4) = .09;
set(h2, 'Position', pos2);

TF = TF_ripple_aligned;
% plot data
to_plot_ripple = [];
for isubject = 1:numel(TF)
    to_plot_ripple(isubject,:,:) = squeeze(mean(TF{isubject}.all));
end

to_plot_ripple = squeeze(mean(to_plot_ripple));

timeaxis   = TF{1}.time;
freqaxis    = TF{1}.freq;

% find peak freq
peak_freq = [];
for isubject = 1:numel(TF)
    [~,idx] = max(squeeze(mean(TF{isubject}.all(:,:,51))));
    peak_freq(isubject,:) = freqaxis(idx);
end

mean(peak_freq)
std(peak_freq)/sqrt(numel(TF))


imagesc(...
    timeaxis,...
    freqaxis,...
    to_plot_ripple);
colormap(hot)
caxis([0 .5])


axis xy
set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Frequency (Hz)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');

%% density ripples

gcf = subplot(2,2,2);

load ripples_characteristics


all_density         = [];

for isubject=1:numel(ripples_characteristics)

   
    all_density(isubject,:)         = ripples_characteristics{isubject}.density;
   
end

N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector


scattercorrect = scatter(.9, all_density(:, 4), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;
scatterincorrect = scatter(1.1, all_density(:, 3), 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color1, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);


[f, xi] = ksdensity(all_density(:, 4));
fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2) 

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(all_density(:, 3)); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;
    
    
  
    line([.9 1.1], [all_density(i, 4), all_density(i, 3)], 'Color', [0.5 0.5 0.5]);
end

[~,p_value_density,~,stats_density] = ttest(all_density(:, 4), all_density(:, 3));
mean_correct = mean(all_density(:, 4));
std_correct = std(all_density(:, 4))/sqrt(12);
mean_incorrect = mean(all_density(:, 3));
std_incorrect = std(all_density(:, 3))/sqrt(12);


hold on
set(gca,'TickDir','out')
ylabel('Ripple density (Hz)')
xticks([.8 1.2])
xticklabels({'AM+', 'AM-'})
xlabel('Trial type')
set(gca,'FontSize',22)
box off




%% ripple rate retrieval STIM ON

gcf = subplot(2,2,3)
load ripple_rate_time_resolved_STIM_ON

time_bnd   = [-.5 3];
time_steps      = .001;

TOI = time_bnd(1):time_steps:time_bnd(2);

TOI = TOI(1:end-1);

TOI_idx = nearest(TOI,-.2):nearest(TOI,3);

TOI = TOI(TOI_idx);


ripple_latency_correct = ripple_rate_time_resolved_STIM_ON.retrieval.correct;
ripple_latency_incorrect = ripple_rate_time_resolved_STIM_ON.retrieval.incorrect;


ripple_rate_correct = [];
ripple_rate_incorrect = [];
for isubject = 1:numel(ripple_latency_correct)
    ripple_rate_correct(isubject,:) = ripple_latency_correct{1, isubject}.rate(TOI_idx);
    ripple_rate_incorrect(isubject,:) = ripple_latency_incorrect{1, isubject}.rate(TOI_idx);
end


% stats
correct_incorrect = {};

for isubject = 1:numel(ripple_latency_correct)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = TOI;

    correct_incorrect{1}.individual(isubject,1,:)   = ripple_rate_correct(isubject,:);
    correct_incorrect{1}.dimord              = 'subj_chan_time';



end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);


correct_incorrect{2}                = correct_incorrect{1};

for isubject =1:numel(ripple_latency_correct)

    correct_incorrect{2}.individual(isubject,1,:)           = ripple_rate_incorrect(isubject,:);
   
end



cfg                     = [];
cfg.latency             = [-.2 3];

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


% plot significant vals

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

d = squeeze(correct_incorrect{1}.individual);

m       = nanmean(d);
max_m   = max(m);
s       = nanstd(d)./sqrt(size(d,1));
max_s   = max(s);

boundedline(correct_incorrect{1, 1}.time(stats_time)  ,m,s, 'cmap', settings.colour_scheme_1(23,:),'alpha');
plot(correct_incorrect{1, 1}.time(stats_time)  ,m,'Color',settings.colour_scheme_1(23,:),'linewidth',2);
hold on

d = squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(correct_incorrect{1, 1}.time(stats_time)  ,m,s, 'cmap', settings.colour_scheme_1(7,:),'alpha');
plot(correct_incorrect{1, 1}.time(stats_time)  ,m,'Color',settings.colour_scheme_1(7,:),'linewidth',2);
hold on

sigline   = nan(1,numel(correct_incorrect{1}.time));
sigline(stats_time(Fieldtripstats.mask==1)) = max_m+max_s;

plot(correct_incorrect{1}.time(stats_time),sigline,'k','linewidth',4);

set(gca,'FontSize',20,'FontName','Arial')
xlabel('Time after cue (sec)')
ylabel('Z-scored ripple density (Hz)')
set(gca,'TickDir','out')


axis([-1 1 -.1 max_m+max_s+.02])
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])

hline(0)

%% ripple rate RT

gcf = subplot(2,2,4)

load ripple_rate_time_resolved_RT


ripple_latency_correct = ripple_rate_time_resolved_RT.retrieval.correct;
ripple_latency_incorrect = ripple_rate_time_resolved_RT.retrieval.incorrect;

ripple_rate_correct = [];
ripple_rate_incorrect = [];
for isubject = 1:numel(ripple_latency_correct)
    ripple_rate_correct(isubject,:) = ripple_latency_correct{1, isubject}.rate;
    ripple_rate_incorrect(isubject,:) = ripple_latency_incorrect{1, isubject}.rate;
end


TOI = ripple_latency_correct{1, isubject}.time(1,1:end-1);

correct_incorrect = {};

for isubject = 1:numel(ripple_latency_correct)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = TOI;

    correct_incorrect{1}.individual(isubject,1,:)   = ripple_rate_correct(isubject,:);
    correct_incorrect{1}.dimord              = 'subj_chan_time';



end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);


correct_incorrect{2}                = correct_incorrect{1};

for isubject =1:numel(ripple_latency_correct)

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


% plot significant vals



d = squeeze(correct_incorrect{1}.individual);

m       = nanmean(d);
max_m   = max(m);
s       = nanstd(d)./sqrt(size(d,1));
max_s   = max(s);

boundedline(correct_incorrect{1, 1}.time  ,m,s, 'cmap', settings.colour_scheme_1(23,:),'alpha');
plot(correct_incorrect{1, 1}.time  ,m,'Color',settings.colour_scheme_1(23,:),'linewidth',2);
hold on

d = squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(correct_incorrect{1, 1}.time  ,m,s, 'cmap', settings.colour_scheme_1(7,:),'alpha');
plot(correct_incorrect{1, 1}.time  ,m,'Color',settings.colour_scheme_1(7,:),'linewidth',2);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = max_m+max_s;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',20,'FontName','Arial')
xlabel('Time pre and post button press (sec)')
ylabel('Z-scored ripple density (Hz)')
set(gca,'TickDir','out')


axis([-1 1 -.1 max_m+max_s+.02])
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])


%%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------FIGURE 2--------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------


figure('units','normalized','outerposition',[0 0 1 1]);
gcf = subplot('Position',[.06 .55 .22 .3]);

%% DECODING ripple-locked - fine

load tvals_fine_enc_ripple
load mask_t_vals_fine_enc_ripple

timeaxis = linspace(-1,1, size(tvals_fine_enc_ripple,2));
freqaxis = linspace(-.2,3, size(tvals_fine_enc_ripple,1));


imagesc(...
    timeaxis,...
    freqaxis,...
    tvals_fine_enc_ripple);
colormap(flipud(settings.colour_scheme_1))
caxis([-5 5])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
ylabel(hcb, 't-values', 'FontSize', 20); % Add colorbar label


mask = mask_t_vals_fine_enc_ripple;
hold on;
alphaData = zeros(size(mask_t_vals_fine_enc_ripple)); 
alphaData(~mask) = 0.5; 
overlay = imagesc(timeaxis, freqaxis, zeros(size(mask_t_vals_fine_enc_ripple)), 'AlphaData', alphaData);
set(overlay, 'AlphaData', alphaData);
set(overlay, 'CData', zeros(size(mask_t_vals_fine_enc_ripple))); 
set(overlay, 'AlphaDataMapping', 'none');


plot_contour(timeaxis,freqaxis,double(mask_t_vals_fine_enc_ripple))

% to later correlate
time_dec = linspace(-.5,3,351);
idx_time = nearest(time_dec,-.2):nearest(time_dec,3);


load correct_dec_fine_all
load incorrect_dec_fine_all


dec_to_corr = [];
for isubject = 1:size(correct_dec_fine_all,1)
    tmp = squeeze((correct_dec_fine_all(isubject,idx_time,:)-incorrect_dec_fine_all(isubject,idx_time,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr(isubject,:) = mean(mean(tmp));
end



%% Dimensionality - ripple-locked

hold on
gcf = subplot('Position',[.35 .55 .22 .3]);


load perf_dimensionality

perf = perf_dimensionality;

RT_correct = cellfun(@mean, perf{1}.RT.correct);
RT_max_correct = cellfun(@max, perf{1}.RT.correct);
RT_min_correct = cellfun(@min, perf{1}.RT.correct);

RT_incorrect = cellfun(@mean, perf{1}.RT.incorrect);
RT_max_incorrect = cellfun(@max, perf{1}.RT.incorrect);
RT_min_incorrect = cellfun(@min, perf{1}.RT.incorrect);

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

set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Dimensionality')
set(gca,'TickDir','out')
xticks([-.95 0 .95])
xticklabels({'-1', '0', '1'})
yticks([2 2.5 3 3.5 4])
yticklabels({'2', '2.5', '3', '3.5' '4'})

ylim([2 4])
vline(0)


xlim([cfg.latency(1), cfg.latency(end)])

d = squeeze(correct_incorrect{1}.individual-correct_incorrect{2}.individual);
dim_to_correlate = mean(d(:,sigline==3.7),2);
d = squeeze(correct_incorrect{1}.individual);
dim_corr_to_correlate = mean(d(:,sigline==3.7),2);

% correlation dimensionality and decoding
[rho_dim_dec,p_dim_dec] = corr(dim_to_correlate,dec_to_corr,'type','spearman');

cohens_d = 2 * rho_dim_dec / sqrt(1 - rho_dim_dec^2);


fprintf('correlation dimensionality and decoding, rho %.3f, p-value %.3f, cohens d %.3f \n',rho_dim_dec,p_dim_dec,cohens_d);

%% RT correlation

hold on

h1 = subplot('Position', [.75 .76 .08 .08]);
pos = get(h1, 'Position'); 
pos(1) = .37;
set(h1, 'Position', pos);


d = squeeze(correct_incorrect{1}.individual);
RT = RT_correct;

dimensionality_change = mean(d(:, sigline == 3.7), 2);

h_scat = scatter(RT, dimensionality_change, 'filled', 'MarkerFaceColor', '#0072BD', 'SizeData', 75);
xlabel('RT', 'FontSize', 14);
ylabel('Dimensionality', 'FontSize', 14);
set(gca, 'FontSize', 14)
grid on;
box on;
hold on;

% Fit a linear regression line
p = polyfit(RT, dimensionality_change, 1);
ylim([2 5])
xlim([1 3])


x_limits = xlim;
x_fit = linspace(x_limits(1), x_limits(2), 100);
y_fit = polyval(p, x_fit);
plot(x_fit, y_fit, 'k-', 'LineWidth', 1.5);

set(gcf, 'Color', 'w'); % Set background color of the figure to white

yresid = dimensionality_change - polyval(p, RT);
SSresid = sum(yresid.^2);
SStotal = (length(dimensionality_change) - 1) * var(dimensionality_change);
rsq = 1 - SSresid / SStotal;
disp(['R-squared: ', num2str(rsq)]);

[rho_rt_dim, p_rt_dim] = corr(RT', dimensionality_change, 'type', 'spearman');

cohens_d = 2 * rho_rt_dim / sqrt(1 - rho_rt_dim^2);

fprintf('correlation dimensionality and RT, rho %.3f, p-value %.3f, cohens d %.3f \n', rho_rt_dim, p_rt_dim, cohens_d);

%% Dimensionality - ripple-locked LME

hold on
gcf = subplot('Position',[.65 .55 .22 .3]);
settings.colour_scheme_1 = brewermap(30,'RdBu');

load perf_dimensionality_LME
perf = perf_dimensionality_LME;


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



subjects = unique(dataTable.SubjectID);

% Preallocate
meanCorrectFirst  = zeros(numel(subjects), 1);
meanCorrectSecond = zeros(numel(subjects), 1);
meanIncorrectFirst  = zeros(numel(subjects), 1);
meanIncorrectSecond = zeros(numel(subjects), 1);

for i = 1:numel(subjects)
    subj = subjects(i);
    
    % Logical indexing for each condition
    idx_CF = dataTable.SubjectID == subj & dataTable.Condition == 1 & dataTable.Half == 1;
    idx_CS = dataTable.SubjectID == subj & dataTable.Condition == 1 & dataTable.Half == 2;
    idx_IF = dataTable.SubjectID == subj & dataTable.Condition == 2 & dataTable.Half == 1;
    idx_IS = dataTable.SubjectID == subj & dataTable.Condition == 2 & dataTable.Half == 2;

    meanCorrectFirst(i)  = mean(dataTable.Accuracy(idx_CF));
    meanCorrectSecond(i) = mean(dataTable.Accuracy(idx_CS));
    meanIncorrectFirst(i)  = mean(dataTable.Accuracy(idx_IF));
    meanIncorrectSecond(i) = mean(dataTable.Accuracy(idx_IS));
end


% Correct: First vs Second Half
[~, p_correct, ~, stats_correct] = ttest(meanCorrectFirst, meanCorrectSecond);
fprintf('Correct trials: t(%d) = %.2f, p = %.4f\n', stats_correct.df, stats_correct.tstat, p_correct);

% Incorrect: First vs Second Half
[~, p_incorrect, ~, stats_incorrect] = ttest(meanIncorrectFirst, meanIncorrectSecond);
fprintf('Incorrect trials: t(%d) = %.2f, p = %.4f\n', stats_incorrect.df, stats_incorrect.tstat, p_incorrect);



%% PLOT
uniqueSubjects = unique(dataTable.SubjectID);

means = zeros(numel(uniqueSubjects), 2, 2);  % Dimensions: subject x condition x half
sems = zeros(numel(uniqueSubjects), 2, 2);

for isubject = 1:numel(uniqueSubjects)
    for condition = 1:2
        for half = 1:2
            subset = dataTable(dataTable.SubjectID == uniqueSubjects(isubject) & ...
                               dataTable.Condition == condition & ...
                               dataTable.Half == half, :);
            means(isubject, condition, half) = mean(subset.Accuracy);
            sems(isubject, condition, half) = std(subset.Accuracy) / sqrt(height(subset));
        end
    end
end


conditionLabels = {'Correct', 'Incorrect'};
halfLabels = {'First Half', 'Second Half'};


mean_correct    = squeeze(means(:,1,:));
mean_incorrect  = squeeze(means(:,2,:));



data_to_plot = [mean_correct(:,2)-mean_correct(:,1),mean_incorrect(:,2)-mean_incorrect(:,1)];



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



[a,p_value_lme.correct,c,d] = ttest(data_to_plot(:,1),0);
stats_lme.correct = d.tstat;
[a,p_value_lme.incorrect,c,d] = ttest(data_to_plot(:,2),0);
stats_lme.incorrect = d.tstat;

hold on
set(gca, 'TickDir', 'out');
ylabel(sprintf('Dimensionality change\n(post-pre ripple)')); 
xticks([.8 1.2])
xticklabels({'AM+', 'AM-'}); 
xlabel('Trial type');
set(gca, 'FontSize', 20);
ylim([-.4 .7])
yticks([-.4 -.2, 0, .2, .4, .6])
yticklabels({'-.4', '-.2 ', '0', '.2', '.4', '.6'})

hline(0)
box off


range=axis;


%%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------------------FIGURE 3--------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------

figure('units','normalized','outerposition',[0 0 1 1]);

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

mean(max_theta)
std(max_theta)./sqrt(length(max_theta))

mean(max_gamma)
std(max_gamma)./sqrt(length(max_gamma))



subplot('Position',[.10 .6 .2 .3]);


scattercorrect = scatter(.9, max_gamma, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', settings.colour_scheme_1(23,:), 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);

hold on
[f, xi] = ksdensity(max_gamma,'Bandwidth',4); 
fill(.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', settings.colour_scheme_1(23,:))
hold on;

set(gca, 'TickDir', 'out');
set(gca, 'FontSize', 20);
hold on;


scattercorrect = scatter(1.3, max_theta, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', settings.colour_scheme_1(7,:), 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on
[f, xi] = ksdensity(max_theta,'Bandwidth',4); 
fill(1.4 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', settings.colour_scheme_1(7,:))


set(gca, 'TickDir', 'out');
ylabel('Frequency (Hz)'); 
set(gca, 'FontSize', 20);


set(gca, 'TickDir', 'out');
ylim([-6 100])

xticks([.8 1.4])
xticklabels({'Gamma', 'Theta'})





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
mask = double(squeeze(Fieldtripstats.mask));

imagesc(...
    timeaxis,...
    freqaxis,...
    tvals);
colormap(flipud(settings.colour_scheme_1))
caxis([-4 4])


axis xy
set(gca,'FontSize',20)
xlabel('Frequency phase (Hz)')
ylabel('Frequency power (Hz)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
ylim([min(freqaxis(stats_freq)) max(freqaxis(stats_freq))])
xlim([min(timeaxis(stats_time)) max(timeaxis(stats_time))])

hold on;
alphaData = zeros(size(tvals)); 
alphaData(~mask) = 0.5; 
overlay = imagesc(timeaxis(stats_time), freqaxis(stats_freq), zeros(size(tvals)), 'AlphaData', alphaData);
set(overlay, 'AlphaData', alphaData);
set(overlay, 'CData', zeros(size(tvals))); 
set(overlay, 'AlphaDataMapping', 'none');


%% PAC across time

subplot('Position',[.6 .7 .2 .2]);
load PAC_across_time_centred_power

PAC = PAC_across_time_centred_power;


freqaxis = -5:5;

timeaxis = PAC{1, 1}.time;

correct_PAC = {};
baseline_PAC = {};

correct_all = [];
baseline_all = [];

for isubject = 1:numel(PAC)
    
    % correct
    correct_PAC{1,isubject}                                 = struct;
    correct_PAC{1,isubject}.label                           = {'chan'};
    correct_PAC{1,isubject}.dimord                          = 'subj_chan_freq_time';
    correct_PAC{1,isubject}.freq                            = freqaxis;
    correct_PAC{1,isubject}.time                            = timeaxis;
    
    tmp = squeeze(mean(PAC{1,isubject}.correct(:,:,1:3),3));
    tmp = tmp';

    correct_PAC{1,isubject}.powspctrm(1,1,:,:)                = tmp;
    correct_all(isubject,:,:) = tmp;
    % baseline
    baseline_PAC{1,isubject}                                = correct_PAC{1,isubject};
    tmp = squeeze(mean(PAC{1,isubject}.permuted(:,:,1:3),3));
    tmp = tmp';

    baseline_PAC{1,isubject}.powspctrm(1,1,:,:)               = tmp;
    baseline_all(isubject,:,:) = tmp; 
    
end

% Run stats

cfg                     = [];
cfg.latency             = [timeaxis(10) timeaxis(end)];%[4 8];
cfg.frequency           = [freqaxis(1) freqaxis(end)]; %[10 20]
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

stats_time = nearest(timeaxis,cfg.latency(1)):nearest(timeaxis,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

% tvalues
tvals = squeeze(Fieldtripstats.stat);

imagesc(...
    timeaxis(stats_time),...
    freqaxis(stats_freq),...
    tvals);
colormap(flipud(settings.colour_scheme_1))
caxis([-4 4])


axis xy
set(gca,'FontSize',20)
% xlabel('Ripple time (sec)')
ylabel('Frequency power (Hz)')
set(gca,'TickDir','out')
yticks([-5 0 5])
yticklabels({'-5', '0', '5'})
xticks([-.5 0 .5 1])
xticklabels({' ', ' ', ' ', ' '})
hcb = colorbar('Location','EastOutside');


% overlay with transparent mask to highlight cluster
mask = double(squeeze(Fieldtripstats.mask));
hold on;
alphaData = zeros(size(tvals)); 
alphaData(~mask) = 0.5; 
overlay = imagesc(timeaxis(stats_time), freqaxis(stats_freq), zeros(size(tvals)), 'AlphaData', alphaData);
set(overlay, 'AlphaData', alphaData);
set(overlay, 'CData', zeros(size(tvals))); 
set(overlay, 'AlphaDataMapping', 'none');


%% plot PAC across time in lineplot

settings.colour_scheme_1 = brewermap(30,'RdBu');

subplot('Position',[.6 .6 .18 .07]);

timeaxis = PAC{1, 1}.time;

d = squeeze(mean(correct_all(:,:,:),2))-squeeze(mean(baseline_all(:,:,:),2));


m = mean(d);

PAC_to_cross_corr = d;

s = std(d)./sqrt(size(d,1));

boundedline(timeaxis,m',s','k','cmap', settings.colour_scheme_1(26,:));
plot(timeaxis,m,'k','linewidth',2);
hold on

stats_time = nearest(timeaxis,cfg.latency(1)):nearest(timeaxis,cfg.latency(2));

sigline   = nan(1,numel(timeaxis));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(squeeze(mean(Fieldtripstats.mask))>0)) = 0;

plot(timeaxis,sigline,'r','linewidth',4);

set(gca,'FontSize',20,'FontName','Arial')
xlabel('Ripple time (sec)')
ylabel('PAC')
xticks([-.99 0 .99])
xticklabels({'-1', '0', '1'})
set(gca,'TickDir','out')

axis tight
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])

% correlate dimensionality and PAC across time

PAC_across_time_to_corr = mean(d(:,sigline==0),2);

[rho_dim_PAC_time, p_dim_PAC_time] = corr(dim_corr_to_correlate, PAC_across_time_to_corr,'type','spearman');

cohens_d = 2 * rho_dim_PAC_time / sqrt(1 - rho_dim_PAC_time^2);

fprintf('correlation PAC across time and dimensionality, rho %.2f, p-value %.2f, cohens d %.2f \n',rho_dim_PAC_time,p_dim_PAC_time,cohens_d);


% check gamma power increase correlates with PAC
load perf_gamma_power
gamma_all = [];
for isubject = 1:numel(PAC)
    
 
    gamma_all(isubject,:) = squeeze(mean(mean(perf_gamma_power{isubject}.gamma(:,sigline==0))));
    % baseline
   
    
end

[rho_dim_PAC_gamma, p_dim_PAC_gamma] = corr(gamma_all, PAC_across_time_to_corr,'type','spearman');




%% correlation PAC and dimensionality

gcf = subplot('Position',[.85 .8 .1 .1]);
scatter(PAC_across_time_to_corr, dim_corr_to_correlate, 'filled', 'MarkerFaceColor', '#0072BD');
xlabel('PAC', 'FontSize', 20);
ylabel('Dimensionality', 'FontSize', 20);
xlim([min(PAC_across_time_to_corr)-.001 max(PAC_across_time_to_corr)+.001])
ylim([min(dim_corr_to_correlate)-.3 max(dim_corr_to_correlate)+.3])
set(gca,'FontSize',20)
grid on;
box on;
hold on;

% Fit a linear regression line
p = polyfit(PAC_across_time_to_corr, dim_corr_to_correlate, 1);
f = polyval(p, PAC_across_time_to_corr);

plot(PAC_across_time_to_corr, f, 'r-', 'LineWidth', 1.5);

% outlier_index = 6;
% scatter(RT(outlier_index), dimensionality_change(outlier_index), 150, 'ro', 'LineWidth', 2);


% Add legend
% legend('Data', 'Linear Fit', 'Location', 'best');

% Customize the plot appearance
set(gcf, 'Color', 'w');    % Set background color of the figure to white


% Can see that subject 6 is off the radar a bit.

yfit = polyval(p, PAC_across_time_to_corr);
yresid = dim_corr_to_correlate - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(dim_corr_to_correlate)-1) * var(dim_corr_to_correlate);
rsq = 1 - SSresid/SStotal;
disp(['R-squared: ', num2str(rsq)]);




%% cross-correlate PAC across time and dimensionality

% do it per participant
% since they have different sample rates, let us resample the PAC_across
% time to match that of dimensionality

subplot('Position',[.85 .6 .1 .1]);


load perf_dimensionality

perf = perf_dimensionality;



correct_incorrect = {};

% Store behavioral accuracy time series
for isubject = 1:numel(perf)
    correct_incorrect{1}.label = {'Channels'};
    correct_incorrect{1}.time = perf{1,1}.time_train;
    correct_incorrect{1}.individual(isubject,1,:) = perf{isubject}.correct.accuracy;
end
correct_incorrect{1}.dimord = 'subj_chan_time';
correct_incorrect{1}.avg = squeeze(correct_incorrect{1}.individual);

% Store incorrect accuracy
correct_incorrect{2} = correct_incorrect{1};
for isubject = 1:numel(perf)
    correct_incorrect{2}.individual(isubject,1,:) = perf{isubject}.incorrect.accuracy;
end
correct_incorrect{2}.avg = squeeze(correct_incorrect{2}.individual);

% Extract time window [-1, 1]
idx_range = nearest(perf{1,1}.time_train, -1):nearest(perf{1,1}.time_train, 1);
t_dim = perf{1,1}.time_train(idx_range); % behavioral time vector

% Accuracy difference signal
% dim_to_cross_corr = squeeze(correct_incorrect{1}.individual(:,:,idx_range))-squeeze(correct_incorrect{2}.individual(:,:,idx_range));
dim_to_cross_corr = squeeze(correct_incorrect{1}.individual(:,:,idx_range));

% PAC difference signal (already computed externally)
PAC_to_cross_corr = squeeze(mean(correct_all(:,:,:),2)) - squeeze(mean(baseline_all(:,:,:),2));

% PAC time vector
t_PAC = PAC_across_time_centred_power{1,1}.time;

% Sampling rates
fs_PAC = 1 / mean(diff(t_PAC));
fs_dim = 1 / mean(diff(t_dim));

% Common time vector over overlapping range
t_start = max(t_dim(1), t_PAC(1));
t_end   = min(t_dim(end), t_PAC(end));
fs_target = max(fs_PAC, fs_dim);
t_common = t_start:1/fs_target:t_end;

% Interpolation
nSubjects = size(dim_to_cross_corr, 1);
PAC_to_cross_corr_interp = nan(nSubjects, length(t_common));
dim_to_cross_corr_interp = nan(nSubjects, length(t_common));

for isubject = 1:nSubjects
    PAC_to_cross_corr_interp(isubject, :) = interp1(t_PAC, PAC_to_cross_corr(isubject, :), t_common, 'linear', NaN);
    dim_to_cross_corr_interp(isubject, :) = interp1(t_dim, dim_to_cross_corr(isubject, :), t_common, 'linear', NaN);
end

% Remove time points with NaNs in any subject
valid_mask = all(~isnan(PAC_to_cross_corr_interp), 1) & all(~isnan(dim_to_cross_corr_interp), 1);
PAC_to_cross_corr_interp = PAC_to_cross_corr_interp(:, valid_mask);
dim_to_cross_corr_interp = dim_to_cross_corr_interp(:, valid_mask);
t_common = t_common(valid_mask);




r_perm              = [];
lags_perm           = [];
r_all               = [];
lags_all            = [];
zscore_all          = [];
num_permutations    = 500;
for isubject = 1:size(dim_to_cross_corr,1)
    isubject
    [r,lags] = xcorr(dim_to_cross_corr_interp(isubject,:),PAC_to_cross_corr_interp(isubject,:));
    r_all(isubject,:) = r;
    lags_all(isubject,:) = lags;


    for iperm = 1:num_permutations
       
        %         rand_trials         = randsample(size(dim_to_cross_corr,2), size(dim_to_cross_corr,2));
        %         [r_tmp,lags_tmp]    = xcorr(dim_to_cross_corr(isubject,rand_trials),PAC_to_cross_corr(isubject,:));

        surrogate = generate_surrogate_iaaft(dim_to_cross_corr_interp(isubject,:), 'M', 1, 'detrend', false, 'verbose', false);
        [r_tmp, lags_tmp] = xcorr(surrogate, PAC_to_cross_corr_interp(isubject,:));



        r_perm(isubject,iperm,:)     = r_tmp;
        lags_perm(isubject,iperm,:)  = lags_tmp;

        
    end

    m_perm = mean(squeeze(r_perm(isubject,:,:)));
    std_perm = std(squeeze(r_perm(isubject,:,:)));

    zscore_all(isubject,:) = (r_all(isubject,:)-m_perm)./std_perm;

end

mean_r_perm = squeeze(mean(mean(r_perm, 1), 2));
sem_r_perm = squeeze(std(mean(r_perm, 1), 0, 2));

mean_r_all = mean(r_all, 1);
sem_r_all = std(r_all, 1)/sqrt(size(r_all,1));

% plot the zscore

zscore_mean = (mean_r_all'-mean_r_perm)./sem_r_perm;

time_lags = lags * (1 / fs_target);



mean_r_all = mean(r_all, 1);                          % size: [1 x time]
r_perm_avg = squeeze(mean(r_perm, 1));                % size: [nPermutations x time]

m_perm = mean(r_perm_avg, 1);                         % mean across permutations
std_perm = std(r_perm_avg, 0, 1);                     % std across permutations

% Step 3: Z-score the real average relative to permuted
z_real = (mean_r_all - m_perm) ./ std_perm;

% Step 4: Threshold and find clusters
z_thresh = 1.96;  % or use prctile of null
binary_real = z_real > z_thresh;
[clusts_real, n_clusts_real] = bwlabel(binary_real);

% Step 5: Get cluster-level stats
real_cluster_stats = zeros(1, n_clusts_real);
for c = 1:n_clusts_real
    real_cluster_stats(c) = sum(z_real(clusts_real == c));
end
max_real_cluster_stat = max(real_cluster_stats);

% Step 6: Get max cluster stat for each permutation
max_cluster_stats_perm = zeros(1, size(r_perm_avg,1));
for iperm = 1:size(r_perm_avg,1)
    z_perm = (r_perm_avg(iperm,:) - m_perm) ./ std_perm;

    binary_perm = z_perm > z_thresh;
    [clusts, num_clusts] = bwlabel(binary_perm);

    cluster_stats = zeros(1, num_clusts);
    for c = 1:num_clusts
        cluster_stats(c) = sum(z_perm(clusts == c));
    end

    if num_clusts > 0
        max_cluster_stats_perm(iperm) = max(cluster_stats);
    else
        max_cluster_stats_perm(iperm) = 0;
    end
end

% Step 7: Compute p-value
p_value = mean(max_cluster_stats_perm >= max_real_cluster_stat)


% Plot empirical data using stem plot


% shadedErrorBar(time_lags, mean_r_all, sem_r_all, {'--s','Color', settings.colour_scheme_1(23, :)},.7);
% hold on;
% plot(time_lags, mean_r_all,'LineWidth', 3, 'Color', settings.colour_scheme_1(23, :));
% 
% % Plot mean and SEM for permuted data
% 
% shadedErrorBar(time_lags, mean_r_perm, sem_r_perm);
% 
% 
% grid on;
% 
% set(gca,'FontSize',20,'FontName','Arial')
% set(gca,'TickDir','out')
% 
% axis tight
% 
% 
% figure;
% 
% plot(time_lags, zscore_mean,'LineWidth', 3, 'Color', settings.colour_scheme_1(23, :));
% 
% hold on;
% yline(1.96, '--r');   % upper 95% threshold
% yline(-1.96, '--r');  % lower 95% threshold
% 
% sig_pos = find(zscore_mean > 1.96);
% sig_neg = find(zscore_mean < -1.96);
% 
% plot(time_lags(sig_pos), zscore_mean(sig_pos), 'ro');  % tau = lag axis
% plot(time_lags(sig_neg), zscore_mean(sig_neg), 'bo');
% 
% axis([-.5 1 -1 3])
% vline(0)
% hline(0)
% 
% 
% figure;
% 
% m = mean(zscore_all);
% s = mean(zscore_all)./sqrt(size(zscore_all,1));
% 
% boundedline(time_lags  ,m,s, 'cmap', settings.colour_scheme_1(23,:),'alpha');
% plot(time_lags  ,m,'Color',settings.colour_scheme_1(23,:),'linewidth',2);
% hold on
% 
% axis([-.5 1 -.5 .5])
% vline(0)
% hline(0)


% permutation test

correct_incorrect = {};

% Real data (cross-correlation)
correct_incorrect{1}.label = {'chan'};
correct_incorrect{1}.time = time_lags;
correct_incorrect{1}.dimord = 'subj_chan_time';
correct_incorrect{1}.individual = zeros(size(r_all,1),1,length(time_lags));

% Baseline (shuffled) data
correct_incorrect{2} = correct_incorrect{1};

for isubject = 1:size(r_all,1)
    correct_incorrect{1}.individual(isubject,1,:) = r_all(isubject,:);
    correct_incorrect{2}.individual(isubject,1,:) = mean(squeeze(r_perm(isubject,:,:)), 1); % mean of perms
end

% Averages (not strictly necessary unless plotting)
correct_incorrect{1}.avg = squeeze(mean(correct_incorrect{1}.individual,1));
correct_incorrect{2}.avg = squeeze(mean(correct_incorrect{2}.individual,1));


cfg = [];
cfg.latency             = [time_lags(147) time_lags(end)];
cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'analytic';
cfg.correctm            = 'no';
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.05;
cfg.tail                = 1;
cfg.correcttail         = 'alpha';
cfg.neighbours          = [];         % no spatial clustering
cfg.minnbchan           = 0;          % only one "channel"
cfg.avgovertime         = 'no';
cfg.avgoverchan         = 'no';
cfg.computecritval      = 'yes';
cfg.numrandomization    = 1000;
cfg.clusterstatistic    = 'maxsum';
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'individual';

nSub = size(r_all,1);
design = zeros(2, 2*nSub);
design(1,:) = [1:nSub 1:nSub];
design(2,:) = [ones(1,nSub) 2*ones(1,nSub)];

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

% run stats
[Fieldtripstats] = ft_timelockstatistics(cfg, correct_incorrect{:});

length(find(Fieldtripstats.mask==1))



d = squeeze(correct_incorrect{1}.individual)-squeeze(correct_incorrect{2}.individual);

m       = nanmean(d);
max_m   = max(m);
s       = nanstd(d)./sqrt(size(d,1));
max_s   = max(s);

boundedline(correct_incorrect{1, 1}.time  ,m,s, 'cmap', settings.colour_scheme_1(23,:),'alpha');
plot(correct_incorrect{1, 1}.time  ,m,'Color',settings.colour_scheme_1(23,:),'linewidth',2);
hold on


stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline(stats_time(Fieldtripstats.mask==1)) = max_m+max_s;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',20,'FontName','Arial')
xlabel('Time lags')
ylabel('Cross-correlation')
set(gca,'TickDir','out')


axis([-.5 1 -.1 max_m+max_s+.02])
vline(0)
hline(0)


[~, max_peak] = max(d(:,stats_time)');

adj_time = correct_incorrect{1}.time(stats_time);

average_peak = mean(adj_time(max_peak))
