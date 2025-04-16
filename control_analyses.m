ng %% DECODING STIM ONSET

clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults


settings                    = [];
settings.base_path          = '/Users/kerrenadmin/Desktop/Postdoc/Project_1/';
settings.base_path_castle   = '/Users/kerrenadmin/Desktop/Other_projects/Dimensionality_ripples_Casper_and_Bernhard/'; % '/castles/nr/projects/w/wimberm-ieeg-compute/';

addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422/external'))
addpath([settings.base_path_castle,'ripple_project_publication_for_replication/subfunctions'])
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/help_functions'))
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/plotting'))
addpath(genpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/toolbox'))
settings.colour_scheme_1 = brewermap(30,'RdBu');
settings.colour_scheme_1 = settings.colour_scheme_1;

settings.nu_perm = 4096;

settings.performance = 1; % 1 correct vs incorrect, 2 = only correct

%%
load perf_decoding_ripple_aligned
perf = perf_decoding_ripple_aligned;


load perf_decoding_ripple_aligned_correct_shuffled
perf_correct_shuff = perf_decoding_ripple_aligned_correct_shuffled;

load perf_decoding_ripple_aligned_circshift_corr_shuffled
perf_circshift = perf_decoding_ripple_aligned_circshift_corr_shuffled;

load perf_decoding_ripple_aligned_circshift_corr_shuffled_backward
perf_circshift_backward = perf_decoding_ripple_aligned_circshift_corr_shuffled_backward;

correct                     = cell(1,numel(perf_correct_shuff));
incorrect                   = cell(1,numel(perf_correct_shuff));
correct_colour              = cell(1,numel(perf_correct_shuff));
incorrect_colour            = cell(1,numel(perf_correct_shuff));
correct_scene               = cell(1,numel(perf_correct_shuff));
incorrect_scene             = cell(1,numel(perf_correct_shuff));
baseline_all                = cell(1,numel(perf_correct_shuff));
baseline_z                  = cell(1,numel(perf_correct_shuff));
z_all                       = cell(1,numel(perf_correct_shuff));
shuff_all                   = cell(1,numel(perf_correct_shuff));
circshift_all               = cell(1,numel(perf_correct_shuff));
correct_scene_colour        = cell(1,numel(perf_correct_shuff));
incorrect_scene_colour      = cell(1,numel(perf_correct_shuff));
correct_dec_fine_all_zvalue = [];
shuffled_dec_fine_all       = [];
corr_circshift_all          = [];
correct_dec_fine_all        = [];
incorrect_dec_fine_all      = [];

sigma = .5; % STD of 2D gaussian smoothing

timeaxis_idx = nearest(perf{1,1}.time_test, -1):nearest(perf{1,1}.time_test, 1);
timeaxis = perf{1,1}.time_test(timeaxis_idx);

freqaxis_idx = nearest(perf{1,1}.time_train, -.2):nearest(perf{1,1}.time_train, 3);
freqaxis = perf{1,1}.time_train(freqaxis_idx);

for isubject = 1:numel(perf_correct_shuff)

     % correct
    correct{1,isubject}                                 = struct;
    correct{1,isubject}.label                           = {'chan'};
    correct{1,isubject}.dimord                          = 'chan_freq_time';
    correct{1,isubject}.freq                            = freqaxis;
    correct{1,isubject}.time                            = timeaxis;
    correct{1,isubject}.powspctrm(1,:,:)                = perf{isubject}.correct.accuracy(freqaxis_idx,timeaxis_idx);
    
    % incorrect
    incorrect{1,isubject}                               = correct{1,isubject};
    incorrect{1,isubject}.powspctrm(1,:,:)              = perf{isubject}.incorrect.accuracy(freqaxis_idx,timeaxis_idx);
    
    % correct colour
    correct_colour{1,isubject}                          = correct{1,isubject};
    correct_colour{1,isubject}.powspctrm(1,:,:)         = perf{isubject}.colour.correct.accuracy(freqaxis_idx,timeaxis_idx);
    
    % incorrect colour
    incorrect_colour{1,isubject}                        = correct{1,isubject};
    incorrect_colour{1,isubject}.powspctrm(1,:,:)       = perf{isubject}.colour.incorrect.accuracy(freqaxis_idx,timeaxis_idx);
    
    % correct scene
    correct_scene{1,isubject}                           = correct{1,isubject};
    correct_scene{1,isubject}.powspctrm(1,:,:)          = perf{isubject}.scene.correct.accuracy(freqaxis_idx,timeaxis_idx);
    
    % incorrect scene
    incorrect_scene{1,isubject}                         = correct{1,isubject};
    incorrect_scene{1,isubject}.powspctrm(1,:,:)        = perf{isubject}.scene.incorrect.accuracy(freqaxis_idx,timeaxis_idx);
    
    % correct scene and colour collapsed
    correct_scene_colour{1,isubject}                    = correct{1,isubject};
    correct_scene_colour{1,isubject}.powspctrm(1,:,:)   = mean(cat(1,correct_colour{1,isubject}.powspctrm,correct_scene{1,isubject}.powspctrm));
    
    correct_dec_fine_all(isubject,:,:) = mean(cat(1,correct_colour{1,isubject}.powspctrm,correct_scene{1,isubject}.powspctrm));


    % incorrect scene and colour collapsed
    incorrect_scene_colour{1,isubject}                  = correct{1,isubject};
    incorrect_scene_colour{1,isubject}.powspctrm(1,:,:) = mean(cat(1,incorrect_colour{1,isubject}.powspctrm,incorrect_scene{1,isubject}.powspctrm));
    
    incorrect_dec_fine_all(isubject,:,:) = mean(cat(1,incorrect_colour{1,isubject}.powspctrm,incorrect_scene{1,isubject}.powspctrm));

    % baseline
    baseline_all{1,isubject}                            = correct{1,isubject};
    baseline_all{1,isubject}.powspctrm(1,:,:)           = .5*ones(size(perf{isubject}.correct.accuracy(freqaxis_idx,timeaxis_idx)));
    

    % shuffled
         
    tmp_shuff = (perf_correct_shuff{isubject}.colour.correct.accuracy(:,freqaxis_idx,timeaxis_idx)+perf_correct_shuff{isubject}.scene.correct.accuracy(:,freqaxis_idx,timeaxis_idx))/2;
    tmp_emp = squeeze(mean(cat(1,correct_colour{1,isubject}.powspctrm,correct_scene{1,isubject}.powspctrm)));

    m_shuff     = squeeze(mean(tmp_shuff));
    std_shuff   = squeeze(std(tmp_shuff));

    correct_dec_fine_all_zvalue(isubject,:,:) = (tmp_emp-m_shuff)./std_shuff;

    z_all{1,isubject}                            = correct{1,isubject};
    z_all{1,isubject}.powspctrm(1,:,:)           = (tmp_emp-m_shuff)./std_shuff;

    % baseline
    shuff_all{1,isubject}                            = correct{1,isubject};
    shuff_all{1,isubject}.powspctrm(1,:,:)           = squeeze(mean(tmp_shuff));

    shuffled_dec_fine_all(isubject,:,:) = squeeze(mean(tmp_shuff));

    circshift_all{1,isubject}                            = correct{1,isubject};
    circshift_all{1,isubject}.powspctrm(1,:,:)           = (perf_circshift{isubject}.colour.correct.accuracy(freqaxis_idx,timeaxis_idx)+perf_circshift{isubject}.scene.correct.accuracy(freqaxis_idx,timeaxis_idx)+perf_circshift_backward{isubject}.colour.correct.accuracy(freqaxis_idx,timeaxis_idx)+perf_circshift_backward{isubject}.scene.correct.accuracy(freqaxis_idx,timeaxis_idx))/4; 

    corr_circshift_all(isubject,:,:) = (perf_circshift{isubject}.colour.correct.accuracy(freqaxis_idx,timeaxis_idx)+perf_circshift{isubject}.scene.correct.accuracy(freqaxis_idx,timeaxis_idx)+perf_circshift_backward{isubject}.colour.correct.accuracy(freqaxis_idx,timeaxis_idx)+perf_circshift_backward{isubject}.scene.correct.accuracy(freqaxis_idx,timeaxis_idx))/4;

    baseline_z{1,isubject}                            = correct{1,isubject};
    baseline_z{1,isubject}.powspctrm(1,:,:)           = zeros(size(perf{isubject}.correct.accuracy(freqaxis_idx,timeaxis_idx)));
    
    
end



%% stats

timeaxis = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);
timeaxis = correct{1, 1}.time(timeaxis);



contrast = 'shuffled'; % baseline shuffled shuff_base, corr_circshift

cfg                     = [];
cfg.latency             = [timeaxis(1) timeaxis(end)]; % ripple time
cfg.frequency           = [-.2 freqaxis(end)]; % encoding time  [freqaxis(1) freqaxis(end)]
cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'no'; % 'no', cluster, bonferroni, fdr, holm;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 1;
cfg.correcttail         = 'alpha'; % alpha prob no
cfg.neighbours          = [];
cfg.minnbchan           = 0;
cfg.computecritval      = 'yes';

cfg.numrandomization    = 500;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = numel(perf_correct_shuff);
% set up design matrix

% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];
cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats

[Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, shuff_all{:});

cfg.correctm            = 'cluster'; % 'no', cluster, bonferroni, fdr, holm;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.numrandomization    = 500;%1000;%'all';
cfg.clustertail         = cfg.tail;

[Fieldtripstats_incorrect] = ft_freqstatistics(cfg, correct_scene_colour{:}, incorrect_scene_colour{:});
length(find(Fieldtripstats_incorrect.mask==1))



[Fieldtripstats_shuff] = ft_freqstatistics(cfg, shuff_all{:}, baseline_all{:});



contrast = 'corr_circshift'; % baseline shuffled shuff_base, corr_circshift

cfg                     = [];
cfg.latency             = [timeaxis(1) timeaxis(end)]; % ripple time
cfg.frequency           = [-.2 freqaxis(end)]; % encoding time  [freqaxis(1) freqaxis(end)]
cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'no'; % 'no', cluster, bonferroni, fdr, holm;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 1;
cfg.correcttail         = 'alpha'; % alpha prob no
cfg.neighbours          = [];
cfg.minnbchan           = 0;
cfg.computecritval      = 'yes';

cfg.numrandomization    = 500;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = numel(perf_correct_shuff);
% set up design matrix

% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];
cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
[Fieldtripstats_circ] = ft_freqstatistics(cfg, correct_scene_colour{:}, circshift_all{:});


cfg.latency             = [timeaxis(1) timeaxis(end)]; % ripple time
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'no'; % 'no', cluster, bonferroni, fdr, holm;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 1;
cfg.clustertail         = cfg.tail;
[Fieldtripstats_z] = ft_freqstatistics(cfg, z_all{:}, baseline_z{:});
any(Fieldtripstats_z.mask(:))
length(find(Fieldtripstats_z.mask==1))


length(find(Fieldtripstats_circ.mask==1))


%% plot (significant) t vals


freq_idx = nearest(Fieldtripstats.freq,0):nearest(Fieldtripstats.freq,3);

emp_data    = squeeze(mean(correct_dec_fine_all(:,freq_idx,:)-incorrect_dec_fine_all(:,freq_idx,:),2));
shuff_data  = squeeze(mean(correct_dec_fine_all(:,freq_idx,:)-shuffled_dec_fine_all(:,freq_idx,:),2));
adj_data    = squeeze(mean(correct_dec_fine_all(:,freq_idx,:)-corr_circshift_all(:,freq_idx,:),2));

load mask_t_vals_fine_enc_ripple

figure('units','normalized','outerposition',[0 0 1 1]);


stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));


gcf = subplot(2,2,1);
ax1 = gcf;
ax1.Position(2) = ax1.Position(2) + 0.03; 


imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    squeeze(mean(shuffled_dec_fine_all)));
colormap(flipud(settings.colour_scheme_1))
caxis([0.45 .65])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title('shuffled trials');

% mask = mask_t_vals_fine_enc_ripple;
% hold on;
% alphaData = zeros(size(mask_t_vals_fine_enc_ripple)); 
% alphaData(~mask) = 0.3; 
% overlay = imagesc(timeaxis, freqaxis, zeros(size(mask_t_vals_fine_enc_ripple)), 'AlphaData', alphaData);
% set(overlay, 'AlphaData', alphaData);
% set(overlay, 'CData', zeros(size(mask_t_vals_fine_enc_ripple))); 
% set(overlay, 'AlphaDataMapping', 'none');

plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(mask_t_vals_fine_enc_ripple)))



dec_to_corr_real        = [];
dec_to_corr_shuff       = [];
dec_to_corr_circshift   = [];
for isubject = 1:size(correct_dec_fine_all,1)
    tmp = squeeze((correct_dec_fine_all(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_real(isubject,:) = mean(mean(tmp));
    tmp = squeeze((shuffled_dec_fine_all(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_shuff(isubject,:) = mean(mean(tmp));
    tmp = squeeze((corr_circshift_all(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_circshift(isubject,:) = mean(mean(tmp));
end


gcf = subplot(2,2,3);
ax1 = gcf;
ax1.Position(2) = ax1.Position(2) + 0.03; 


imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    squeeze(mean(corr_circshift_all)));
colormap(flipud(settings.colour_scheme_1))
caxis([0.45 .65])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title('adjacent trials');


plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(mask_t_vals_fine_enc_ripple)))




gcf = subplot(2,2,2);
ax1 = gcf;
ax1.Position(2) = ax1.Position(2) + 0.03;  % Shift downwards


% data_dec = {};
% data_dec{1, 1} = dec_to_corr_real;
% data_dec{2, 1} = dec_to_corr_shuff;
% 
% 
% glm_fit = [];
% yfit = [];
% 
% for isubject = 1:numel(perf)
%     x = [data_dec{1, 1}(isubject,1),1;data_dec{2, 1}(isubject,1),2];
%     [bb,dev,~] = glmfit(x(:,2),x(:,1),'normal');
%     glm_fit(isubject,:) = bb(2,1);
%     yfit(isubject,:) = polyval([bb(2,1),bb(1,1)],[1,2]);
% end
% 
% % gcf = subplot(3,3,9);
% 
% 
% h   = rm_raincloud(data_dec, [settings.colour_scheme_1(15,:)],0,'ks',[],settings.colour_scheme_1,yfit);
% h.p{1, 1}.FaceColor = settings.colour_scheme_1(23,:);
% h.s{1, 1}.MarkerFaceColor = settings.colour_scheme_1(23,:);
% h.m(1, 1).MarkerFaceColor = settings.colour_scheme_1(23,:);
% h.p{2, 1}.FaceColor = settings.colour_scheme_1(7,:);
% h.s{2, 1}.MarkerFaceColor = settings.colour_scheme_1(7,:);
% h.m(2, 1).MarkerFaceColor = settings.colour_scheme_1(7,:);


N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector

scattercorrect = scatter(.9, dec_to_corr_real, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;
scatterincorrect = scatter(1.1, dec_to_corr_shuff, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color1, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);


[f, xi] = ksdensity(dec_to_corr_real); 
fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2)

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(dec_to_corr_shuff); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;
    
    
    
    line([.9 1.1], [dec_to_corr_real, dec_to_corr_shuff], 'Color', [0.5 0.5 0.5]);
end



hold on
set(gca,'TickDir','out')
ylabel('Deoding perf. (arb. val.)')

xticks([.8 1.2])
xticklabels({'AM+', 'AM+ shuffled'})
xlabel('Condition')
set(gca,'FontSize',20)
box off

range=axis;



gcf = subplot(2,2,4);
ax1 = gcf;
ax1.Position(2) = ax1.Position(2) + 0.03;  % Shift downwards



N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector

scattercorrect = scatter(.9, dec_to_corr_real, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;
scatterincorrect = scatter(1.1, dec_to_corr_circshift, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color1, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);


[f, xi] = ksdensity(dec_to_corr_real); 
fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2)

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(dec_to_corr_circshift); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;
    
    
    
    line([.9 1.1], [dec_to_corr_real, dec_to_corr_circshift], 'Color', [0.5 0.5 0.5]);
end


hold on
set(gca,'TickDir','out')
ylabel('Deoding perf. (arb. val.)')

xticks([.8 1.2])
xticklabels({'AM+', 'AM+ adjacent'})
xlabel('Condition')
set(gca,'FontSize',20)
box off

range=axis;






stats = [];
[~,pvalue.dec.real_shuff, ~,d] = ttest(dec_to_corr_real,dec_to_corr_shuff);
stats.dec.real_shuff = d.tstat;
[~,pvalue.dec.real_circshift, ~,d] = ttest(dec_to_corr_real,dec_to_corr_circshift);
stats.dec.real_circshift = d.tstat;

pvalue.dec
stats.dec


%% Dimensionality


load perf_dimensionality

perf = perf_dimensionality;


load perf_dimensionality_shuffled_correct

perf_shuf_corr = perf_dimensionality_shuffled_correct;

load perf_dimensionality_circshift_shuffled_correct

perf_circshift_corr = perf_dimensionality_circshift_shuffled_correct;

load perf_dimensionality_circshift_shuffled_correct_backwards

perf_circshift_corr_back = perf_dimensionality_circshift_shuffled_correct_backwards;


correct_incorrect = {};
xlimits = nearest(perf{1,1}.time_train, -1):nearest(perf{1,1}.time_train, 1);

correct_incorrect = {}; % correct_incorrect{1} = empirical, {2} = shuffled. {3} = zvalue, {4} = empirical incorrect, {5} = circshift correct

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train(xlimits);

    correct_incorrect{1}.individual(isubject,1,:)           = perf{1,isubject}.correct.accuracy(xlimits);

    correct_incorrect{1}.dimord              = 'subj_chan_time';

end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};
correct_incorrect{3}                = correct_incorrect{1};
correct_incorrect{4}                = correct_incorrect{1};
correct_incorrect{5}                = correct_incorrect{1};
correct_incorrect{6}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = mean(perf_shuf_corr{1,isubject}.correct.accuracy(:,xlimits));

    tmp = perf{isubject}.correct.accuracy(xlimits);
    m =  mean(perf_shuf_corr{isubject}.correct.accuracy(:,xlimits));
    s = std(perf_shuf_corr{isubject}.correct.accuracy(:,xlimits));

    correct_incorrect{3}.individual(isubject,1,:)   = (tmp-m)/s;
    correct_incorrect{4}.individual(isubject,1,:)   = perf{1,isubject}.incorrect.accuracy(xlimits);

    correct_incorrect{5}.individual(isubject,1,:)           = (perf_circshift_corr{1,isubject}.correct.accuracy(xlimits)+perf_circshift_corr_back{1,isubject}.correct.accuracy(xlimits))/2;
  
    correct_incorrect{6}.individual(isubject,1,:)   = zeros(size(perf{1,isubject}.incorrect.accuracy(xlimits)));
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);
correct_incorrect{3}.avg            = squeeze(correct_incorrect{3}.individual);
correct_incorrect{4}.avg            = squeeze(correct_incorrect{4}.individual);
correct_incorrect{5}.avg            = squeeze(correct_incorrect{5}.individual);
correct_incorrect{6}.avg            = squeeze(correct_incorrect{6}.individual);




%% stats

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
[Fieldtripstats_corr_inc] = ft_timelockstatistics(cfg, correct_incorrect{1,[1,4]});
length(find(Fieldtripstats_corr_inc.mask))

cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.tail                = 0;
cfg.clustertail         = cfg.tail;

[Fieldtripstats_shuff] = ft_timelockstatistics(cfg, correct_incorrect{1,[1,2]});
length(find(Fieldtripstats_shuff.mask))

cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.tail                = 1;
cfg.clustertail         = cfg.tail;

[Fieldtripstats_shuf_inc] = ft_timelockstatistics(cfg, correct_incorrect{1,[2,4]});
length(find(Fieldtripstats_shuf_inc.mask))

cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.tail                = 1;
cfg.clustertail         = cfg.tail;

[Fieldtripstats_circshift] = ft_timelockstatistics(cfg, correct_incorrect{1,[1,5]});
length(find(Fieldtripstats_circshift.mask))

cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.tail                = 0;
cfg.clustertail         = cfg.tail;

[Fieldtripstats_circshift_vs_inc] = ft_timelockstatistics(cfg, correct_incorrect{1,[5,4]});
length(find(Fieldtripstats_circshift_vs_inc.mask))


cfg.tail                = 1;
cfg.clustertail         = cfg.tail;
[Fieldtripstats_z] = ft_timelockstatistics(cfg, correct_incorrect{1,[3,6]});
length(find(Fieldtripstats_z.mask))
%% plot significant values

figure('units','normalized','outerposition',[0 0 1 1]);

d = squeeze(correct_incorrect{1}.individual);

dim_to_correlate_corr = d;

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline_emp   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline_emp(stats_time(Fieldtripstats_corr_inc.mask==1)) = 2;

dim_to_correlate_corr = mean(dim_to_correlate_corr(:,sigline_emp==2),2);


subplot(2,2,1)

d = squeeze(correct_incorrect{2}.individual);

dim_to_correlate_shuff_corr = mean(d(:,sigline_emp == 2), 2);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));


boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(7,:), 'alpha');
hold on
plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 1);
sig_times_2 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 2);


if ~isempty(sig_times_1)
    fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

if ~isempty(sig_times_2)
    fill([sig_times_2(1) sig_times_2(end) sig_times_2(end) sig_times_2(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');  
end

set(gca, 'FontSize', 20, 'FontName', 'Arial')
xlabel('Ripple time (s)')
ylabel('Dimensionality')
title('Shuffled')
set(gca, 'TickDir', 'out')
ylim([2 3.7])
vline(0)


subplot(2,2,2)



N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector

scattercorrect = scatter(.9, dim_to_correlate_corr, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;
scatterincorrect = scatter(1.1, dim_to_correlate_shuff_corr, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color1, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);


[f, xi] = ksdensity(dim_to_correlate_corr); 
fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2)

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(dim_to_correlate_shuff_corr); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;

    line([.9 1.1], [dim_to_correlate_corr, dim_to_correlate_shuff_corr], 'Color', [0.5 0.5 0.5]);
end


hold on
set(gca,'TickDir','out')
ylabel('Dimensionality')

xticks([.8 1.2])
xticklabels({'AM+', 'AM+ shuffled'})
xlabel('Trial type')
set(gca,'FontSize',20)
box off

range=axis;


subplot(2,2,3)

d = squeeze(correct_incorrect{5}.individual);

dim_to_correlate_ajdacent_corr = mean(d(:,sigline_emp == 2), 2);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));


boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(7,:), 'alpha');
hold on
plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);


sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 1);
sig_times_2 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 2);

if ~isempty(sig_times_1)
    fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end


if ~isempty(sig_times_2)
    fill([sig_times_2(1) sig_times_2(end) sig_times_2(end) sig_times_2(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none'); 
end

set(gca, 'FontSize', 20, 'FontName', 'Arial')
xlabel('Ripple time (s)')
ylabel('Dimensionality')
title('Adjacent')
set(gca, 'TickDir', 'out')
ylim([2 3.7])
vline(0)





subplot(2,2,4)




data_to_plot = [dim_to_correlate_corr,dim_to_correlate_ajdacent_corr];


N = 12;
color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector

scattercorrect = scatter(.9, dim_to_correlate_corr, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);
hold on;
scatterincorrect = scatter(1.1, dim_to_correlate_ajdacent_corr, 75, 'MarkerEdgeColor', 'black',...
    'MarkerFaceColor', color1, 'MarkerFaceAlpha',.3,...
    'MarkerEdgeAlpha',.8);


[f, xi] = ksdensity(dim_to_correlate_corr); 
fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color2)

% Second kernel density estimation and fill, mirrored
[f, xi] = ksdensity(dim_to_correlate_ajdacent_corr); 
fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
    'EdgeColor', 'none', 'FaceColor', color1) 

for i = 1:N
    hold on;

    line([.9 1.1], [dim_to_correlate_corr, dim_to_correlate_ajdacent_corr], 'Color', [0.5 0.5 0.5]);
end


hold on
set(gca,'TickDir','out')
ylabel('Dimensionality')

xticks([.8 1.2])
xticklabels({'AM+', 'AM+ adjacent'})
xlabel('Trial type')
set(gca,'FontSize',20)
box off

range=axis;


%%

[~,pvalue.dim.shuff_corr, ~,d] = ttest(dim_to_correlate_corr,dim_to_correlate_shuff_corr)
stats.dim.shuff_corr = d.tstat
[~,pvalue.dim.ajd_corr, ~,d] = ttest(dim_to_correlate_corr,dim_to_correlate_ajdacent_corr)
stats.dim.ajd_corr = d.tstat


