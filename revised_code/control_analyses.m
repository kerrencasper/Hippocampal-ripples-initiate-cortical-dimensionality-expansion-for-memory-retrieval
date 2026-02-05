%% Controls decoding and dimensionality

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

load 'perf_dec_shuffle_time_of_ripple'
perf_shuffle_retr_time = perf_shuffle_time_of_ripple;

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
shuff_retrieval             = cell(1,numel(perf_correct_shuff));
correct_dec_fine_all_zvalue = [];
shuffled_dec_fine_all       = [];
corr_circshift_all          = [];
correct_dec_fine_all        = [];
incorrect_dec_fine_all      = [];
all_suffled                 = [];
z_value_shuff_retrieval     = [];
shuffled_retrieval_all      = [];

     

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

    % compared to shuffled time retrieval

    mean_shuffle_all = squeeze(mean((perf_shuffle_retr_time{isubject}.colour.correct+perf_shuffle_retr_time{isubject}.scene.correct)/2));
    std_shuffle_all = squeeze(std((perf_shuffle_retr_time{isubject}.colour.correct+perf_shuffle_retr_time{isubject}.scene.correct)));

    shuffled_retrieval_all(isubject,:,:) = mean_shuffle_all(freqaxis_idx,:);

    z_value_shuff_retrieval(isubject,:,:) = (squeeze(correct_dec_fine_all(isubject,:,:))-mean_shuffle_all(freqaxis_idx,:))./std_shuffle_all(freqaxis_idx,:);

    shuff_retrieval{1,isubject}                          = correct{1,isubject};
    shuff_retrieval{1,isubject}.powspctrm(1,:,:)         = (squeeze(correct_dec_fine_all(isubject,:,:))-mean_shuffle_all(freqaxis_idx,:))./std_shuffle_all(freqaxis_idx,:);
    

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
    
    all_suffled(isubject,:,:,:) = tmp_shuff;
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

figure;
z_on_group = [];
for isubject = 1:numel(perf_correct_shuff)
    subplot(4,3,isubject)
   m_shuff     = squeeze(mean(all_suffled(isubject,:,:,:)));
    std_shuff   = squeeze(std(all_suffled(isubject,:,:,:)));


    z_on_group(isubject,:,:) = (squeeze(correct_dec_fine_all(isubject,:,:))-m_shuff)./std_shuff;

    imagesc(timeaxis,freqaxis,squeeze(z_on_group(isubject,:,:)))
    axis xy
    colormap(jet)
    caxis([0 2])
    vline(0)
    hline(0)


end

m_shuff     = squeeze(mean(mean(all_suffled)));
std_shuff   = squeeze(std(mean(all_suffled)));
z_on_group_group = (squeeze(mean(correct_dec_fine_all))-m_shuff)./std_shuff;
imagesc(z_on_group_group)
axis xy


%% Stelzer method

load mask_t_vals_fine_enc_ripple
mask_idx = find(mask_t_vals_fine_enc_ripple);

Nsub = numel(perf_correct_shuff);

emp_sub = zeros(Nsub,1);

for isubject = 1:Nsub
    emp_map = squeeze(mean(cat(1,correct_colour{1,isubject}.powspctrm,correct_scene{1,isubject}.powspctrm)));

    emp_sub(isubject) = mean(emp_map(mask_idx));  % masked mean (no zeros)
end

emp_group_fine_grained_correct = mean(emp_sub);


B = 10000;                  % iterations for group-level null
K = 1000;                   % surrogates per subject (trials)

null_group = nan(B,1);

for b = 1:B
    sub_vals = nan(Nsub,1);

    for s = 1:Nsub
        k = randi(K);  % pick one surrogate index for this subject

        shuf_map = ( perf_correct_shuff{s}.colour.correct.accuracy(k,freqaxis_idx,timeaxis_idx) + ...
                     perf_correct_shuff{s}.scene.correct.accuracy(k,freqaxis_idx,timeaxis_idx) )/2;

        shuf_map = squeeze(shuf_map);
        sub_vals(s) = mean(shuf_map(mask_idx));  % masked mean
    end

    null_group(b) = mean(sub_vals);
end

null_group_fine_grained = null_group;

p_emp_vs_random_correct = (sum(abs(null_group_fine_grained - mean(null_group_fine_grained)) >= abs(emp_group_fine_grained_correct - mean(null_group_fine_grained))) + 1) / (B + 1);

z_emp_vs_random_correct = (emp_group_fine_grained_correct - mean(null_group_fine_grained)) / std(null_group_fine_grained);


% now against shuffled retrieval

B = 10000;                  % iterations for group-level null
K = 40;                   % surrogates per subject (trials)

null_group = nan(B,1);

for b = 1:B
    sub_vals = nan(Nsub,1);

    for s = 1:Nsub
        k = randi(K);  % pick one surrogate index for this subject

        shuf_map = ( perf_shuffle_retr_time{s}.colour.correct(k,freqaxis_idx,:) + ...
                     perf_shuffle_retr_time{s}.scene.correct(k,freqaxis_idx,:) )/2;
        shuf_map = squeeze(shuf_map);
        sub_vals(s) = mean(shuf_map(mask_idx));  % masked mean
    end

    null_group(b) = mean(sub_vals);
end

null_group_retrieval = null_group;

p_emp_vs_retrieval_correct = (sum(abs(null_group - mean(null_group)) >= abs(emp_group_fine_grained_correct - mean(null_group))) + 1) / (B + 1);

z_emp_vs_retrieval_correct = (emp_group_fine_grained_correct - mean(null_group_retrieval)) / std(null_group_retrieval);



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

% empirical vs shuffled

[Fieldtripstats] = ft_freqstatistics(cfg, correct_scene_colour{:}, shuff_all{:});



% empirical correct vs empirical incorrect
cfg.correctm            = 'cluster'; % 'no', cluster, bonferroni, fdr, holm;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.numrandomization    = 500;%1000;%'all';
cfg.clustertail         = cfg.tail;

[Fieldtripstats_incorrect] = ft_freqstatistics(cfg, correct_scene_colour{:}, incorrect_scene_colour{:});
length(find(Fieldtripstats_incorrect.mask==1))

% shuffled retrieval vs baseline
cfg.correctm            = 'cluster'; % 'no', cluster, bonferroni, fdr, holm;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 1;
cfg.numrandomization    = 500;%1000;%'all';
cfg.clustertail         = cfg.tail;
% shuffled retrieval data
[Fieldtripstats_shuff_retri] = ft_freqstatistics(cfg, shuff_retrieval{:}, baseline_z{:});
length(find(Fieldtripstats_shuff_retri.mask==1))


% empirical vs cirular shift
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
length(find(Fieldtripstats_circ.mask==1))


% empirical z-value vs baseline
cfg.latency             = [timeaxis(1) timeaxis(end)]; % ripple time
cfg.method              = 'analytic'; % 'montecarlo' 'analytic';
cfg.correctm            = 'no'; % 'no', cluster, bonferroni, fdr, holm;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 1;
cfg.clustertail         = cfg.tail;
[Fieldtripstats_z] = ft_freqstatistics(cfg, z_all{:}, baseline_z{:});
any(Fieldtripstats_z.mask(:))
length(find(Fieldtripstats_z.mask==1))


length(find(Fieldtripstats_z.mask==1))


%% plot (significant) t vals


freq_idx = nearest(Fieldtripstats.freq,0):nearest(Fieldtripstats.freq,3);

emp_data    = squeeze(mean(correct_dec_fine_all(:,freq_idx,:)-incorrect_dec_fine_all(:,freq_idx,:),2));
shuff_data  = squeeze(mean(correct_dec_fine_all(:,freq_idx,:)-shuffled_dec_fine_all(:,freq_idx,:),2));
shuff_z_data  = squeeze(mean(correct_dec_fine_all_zvalue(:,freq_idx,:),2));
adj_data    = squeeze(mean(correct_dec_fine_all(:,freq_idx,:)-corr_circshift_all(:,freq_idx,:),2));
shuffle_retrieval_data    = squeeze(mean(shuffled_retrieval_all(:,freq_idx,:),2));


load mask_t_vals_fine_enc_ripple

figure('units','normalized','outerposition',[0 0 1 1]);


stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));


% shuffled retrieval data

gcf = subplot(2,3,1);
ax1 = gcf;
ax1.Position(2) = ax1.Position(2) + 0.03; 


imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    squeeze(mean(shuffled_retrieval_all)));
colormap(flipud(settings.colour_scheme_1))
caxis([0.4 .6])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title('random "ripples" retrieval trials');


plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(mask_t_vals_fine_enc_ripple)))


gcf             = subplot('Position',[.5 .85 .10 .10]);
ax1 = gcf;
ax1.Position(1) = ax1.Position(1) - 0.31; 

tvals = squeeze(Fieldtripstats_shuff_retri.stat);
imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    squeeze(tvals));
colormap(flipud(settings.colour_scheme_1))
% caxis([0.45 .65])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',14)
xlabel('Ripple time (sec)')
ylabel('Encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');


plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(Fieldtripstats_shuff_retri.mask)))




gcf = subplot(2,3,2);
ax1 = gcf;
ax1.Position(2) = ax1.Position(2) + 0.03; 


imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    squeeze(mean(shuffled_dec_fine_all)));
colormap(flipud(settings.colour_scheme_1))
caxis([0.4 .6])
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



gcf = subplot(2,3,3);
ax1 = gcf;
ax1.Position(2) = ax1.Position(2) + 0.03; 


imagesc(...
    correct{1, 1}.time(stats_time),...
    freqaxis(stats_freq),...
    squeeze(mean(corr_circshift_all)));
colormap(flipud(settings.colour_scheme_1))
caxis([0.4 .6])
vline(0)
hline(0)

axis xy
set(gca,'FontSize',20)
xlabel('Ripple time (sec)')
ylabel('Encoding time (sec)')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');
title('adjacent trials (circular shift)');


plot_contour(correct{1, 1}.time(stats_time),freqaxis(stats_freq),double(squeeze(mask_t_vals_fine_enc_ripple)))



dec_to_corr_real        = [];
dec_to_corr_shuff       = [];
dec_to_corr_circshift   = [];
dec_to_corr_shuff_retrieval = [];
dec_to_corr_shuff_z = [];
dec_to_corr_shuff_retrieval_z = [];
for isubject = 1:size(correct_dec_fine_all,1)
    tmp = squeeze((correct_dec_fine_all(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_real(isubject,:) = mean(mean(tmp));
    tmp = squeeze((shuffled_dec_fine_all(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_shuff(isubject,:) = mean(mean(tmp));
    tmp = squeeze((corr_circshift_all(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_circshift(isubject,:) = mean(mean(tmp));
    tmp = squeeze((shuffled_retrieval_all(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_shuff_retrieval(isubject,:) = mean(mean(tmp));

    % z-value shuffled
    tmp = squeeze((correct_dec_fine_all_zvalue(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_shuff_z(isubject,:) = mean(mean(tmp));
    % z-value shuffled retrieval
    tmp = squeeze((z_value_shuff_retrieval(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_shuff_retrieval_z(isubject,:) = mean(mean(tmp));

end


% 
% % shuffled retrieval trials
% 
% gcf = subplot(2,3,4);
% ax1 = gcf;
% ax1.Position(2) = ax1.Position(2) + 0.03;  % Shift downwards
% 
% 
% N = 12;
% color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
% color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector
% 
% scattercorrect = scatter(.9, dec_to_corr_shuff_retrieval_z, 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% hold on;
% 
% [f, xi] = ksdensity(dec_to_corr_shuff_retrieval_z); 
% fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color2)
% 
% 
% 
% 
% hold on
% set(gca,'TickDir','out')
% ylabel('z-values within cluster')
% 
% xticks([.8])
% xticklabels({'AM+ z-value'})
% xlabel('Condition')
% set(gca,'FontSize',20)
% hline(0)
% box off
% 
% range=axis;
% 
% 
% % shuffled ripples ripple aligned trials
% 
% gcf = subplot(2,3,5);
% ax1 = gcf;
% ax1.Position(2) = ax1.Position(2) + 0.03;  % Shift downwards
% 
% 
% N = 12;
% color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
% color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector
% 
% scattercorrect = scatter(.9, dec_to_corr_shuff_z, 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% hold on;
% 
% [f, xi] = ksdensity(dec_to_corr_shuff_z); 
% fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color2)
% 
% 
% 
% 
% hold on
% set(gca,'TickDir','out')
% ylabel('z-values within cluster')
% 
% xticks([.8])
% xticklabels({'AM+ z-value'})
% xlabel('Condition')
% set(gca,'FontSize',20)
% hline(0)
% box off
% 
% range=axis;
% 
% 
% % adjacent trials
% 
% gcf = subplot(2,3,6);
% ax1 = gcf;
% ax1.Position(2) = ax1.Position(2) + 0.03;  % Shift downwards
% 
% 
% 
% N = 12;
% color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
% color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector
% 
% scattercorrect = scatter(.9, dec_to_corr_real, 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% hold on;
% scatterincorrect = scatter(1.1, dec_to_corr_circshift, 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color1, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% 
% 
% [f, xi] = ksdensity(dec_to_corr_real); 
% fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color2)
% 
% % Second kernel density estimation and fill, mirrored
% [f, xi] = ksdensity(dec_to_corr_circshift); 
% fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color1) 
% 
% for i = 1:N
%     hold on;
%     
%     
%     
%     line([.9 1.1], [dec_to_corr_real, dec_to_corr_circshift], 'Color', [0.5 0.5 0.5]);
% end
% 
% 
% hold on
% set(gca,'TickDir','out')
% ylabel('Deoding acc. within cluster')
% 
% xticks([.8 1.2])
% xticklabels({'AM+', 'AM+ adjacent'})
% xlabel('Condition')
% set(gca,'FontSize',20)
% box off
% 
% range=axis;






stats = [];
[~,pvalue.dec.real_shuffl_retrieval, ~,d] = ttest(dec_to_corr_real,dec_to_corr_shuff_retrieval)
stats.dec.real_shuffl_retrieval = d.tstat;
[~,pvalue.dec.real_shuff, ~,d] = ttest(dec_to_corr_real,dec_to_corr_shuff)
stats.dec.real_shuff = d.tstat;
[~,pvalue.dec.real_circshift, ~,d] = ttest(dec_to_corr_real,dec_to_corr_circshift)
stats.dec.real_circshift = d.tstat;




% also calculate the z-value compared to 0 in cluster

dec_to_corr_real        = [];
shuffled_retrieval_z_mask = [];

for isubject = 1:size(correct_dec_fine_all,1)
    tmp = squeeze((z_on_group(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    dec_to_corr_real(isubject,:) = mean(mean(tmp));


    % and for shuffled "ripples" retrieval
    tmp = squeeze((z_value_shuff_retrieval(isubject,:,:))).*mask_t_vals_fine_enc_ripple;
    shuffled_retrieval_z_mask(isubject,:) = mean(mean(tmp));
    
    
end


[~,pvalue.dec.real_shuff_z, ~,d] = ttest(dec_to_corr_real);
stats.dec.real_shuff_z = d.tstat;

[~,pvalue.dec.real_shuff_retrieval_z, ~,d] = ttest(shuffled_retrieval_z_mask);
stats.dec.real_shuff_retrieval_z = d.tstat;

pvalue.dec
stats.dec


subplot(2,3,4)

pvals = [p_emp_vs_retrieval_correct; ...
         p_emp_vs_random_correct; ...
         pvalue.dec.real_circshift];

y = -log10(pvals);

hold on
bar(y, 'FaceColor', [0.7 0.7 0.7])
yline(-log10(0.05), '--k', 'p = 0.05');
yline(-log10(0.01), ':k', 'p = 0.01');

set(gca, 'XTick', 1:3, ...
         'XTickLabel', {'RetrCorrect','RandCorrect','CircCorrect'});
ylabel('-log_{10}(p)');
box off

for i = 1:numel(pvals)
    text(i, y(i)+0.1, sprintf('p = %.3f', pvals(i)), ...
         'HorizontalAlignment','center');
end

set(gca,'FontSize',20)

%% Dimensionality
clear
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

load perf_dimensionality

perf = perf_dimensionality;

load perf_dim_shuffle_time_of_ripples

perf_shuffle_retr = perf_dim_shuffle_time_of_ripples;

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
correct_incorrect{7}                = correct_incorrect{1};
correct_incorrect{8}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = mean(perf_shuf_corr{1,isubject}.correct.accuracy(:,xlimits));

    tmp = perf{isubject}.correct.accuracy(xlimits);
    m =  mean(perf_shuf_corr{isubject}.correct.accuracy(:,xlimits));
    s = std(perf_shuf_corr{isubject}.correct.accuracy(:,xlimits));

    correct_incorrect{3}.individual(isubject,1,:)   = (tmp-m)/s;
    correct_incorrect{4}.individual(isubject,1,:)   = perf{1,isubject}.incorrect.accuracy(xlimits);

    correct_incorrect{5}.individual(isubject,1,:)           = (perf_circshift_corr{1,isubject}.correct.accuracy(xlimits)+perf_circshift_corr_back{1,isubject}.correct.accuracy(xlimits))/2;
  
    correct_incorrect{6}.individual(isubject,1,:)   = zeros(size(perf{1,isubject}.incorrect.accuracy(xlimits)));

     tmp = perf{isubject}.correct.accuracy(xlimits);
    m =  mean(perf_shuffle_retr{isubject}.correct.accuracy);
    s = std(perf_shuffle_retr{isubject}.correct.accuracy);
    correct_incorrect{7}.individual(isubject,1,:)   =(tmp-m)/s;
    correct_incorrect{8}.individual(isubject,1,:)   = m;
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);
correct_incorrect{3}.avg            = squeeze(correct_incorrect{3}.individual);
correct_incorrect{4}.avg            = squeeze(correct_incorrect{4}.individual);
correct_incorrect{5}.avg            = squeeze(correct_incorrect{5}.individual);
correct_incorrect{6}.avg            = squeeze(correct_incorrect{6}.individual);
correct_incorrect{7}.avg            = squeeze(correct_incorrect{7}.individual);
correct_incorrect{8}.avg            = squeeze(correct_incorrect{8}.individual);




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

cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.tail                = 1;
cfg.clustertail         = cfg.tail;

[Fieldtripstats_shuffl_retri] = ft_timelockstatistics(cfg, correct_incorrect{1,[7,6]});
length(find(Fieldtripstats_shuffl_retri.mask))





%% Stelzer method
stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline_emp   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline_emp(stats_time(Fieldtripstats_corr_inc.mask==1)) = 1;



Nsub = numel(perf);

emp_sub = zeros(Nsub,1);

for isubject = 1:Nsub
    emp_map = perf{1,isubject}.correct.accuracy(xlimits);

    emp_sub(isubject) = mean(emp_map(sigline_emp==1));  % masked mean (no zeros)
end

emp_group_fine_grained_correct = mean(emp_sub);


B = 10000;                  % iterations for group-level null
K = 1000;                   % surrogates per subject (trials)

null_group = nan(B,1);

for b = 1:B
    sub_vals = nan(Nsub,1);

    for s = 1:Nsub
        k = randi(K);  % pick one surrogate index for this subject

        shuf_map = perf_shuf_corr{1,isubject}.correct.accuracy(k,xlimits);

         sub_vals(s) = mean(shuf_map(sigline_emp==1));  % masked mean
    end

    null_group(b) = mean(sub_vals);
end

null_group_fine_grained = null_group;

p_emp_vs_random_correct = (sum(abs(null_group_fine_grained - mean(null_group_fine_grained)) >= abs(emp_group_fine_grained_correct - mean(null_group_fine_grained))) + 1) / (B + 1)

z_emp_vs_random_correct = (emp_group_fine_grained_correct - mean(null_group_fine_grained)) / std(null_group_fine_grained)


% now against shuffled retrieval

B = 10000;                  % iterations for group-level null
K = 40;                   % surrogates per subject (trials)

null_group = nan(B,1);

for b = 1:B
    sub_vals = nan(Nsub,1);

    for s = 1:Nsub
        k = randi(K);  % pick one surrogate index for this subject

        shuf_map = perf_shuffle_retr{1,isubject}.correct.accuracy(k,:);

         sub_vals(s) = mean(shuf_map(sigline_emp==1));  % masked mean
    end

    null_group(b) = mean(sub_vals);
end



null_group_retrieval = null_group;

p_emp_vs_retrieval_correct = (sum(abs(null_group - mean(null_group)) >= abs(emp_group_fine_grained_correct - mean(null_group))) + 1) / (B + 1);

z_emp_vs_retrieval_correct = (emp_group_fine_grained_correct - mean(null_group_retrieval)) / std(null_group_retrieval);



%% plot significant values

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,3,1)

d = squeeze(correct_incorrect{1}.individual);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));

dim_to_correlate_corr = d;

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline_emp   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline_emp(stats_time(Fieldtripstats_corr_inc.mask==1)) = 2;

dim_to_correlate_corr = mean(dim_to_correlate_corr(:,sigline_emp==2),2);

boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(23,:), 'alpha');
hold on
plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

hold on

d = squeeze(correct_incorrect{8}.individual);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));


boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(7,:), 'alpha');
hold on
plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 1);

if ~isempty(sig_times_1)
    fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 2);

if ~isempty(sig_times_1)
    fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end


set(gca, 'FontSize', 20, 'FontName', 'Arial')
xlabel('Ripple time (sec)')
ylabel('Dimensionality')
title('RetrCorrect')
set(gca, 'TickDir', 'out')
ylim([2 3.7])
vline(0)


d = squeeze(correct_incorrect{7}.individual);

dim_to_correlate_corr_shuff_retr = d;

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline_emp   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline_emp(stats_time(Fieldtripstats_corr_inc.mask==1)) = 2;

dim_to_correlate_corr_shuff_retr = mean(dim_to_correlate_corr_shuff_retr(:,sigline_emp==2),2);




d = squeeze(correct_incorrect{7}.individual);


m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));


% boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(7,:), 'alpha');
% hold on
% plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_shuffl_retri.posclusterslabelmat == 1);


% if ~isempty(sig_times_1)
%     fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
%          [-.7 -.7 2 2], ...
%          'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% end
% 
% 
% set(gca, 'FontSize', 20, 'FontName', 'Arial')
% xlabel('Ripple time (sec)')
% ylabel('Z-values')
% title('Random "ripples" retrieval trials')
% set(gca, 'TickDir', 'out')
% ylim([-.7 2])
% vline(0)



d = squeeze(correct_incorrect{3}.individual);

dim_to_correlate_corr_shuffled = d;

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline_emp   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline_emp(stats_time(Fieldtripstats_corr_inc.mask==1)) = 2;

dim_to_correlate_corr_shuffled = mean(dim_to_correlate_corr_shuffled(:,sigline_emp==2),2);


subplot(2,3,2)


d = squeeze(correct_incorrect{1}.individual);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));

dim_to_correlate_corr = d;

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline_emp   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline_emp(stats_time(Fieldtripstats_corr_inc.mask==1)) = 2;

dim_to_correlate_corr = mean(dim_to_correlate_corr(:,sigline_emp==2),2);

boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(23,:), 'alpha');
hold on
plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

hold on

d = squeeze(correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));


boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(7,:), 'alpha');
hold on
plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 1);

if ~isempty(sig_times_1)
    fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end

sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_corr_inc.posclusterslabelmat == 2);

if ~isempty(sig_times_1)
    fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
         [2 2 3.7 3.7], ...
         'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
end


set(gca, 'FontSize', 20, 'FontName', 'Arial')
xlabel('Ripple time (sec)')
ylabel('Dimensionality')
title('RandCorrect')
set(gca, 'TickDir', 'out')
ylim([2 3.7])
vline(0)




d = squeeze(correct_incorrect{3}.individual);

dim_to_correlate_shuff_corr = mean(d(:,sigline_emp == 2), 2);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));


% boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(7,:), 'alpha');
% hold on
% plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

sig_times_1 = correct_incorrect{1}.time(Fieldtripstats_z.posclusterslabelmat == 1);

% 
% if ~isempty(sig_times_1)
%     fill([sig_times_1(1) sig_times_1(end) sig_times_1(end) sig_times_1(1)], ...
%          [-.7 -.7 2 2], ...
%          'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% end


% 
% set(gca, 'FontSize', 20, 'FontName', 'Arial')
% xlabel('Ripple time (sec)')
% ylabel('Z-values')
% title('Shuffled correct trials')
% set(gca, 'TickDir', 'out')
% ylim([-.7 2])
% vline(0)



subplot(2,3,3)

d = squeeze(correct_incorrect{1}.individual);

m = nanmean(d);
s = nanstd(d) ./ sqrt(size(d, 1));

dim_to_correlate_corr = d;

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline_emp   = nan(1,numel(correct_incorrect{1}.time));
%         sigline(stats_time(Fieldtripstats.mask==1)) = m(stats_time(Fieldtripstats.mask==1));
sigline_emp(stats_time(Fieldtripstats_corr_inc.mask==1)) = 2;

dim_to_correlate_corr = mean(dim_to_correlate_corr(:,sigline_emp==2),2);

boundedline(correct_incorrect{1, 1}.time, m, s, 'cmap', settings.colour_scheme_1(23,:), 'alpha');
hold on
plot(correct_incorrect{1, 1}.time, m, 'k', 'linewidth', 2);

hold on

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
xlabel('Ripple time (sec)')
ylabel('Dimensionality')
title('CircCorrect')
set(gca, 'TickDir', 'out')
ylim([2 3.7])
vline(0)


subplot(2,3,4)


% 
% N = 12;
% color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
% color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector
% 
% scattercorrect = scatter(.9, dim_to_correlate_corr_shuff_retr, 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% hold on
% [f, xi] = ksdensity(dim_to_correlate_corr_shuff_retr); 
% fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color2)
% 
% 
% hold on
% set(gca,'TickDir','out')
% ylabel('Z-values')
% 
% xticks([.8])
% xticklabels({'AM+ z-values'})
% xlabel('Trial type')
% set(gca,'FontSize',20)
% hline(0)
% box off
% 
% range=axis;
% 
% 
% 
% 
% 
% 
% 
% subplot(2,3,5)
% 
% 
% N = 12;
% color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
% color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector
% 
% scattercorrect = scatter(.9, dim_to_correlate_shuff_corr, 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% hold on;
% 
% [f, xi] = ksdensity(dim_to_correlate_shuff_corr); 
% fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color2)
% 
% 
% hold on
% set(gca,'TickDir','out')
% ylabel('Z-values')
% 
% xticks([.8])
% xticklabels({'AM+ z-values'})
% xlabel('Trial type')
% set(gca,'FontSize',20)
% hline(0)
% box off
% 
% range=axis;
% 
% 
% 
% subplot(2,3,6)
% 
% 
% 
% 
% data_to_plot = [dim_to_correlate_corr,dim_to_correlate_ajdacent_corr];
% 
% 
% N = 12;
% color1 = settings.colour_scheme_1(7, :); % Should be a 1x3 vector
% color2 = settings.colour_scheme_1(23, :); % Should also be a 1x3 vector
% 
% scattercorrect = scatter(.9, data_to_plot(:,1), 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color2, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% hold on;
% scattercorrect = scatter(1.1, data_to_plot(:,2), 75, 'MarkerEdgeColor', 'black',...
%     'MarkerFaceColor', color1, 'MarkerFaceAlpha',.3,...
%     'MarkerEdgeAlpha',.8);
% hold on;
% 
% [f, xi] = ksdensity(data_to_plot(:,1)); 
% fill(0.8 - (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color2)
% hold on;
% 
% [f, xi] = ksdensity(data_to_plot(:,2)); 
% fill(1.2 + (f / max(f) * 0.8), xi, 'k', 'FaceAlpha', .3,...
%     'EdgeColor', 'none', 'FaceColor', color1)
% 
% for i = 1:N
%     hold on;
% 
%     line([.9 1.1], [dim_to_correlate_corr, dim_to_correlate_ajdacent_corr], 'Color', [0.5 0.5 0.5]);
% end
% 
% 
% hold on
% set(gca,'TickDir','out')
% ylabel('Dimensionality')
% 
% xticks([.8 1.2])
% xticklabels({'AM+', 'AM+ adjacent'})
% xlabel('Trial type')
% set(gca,'FontSize',20)
% box off
% 
% range=axis;


[~,pvalue.dim.shuff_corr, ~,d] = ttest(dim_to_correlate_corr,dim_to_correlate_shuff_corr)
stats.dim.shuff_corr = d.tstat
[~,pvalue.dim.ajd_corr, ~,d] = ttest(dim_to_correlate_corr,dim_to_correlate_ajdacent_corr)
stats.dim.ajd_corr = d.tstat
[~,pvalue.dim.shuff_corr_retr_z, ~,d] = ttest(dim_to_correlate_corr_shuff_retr)
stats.dim.shuff_corr_retr_z = d.tstat
[~,pvalue.dim.shuff_z, ~,d] = ttest(dim_to_correlate_corr_shuffled)
stats.dim.shuff_z = d.tstat


pvalue.dim
stats.dim


pvals = [p_emp_vs_retrieval_correct; ...
         p_emp_vs_random_correct; ...
         pvalue.dim.ajd_corr];

y = -log10(pvals);

hold on
bar(y, 'FaceColor', [0.7 0.7 0.7])
yline(-log10(0.05), '--k', 'p = 0.05');
yline(-log10(0.01), ':k', 'p = 0.01');

set(gca, 'XTick', 1:3, ...
         'XTickLabel', {'RetrCorrect','RandCorrect','CircCorrect'});
ylabel('-log_{10}(p)');
box off

for i = 1:numel(pvals)
    text(i, y(i)+0.1, sprintf('p = %.3f', pvals(i)), ...
         'HorizontalAlignment','center');
end

set(gca,'FontSize',20)