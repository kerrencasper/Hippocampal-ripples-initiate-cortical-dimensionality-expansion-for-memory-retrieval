% Figure 3 - Dimensionality transformation are locked to ripple events
%                  Casper Kerren      [kerren@cbs.mpg.de]


clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

% [~,ftpath]=ft_version;

%% path settings
settings                    = [];

settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

cwd = pwd;

addpath(genpath([pwd,'/help_functions']))
addpath([pwd,'/dPCA'])


settings.colour_scheme_1 = brewermap(30,'RdBu');
settings.colour_scheme_1 = settings.colour_scheme_1;

%% LOAD and organise

load data_dPCA.mat


dataToClassifyTest = data_dPCA;


fs          = 100;
timeEvents  = -1:1/fs:1;
TOI         = nearest(timeEvents,-1):nearest(timeEvents,1);
time        = timeEvents(TOI);
timeEvents  = 0;

dataToClassifyTest = dataToClassifyTest(:,:,:,TOI);

N = size(dataToClassifyTest, 3);  % Number of channels (neurons)
S = 4;  % Number of stimuli conditions (blue, red, indoor, outdoor)
D = 2;  % Number of memory (correct, incorrect)
T = size(dataToClassifyTest, 4);  % Number of time points

trials_per_condition    = 46;
trialNum                = repmat(trials_per_condition, [N, S, D]);  % Same number of trials per channel, condition, memory

maxTrialNum = 46;
firingRates = nan(N, S, D, T, maxTrialNum);

for d = 1:D
    for s = 1:S

        data_subset = dataToClassifyTest(d, ((s-1)*trials_per_condition + 1):(s*trials_per_condition), :, :);  % [trials, channels, time]

        % Permute to get [channels, trials, time] -> [N, trials_per_condition, T]
        data_subset = permute(data_subset, [3, 1, 4, 2]);  % Now size is [N, trials, T]


        firingRates(:, s, d, :, 1:trials_per_condition) = data_subset;  % Filling in [N, 1, 1, T, maxTrialNum]
    end
end

firingRatesAverage = nanmean(firingRates, 5);

firingRates(isnan(firingRates))=0;

for icond = 1:D
    for ichan = 1:N
        tmp = squeeze(firingRates(ichan,:,icond,1,:));
        firstZeroIndex = zeros(size(tmp, 1), 1);

        for i = 1:size(tmp, 1)

            zeroIndices = find(tmp(i, :) == 0);

            if ~isempty(zeroIndices)

                firstZeroIndex(i) = zeroIndices(1);
            else

                firstZeroIndex(i) = trials_per_condition;
            end
        end

        trialNum(ichan,:,icond) = firstZeroIndex-1;
        if trialNum(ichan,:,icond) == 45
            idx = find(trialNum(ichan,:,icond) == 45);
            trialNum(ichan,idx,icond) = 46;
        end
    end
end


Xtrial = permute(firingRates, [1, 2, 4, 5, 3]);

combinedParams = {{1,[1 3]}, {2,[2 3]}, {3}, {[1,2],[1,2,3]}};
margNames = {'Stimulus', 'Memory', 'Ind.', 'S x D'};
numMargNames = length(margNames);
selectedColors = settings.colour_scheme_1(5:25, :);
margColours = selectedColors(round(linspace(1, size(selectedColors, 1), numMargNames)), :);
margNameColourMap = containers.Map(margNames, num2cell(margColours, 2));


% Check consistency between trialNum and firingRates
for n = 1:size(firingRates, 1)
    for s = 1:size(firingRates, 2)
        assert(isempty(find(isnan(firingRates(n, s, :, 1:trialNum(n, s))), 1)), 'Something is wrong!');
    end
end

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X, 2));


%% parameter settings

combinedParams = {{1,[1 3]}, {2,[2 3]}, {3}, {[1,2],[1,2,3]}};
margNames = {'Target Assoc.', 'Memory', 'Ind.', 'T x M'};
numMargNames = length(margNames);
selectedColors = settings.colour_scheme_1(5:25, :);
margColours = selectedColors(round(linspace(1, size(selectedColors, 1), numMargNames)), :);
margNameColourMap = containers.Map(margNames, num2cell(margColours, 2));


%% load data

load('accuracy.mat')
load('accuracyShuffle.mat')
load('Cnoise.mat')
load('combinedParams.mat')
load('explVar.mat')
load('firingRates.mat')
load('firingRatesAverage.mat')
load('optimalLambda.mat')
load('options.mat')
load('trialNum.mat')
load('V.mat')
load('W.mat')
load('whichMarg.mat')
load('Z.mat')

whichMarg = options.whichMarg;


options.explainedVar            = explVar;
options.marginalizationNames    = margNames;
options.marginalizationColours  = margColours;
options.whichMarg               = whichMarg;
options.time                    = time;
options.timeEvents              = timeEvents;
options.timeMarginalization     = 3;
options.legendSubplot           = 16;
options_plotting_dpca           = options;


%% show decoding significance - one-sided test against shuffle
componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg,'minChunk',10);

which_comp = whichMarg(1,1:4);

dpca_to_plot = zeros(size(which_comp,2), size(time,2));



alpha = 0.05;
z_threshold = norminv(1 - alpha/2); 

marg        = which_comp(1);
num         = 1;

m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

a = (squeeze(accuracy(marg,num,:))-m)./s';
dpca_to_plot(1,:) = a>z_threshold;

critical_alpha = .05; % alpha for after permutation
alpha = 0.05; % cluster alpha
z_threshold = norminv(1 - alpha/2);

marg = which_comp(1);
num = 1;

% Compute mean and std from shuffled data
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

% Compute z-scored accuracy
a = (squeeze(accuracy(marg,num,:)) - m) ./ s';
observed_data = a;

% Threshold the observed data
dpca_to_plot(1,:) = observed_data > z_threshold;

% Identify clusters in the observed data
observed_clusters = bwlabeln(dpca_to_plot(1,:) > 0);
cluster_labels = unique(observed_clusters);
cluster_labels = cluster_labels(cluster_labels > 0); % exclude zero

% Calculate cluster-level statistics (sum of z-values)
observed_cluster_stats = arrayfun(@(c) sum(observed_data(observed_clusters == c)), cluster_labels);

% Permutation test
num_permutations = 1000;
permuted_max_cluster_stats = zeros(1, num_permutations);

for perm = 1:num_permutations
    % Permute data
    permuted = observed_data(randperm(length(observed_data)));
    perm_clusters = bwlabel(permuted > z_threshold);
    perm_cluster_labels = unique(perm_clusters);
    perm_cluster_labels = perm_cluster_labels(perm_cluster_labels > 0);
    if ~isempty(perm_cluster_labels)
        perm_cluster_stats = arrayfun(@(c) sum(permuted(perm_clusters == c)), perm_cluster_labels);
        permuted_max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        permuted_max_cluster_stats(perm) = 0;
    end
end


% Calculate cluster-level significance
significant_clusters = false(size(cluster_labels));
cluster_threshold = prctile(permuted_max_cluster_stats, 100 * (1 - critical_alpha));

for i = 1:length(cluster_labels)
    significant_clusters(i) = observed_cluster_stats(i) > cluster_threshold;
end
% 

significant_mask = zeros(size(observed_data));
% 
% Fill in the mask for significant clusters
for i = 1:length(cluster_labels)
    cluster_indices = find(observed_clusters == cluster_labels(i));
    if significant_clusters(i)
        significant_mask(cluster_indices) = 1;
    end
end


dpca_to_plot(1,:) = significant_mask;



marg        = which_comp(2);
num         = 1;

m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

a = (squeeze(accuracy(marg,num,:))-m)./s';
dpca_to_plot(2,:) = a>z_threshold;


% Compute mean and std from shuffled data
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

% Compute z-scored accuracy
a = (squeeze(accuracy(marg,num,:)) - m) ./ s';
observed_data = a;

% Threshold the observed data
dpca_to_plot(2,:) = observed_data > z_threshold;

% Identify clusters in the observed data
observed_clusters = bwlabel(dpca_to_plot(2,:) > 0);
cluster_labels = unique(observed_clusters);
cluster_labels = cluster_labels(cluster_labels > 0); % exclude zero

% Calculate cluster-level statistics (sum of z-values)
observed_cluster_stats = arrayfun(@(c) sum(observed_data(observed_clusters == c)), cluster_labels);

% Permutation test
num_permutations = 1000;
permuted_max_cluster_stats = zeros(1, num_permutations);

for perm = 1:num_permutations
    % Permute data
    permuted = observed_data(randperm(length(observed_data)));
    perm_clusters = bwlabel(permuted > z_threshold);
    perm_cluster_labels = unique(perm_clusters);
    perm_cluster_labels = perm_cluster_labels(perm_cluster_labels > 0);
    if ~isempty(perm_cluster_labels)
        perm_cluster_stats = arrayfun(@(c) sum(permuted(perm_clusters == c)), perm_cluster_labels);
        permuted_max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        permuted_max_cluster_stats(perm) = 0;
    end
end

% Calculate cluster-level significance
significant_clusters = false(size(cluster_labels));
cluster_threshold = prctile(permuted_max_cluster_stats, 100 * (1 - critical_alpha));

for i = 1:length(cluster_labels)
    significant_clusters(i) = observed_cluster_stats(i) > cluster_threshold;
end
% 

significant_mask = zeros(size(observed_data));

% Fill in the mask for significant clusters
for i = 1:length(cluster_labels)
    cluster_indices = find(observed_clusters == cluster_labels(i));
    if significant_clusters(i)
        significant_mask(cluster_indices) = 1;
    end
end


dpca_to_plot(2,:) = significant_mask;


marg        = which_comp(3);
num         = 2;

m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

a = (squeeze(accuracy(marg,num,:))-m)./s';
dpca_to_plot(3,:) = a>z_threshold;

% Compute mean and std from shuffled data
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

% Compute z-scored accuracy
a = (squeeze(accuracy(marg,num,:)) - m) ./ s';
observed_data = a;

% Threshold the observed data
dpca_to_plot(3,:) = observed_data > z_threshold;

% Identify clusters in the observed data
observed_clusters = bwlabel(dpca_to_plot(3,:) > 0);
cluster_labels = unique(observed_clusters);
cluster_labels = cluster_labels(cluster_labels > 0); % exclude zero

% Calculate cluster-level statistics (sum of z-values)
observed_cluster_stats = arrayfun(@(c) sum(observed_data(observed_clusters == c)), cluster_labels);

% Permutation test
num_permutations = 1000;
permuted_max_cluster_stats = zeros(1, num_permutations);

for perm = 1:num_permutations
    % Permute data
    permuted = observed_data(randperm(length(observed_data)));
    perm_clusters = bwlabel(permuted > z_threshold);
    perm_cluster_labels = unique(perm_clusters);
    perm_cluster_labels = perm_cluster_labels(perm_cluster_labels > 0);
    if ~isempty(perm_cluster_labels)
        perm_cluster_stats = arrayfun(@(c) sum(permuted(perm_clusters == c)), perm_cluster_labels);
        permuted_max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        permuted_max_cluster_stats(perm) = 0;
    end
end

% Calculate cluster-level significance
significant_clusters = false(size(cluster_labels));
cluster_threshold = prctile(permuted_max_cluster_stats, 100 * (1 - critical_alpha));

for i = 1:length(cluster_labels)
    significant_clusters(i) = observed_cluster_stats(i) > cluster_threshold;
end
% 

significant_mask = zeros(size(observed_data));

% Fill in the mask for significant clusters
for i = 1:length(cluster_labels)
    cluster_indices = find(observed_clusters == cluster_labels(i));
    if significant_clusters(i)
        significant_mask(cluster_indices) = 1;
    end
end


dpca_to_plot(3,:) = significant_mask;



marg        = which_comp(4);
num         = 1;
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

a = (squeeze(accuracy(marg,num,:))-m)./s';
dpca_to_plot(4,:) = a>z_threshold;

% Compute mean and std from shuffled data
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

% Compute z-scored accuracy
a = (squeeze(accuracy(marg,num,:)) - m) ./ s';
observed_data = a;

% Threshold the observed data
dpca_to_plot(4,:) = observed_data > z_threshold;

% Identify clusters in the observed data
observed_clusters = bwlabel(dpca_to_plot(4,:) > 0);
cluster_labels = unique(observed_clusters);
cluster_labels = cluster_labels(cluster_labels > 0); % exclude zero

% Calculate cluster-level statistics (sum of z-values)
observed_cluster_stats = arrayfun(@(c) sum(observed_data(observed_clusters == c)), cluster_labels);

% Permutation test
num_permutations = 1000;
permuted_max_cluster_stats = zeros(1, num_permutations);

for perm = 1:num_permutations
    % Permute data
    permuted = observed_data(randperm(length(observed_data)));
    perm_clusters = bwlabel(permuted > z_threshold);
    perm_cluster_labels = unique(perm_clusters);
    perm_cluster_labels = perm_cluster_labels(perm_cluster_labels > 0);
    if ~isempty(perm_cluster_labels)
        perm_cluster_stats = arrayfun(@(c) sum(permuted(perm_clusters == c)), perm_cluster_labels);
        permuted_max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        permuted_max_cluster_stats(perm) = 0;
    end
end

% Calculate cluster-level significance
significant_clusters = false(size(cluster_labels));
cluster_threshold = prctile(permuted_max_cluster_stats, 100 * (1 - critical_alpha));

for i = 1:length(cluster_labels)
    significant_clusters(i) = observed_cluster_stats(i) > cluster_threshold;
end
% 

significant_mask = zeros(size(observed_data));

% Fill in the mask for significant clusters
for i = 1:length(cluster_labels)
    cluster_indices = find(observed_clusters == cluster_labels(i));
    if significant_clusters(i)
        significant_mask(cluster_indices) = 1;
    end
end


dpca_to_plot(4,:) = significant_mask;




marg        = which_comp(2);
num         = 3;
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

a = (squeeze(accuracy(marg,num,:))-m)./s';
dpca_to_plot(5,:) = a>z_threshold;

% Compute mean and std from shuffled data
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

% Compute z-scored accuracy
a = (squeeze(accuracy(marg,num,:)) - m) ./ s';
observed_data = a;

% Threshold the observed data
dpca_to_plot(5,:) = observed_data > z_threshold;

% Identify clusters in the observed data
observed_clusters = bwlabel(dpca_to_plot(5,:) > 0);
cluster_labels = unique(observed_clusters);
cluster_labels = cluster_labels(cluster_labels > 0); % exclude zero

% Calculate cluster-level statistics (sum of z-values)
observed_cluster_stats = arrayfun(@(c) sum(observed_data(observed_clusters == c)), cluster_labels);

% Permutation test
num_permutations = 1000;
permuted_max_cluster_stats = zeros(1, num_permutations);

for perm = 1:num_permutations
    % Permute data
    permuted = observed_data(randperm(length(observed_data)));
    perm_clusters = bwlabel(permuted > z_threshold);
    perm_cluster_labels = unique(perm_clusters);
    perm_cluster_labels = perm_cluster_labels(perm_cluster_labels > 0);
    if ~isempty(perm_cluster_labels)
        perm_cluster_stats = arrayfun(@(c) sum(permuted(perm_clusters == c)), perm_cluster_labels);
        permuted_max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        permuted_max_cluster_stats(perm) = 0;
    end
end

% Calculate cluster-level significance
significant_clusters = false(size(cluster_labels));
cluster_threshold = prctile(permuted_max_cluster_stats, 100 * (1 - critical_alpha));

for i = 1:length(cluster_labels)
    significant_clusters(i) = observed_cluster_stats(i) > cluster_threshold;
end
% 

significant_mask = zeros(size(observed_data));

% Fill in the mask for significant clusters
for i = 1:length(cluster_labels)
    cluster_indices = find(observed_clusters == cluster_labels(i));
    if significant_clusters(i)
        significant_mask(cluster_indices) = 1;
    end
end


dpca_to_plot(5,:) = significant_mask;


marg        = which_comp(1); % S X D
num         = 2;
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

a = (squeeze(accuracy(marg,num,:))-m)./s';
dpca_to_plot(6,:) = a>z_threshold;



% Compute mean and std from shuffled data
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

% Compute z-scored accuracy
a = (squeeze(accuracy(marg,num,:)) - m) ./ s';
observed_data = a;

% Threshold the observed data
dpca_to_plot(6,:) = observed_data > z_threshold;

% Identify clusters in the observed data
observed_clusters = bwlabel(dpca_to_plot(6,:) > 0);
cluster_labels = unique(observed_clusters);
cluster_labels = cluster_labels(cluster_labels > 0); % exclude zero

% Calculate cluster-level statistics (sum of z-values)
observed_cluster_stats = arrayfun(@(c) sum(observed_data(observed_clusters == c)), cluster_labels);

% Permutation test
num_permutations = 1000;
permuted_max_cluster_stats = zeros(1, num_permutations);

for perm = 1:num_permutations
    % Permute data
    permuted = observed_data(randperm(length(observed_data)));
    perm_clusters = bwlabel(permuted > z_threshold);
    perm_cluster_labels = unique(perm_clusters);
    perm_cluster_labels = perm_cluster_labels(perm_cluster_labels > 0);
    if ~isempty(perm_cluster_labels)
        perm_cluster_stats = arrayfun(@(c) sum(permuted(perm_clusters == c)), perm_cluster_labels);
        permuted_max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        permuted_max_cluster_stats(perm) = 0;
    end
end

% Calculate cluster-level significance
significant_clusters = false(size(cluster_labels));
cluster_threshold = prctile(permuted_max_cluster_stats, 100 * (1 - critical_alpha));

for i = 1:length(cluster_labels)
    significant_clusters(i) = observed_cluster_stats(i) > cluster_threshold;
end
% 

significant_mask = zeros(size(observed_data));

% Fill in the mask for significant clusters
for i = 1:length(cluster_labels)
    cluster_indices = find(observed_clusters == cluster_labels(i));
    if significant_clusters(i)
        significant_mask(cluster_indices) = 1;
    end
end


dpca_to_plot(6,:) = significant_mask;






marg        = which_comp(4); % Decision
num         = 2;
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

a = (squeeze(accuracy(marg,num,:))-m)./s';
dpca_to_plot(13,:) = a>z_threshold;




% Compute mean and std from shuffled data
m = mean(squeeze(accuracyShuffle(marg,:,1:end)),2);
s = std(squeeze(accuracyShuffle(marg,:,1:end))');

% Compute z-scored accuracy
a = (squeeze(accuracy(marg,num,:)) - m) ./ s';
observed_data = a;

% Threshold the observed data
dpca_to_plot(13,:) = observed_data > z_threshold;

% Identify clusters in the observed data
observed_clusters = bwlabel(dpca_to_plot(13,:) > 0);
cluster_labels = unique(observed_clusters);
cluster_labels = cluster_labels(cluster_labels > 0); % exclude zero

% Calculate cluster-level statistics (sum of z-values)
observed_cluster_stats = arrayfun(@(c) sum(observed_data(observed_clusters == c)), cluster_labels);

% Permutation test
num_permutations = 1000;
permuted_max_cluster_stats = zeros(1, num_permutations);

for perm = 1:num_permutations
    % Permute data
    permuted = observed_data(randperm(length(observed_data)));
    perm_clusters = bwlabel(permuted > z_threshold);
    perm_cluster_labels = unique(perm_clusters);
    perm_cluster_labels = perm_cluster_labels(perm_cluster_labels > 0);
    if ~isempty(perm_cluster_labels)
        perm_cluster_stats = arrayfun(@(c) sum(permuted(perm_clusters == c)), perm_cluster_labels);
        permuted_max_cluster_stats(perm) = max(perm_cluster_stats);
    else
        permuted_max_cluster_stats(perm) = 0;
    end
end

% Calculate cluster-level significance
significant_clusters = false(size(cluster_labels));
cluster_threshold = prctile(permuted_max_cluster_stats, 100 * (1 - critical_alpha));

for i = 1:length(cluster_labels)
    significant_clusters(i) = observed_cluster_stats(i) > cluster_threshold;
end
% 

significant_mask = zeros(size(observed_data));

% Fill in the mask for significant clusters
for i = 1:length(cluster_labels)
    cluster_indices = find(observed_clusters == cluster_labels(i));
    if significant_clusters(i)
        significant_mask(cluster_indices) = 1;
    end
end


dpca_to_plot(13,:) = significant_mask;




componentsSignif_95CI       = componentsSignif;
componentsSignif_95CI(1,:)  = dpca_to_plot(1,:);
componentsSignif_95CI(2,:)  = dpca_to_plot(2,:);
componentsSignif_95CI(3,:)  = dpca_to_plot(3,:);
componentsSignif_95CI(4,:)  = dpca_to_plot(4,:);
componentsSignif_95CI(5,:)  = dpca_to_plot(5,:);
componentsSignif_95CI(6,:)  = dpca_to_plot(6,:);
componentsSignif_95CI(13,:)  = dpca_to_plot(13,:);


dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16,                ...
    'componentsSignif', componentsSignif_95CI, ...
    'numCompToShow', 50);

options.componentsSignif        = componentsSignif_95CI;

%%
hFig = figure('units','normalized','outerposition',[0 0 1 1]);


load data_dPCA.mat


dataToClassifyTest = data_dPCA;


fs          = 100;
timeEvents  = -1:1/fs:1;
TOI         = nearest(timeEvents,-1):nearest(timeEvents,1);
time        = timeEvents(TOI);
timeEvents  = 0;

dataToClassifyTest = dataToClassifyTest(:,:,:,TOI);

N = size(dataToClassifyTest, 3);  % Number of channels (neurons)
S = 4;  % Number of stimuli conditions (blue, red, indoor, outdoor)
D = 2;  % Number of decisions (correct, incorrect)
T = size(dataToClassifyTest, 4);  % Number of time points

trials_per_condition    = 46;
trialNum                = repmat(trials_per_condition, [N, S, D]);  % Same number of trials per channel, condition, decision

maxTrialNum = 46;
firingRates = nan(N, S, D, T, maxTrialNum);

for d = 1:D
    for s = 1:S

        data_subset = dataToClassifyTest(d, ((s-1)*trials_per_condition + 1):(s*trials_per_condition), :, :);  % [trials, channels, time]

        % Permute to get [channels, trials, time] -> [N, trials_per_condition, T]
        data_subset = permute(data_subset, [3, 1, 4, 2]);  % Now size is [N, trials, T]


        firingRates(:, s, d, :, 1:trials_per_condition) = data_subset;  % Filling in [N, 1, 1, T, maxTrialNum]
    end
end

firingRatesAverage = nanmean(firingRates, 5);

firingRates(isnan(firingRates))=0;

for icond = 1:D
    for ichan = 1:N
        tmp = squeeze(firingRates(ichan,:,icond,1,:));
        firstZeroIndex = zeros(size(tmp, 1), 1);

        for i = 1:size(tmp, 1)

            zeroIndices = find(tmp(i, :) == 0);

            if ~isempty(zeroIndices)

                firstZeroIndex(i) = zeroIndices(1);
            else

                firstZeroIndex(i) = trials_per_condition;
            end
        end

        trialNum(ichan,:,icond) = firstZeroIndex-1;
        if trialNum(ichan,:,icond) == 45
            idx = find(trialNum(ichan,:,icond) == 45);
            trialNum(ichan,idx,icond) = 46;
        end
    end
end


Xtrial = permute(firingRates, [1, 2, 4, 5, 3]);

combinedParams = {{1,[1 3]}, {2,[2 3]}, {3}, {[1,2],[1,2,3]}};
margNames = {'Target Assoc.', 'Block', 'Ind.', 'T x B'};
numMargNames = length(margNames);
selectedColors = settings.colour_scheme_1(5:25, :);
margColours = selectedColors(round(linspace(1, size(selectedColors, 1), numMargNames)), :);
margNameColourMap = containers.Map(margNames, num2cell(margColours, 2));


% Check consistency between trialNum and firingRates
for n = 1:size(firingRates, 1)
    for s = 1:size(firingRates, 2)
        assert(isempty(find(isnan(firingRates(n, s, :, 1:trialNum(n, s))), 1)), 'Something is wrong!');
    end
end

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X, 2));


%% parameter settings

combinedParams = {{1,[1 3]}, {2,[2 3]}, {3}, {[1,2],[1,2,3]}};
margNames = {'Target Assoc', 'Memory', 'Ind.', 'T x M'};
numMargNames = length(margNames);
selectedColors = settings.colour_scheme_1(5:25, :);
margColours = selectedColors(round(linspace(1, size(selectedColors, 1), numMargNames)), :);
margNameColourMap = containers.Map(margNames, num2cell(margColours, 2));


%% load data


load('accuracy.mat')
load('accuracyShuffle.mat')
load('Cnoise.mat')
load('combinedParams.mat')
load('explVar.mat')
load('firingRates.mat')
load('firingRatesAverage.mat')
load('optimalLambda.mat')
load('options.mat')
load('trialNum.mat')
load('V.mat')
load('W.mat')
load('whichMarg.mat')
load('Z.mat')

whichMarg = options.whichMarg;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% plot components %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X_1     = firingRatesAverage(:,:)';
Xcen    = bsxfun(@minus, X_1, mean(X_1));
Z       = Xcen * W;
dataDim = size(firingRatesAverage);
Zfull   = reshape(Z(:,15)', [1 dataDim(2:end)]);


options = struct('time',           [], ...
    'whichMarg',      [], ...
    'timeEvents',     [], ...
    'ylims',          [], ...
    'componentsSignif', [], ...
    'timeMarginalization', [], ...
    'legendSubplot',  [], ...
    'marginalizationNames', [], ...
    'marginalizationColours', [], ...
    'explainedVar',   [], ...
    'numCompToShow',  15, ...
    'X_extra',        [], ...
    'showNonsignificantComponents', false);

options.explainedVar            = explVar;
options.marginalizationNames    = margNames;
options.marginalizationColours  = margColours;
options.whichMarg               = whichMarg;
options.time                    = time;
options.timeEvents              = timeEvents;
options.timeMarginalization     = 3;
options.legendSubplot           = 16;
options.componentsSignif        = componentsSignif_95CI;
options_plotting_dpca           = options;


toDisplayMargNames = 0;

% if there are 4 or less marginalizations, split them into rows
if ~isempty(options.whichMarg) && ...
        length(unique(options.whichMarg)) <= 4 && length(unique(options.whichMarg)) > 1

    % time marginalization, if specified, goes on top
    if ~isempty(options.timeMarginalization)
        margRowSeq = [options.timeMarginalization setdiff(1:max(options.whichMarg), options.timeMarginalization)];
    else
        margRowSeq = 1:max(options.whichMarg);
    end

    componentsToPlot = [];
    subplots = [];
    for i=1:length(margRowSeq)
        if ~isempty(options.componentsSignif) && margRowSeq(i) ~= options.timeMarginalization
            % selecting only significant components
            minL = min(length(options.whichMarg), size(options.componentsSignif,1));
            moreComponents = find(options.whichMarg(1:minL) == margRowSeq(i) & ...
                sum(options.componentsSignif(1:minL,:), 2)'~=0, 3);
            if options.showNonsignificantComponents && (length(moreComponents) < 3)
                % Optionally add non-significant components to fill subplots
                moreComponents = [moreComponents setdiff(find(options.whichMarg == margRowSeq(i), 3), moreComponents)];
            end
        else
            moreComponents = find(options.whichMarg == margRowSeq(i), 3);
        end
        componentsToPlot = [componentsToPlot moreComponents];
        subplots = [subplots (i-1)*4+2:(i-1)*4+2 + length(moreComponents) - 1];
    end
else
    % if there are more than 4 marginalizatons

    if isempty(options.whichMarg)
        % if there is no info about marginaliations
        componentsToPlot = 1:12;
    else
        % if there is info about marginaliations, select first 3 in each
        uni = unique(options.whichMarg);
        componentsToPlot = [];
        for u = 1:length(uni)
            componentsToPlot = [componentsToPlot find(options.whichMarg==uni(u), 2)];
        end
        componentsToPlot = sort(componentsToPlot);
        if length(componentsToPlot) > 12
            componentsToPlot = componentsToPlot(1:12);
        end

        toDisplayMargNames = 1;
    end
    subplots = [2 3 4 6 7 8 10 11 12 14 15 16];

    if numCompToShow < 12
        componentsToPlot = componentsToPlot(1:numCompToShow);
        subplots = subplots(1:numCompToShow);
    end
end

Zfull = reshape(Z(:,componentsToPlot)', [length(componentsToPlot) dataDim(2:end)]);

% exclude component 8, 9 and 18 which is number 1, 2, 3 (don't need to show this one)

% componentsToPlot(:,5)   = [];
componentsToPlot(:,3)   = [];
componentsToPlot(:,2)   = [];
componentsToPlot(:,1)   = [];
subplots                = [2 3 4];
% Zfull(5,:,:,:)          = [];
Zfull(3,:,:,:)          = [];
Zfull(2,:,:,:)          = [];
Zfull(1,:,:,:)          = [];


% y-axis spans
if isempty(options.ylims)
    options.ylims = max(abs(Zfull(:))) * 1.1;
end
if length(options.ylims) == 1
    if ~isempty(options.whichMarg)
        options.ylims = repmat(options.ylims, [1 max(options.whichMarg)]);
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot signifiance %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


gcf             = subplot('Position',[.13 .7 .40 .20]);


% dPC 1 (T X M)
acc_resh    = permute(accuracyShuffle,[3,2,1]);
m           = squeeze(nanmean(acc_resh(:,:,which_comp(1))));
s           = squeeze(nanstd(acc_resh(:,:,which_comp(1))));

plot(time,squeeze(accuracy(which_comp(1),1,:))-m,'LineWidth',2,'Color',margColours(4, :))
hold on

sigline   = nan(1,numel(time));
sigline(componentsSignif_95CI(1,:)==1) = -.09;

plot(time,sigline,'linewidth',4,'Color',margColours(4, :));
hline(0)

ylim([-.1 .3])
vline(0)
set(gca,'FontSize',20,'FontName','Arial')
xlabel('ripple time (sec)')
set(gca, 'XTickLabel', [])
ylabel('decoding perf')
set(gca,'TickDir','out')
box off;


% dPC 2 (Target assoc.)
hold on
m           = squeeze(nanmean(acc_resh(:,:,which_comp(2))));
s           = squeeze(nanstd(acc_resh(:,:,which_comp(2))));

plot(time,squeeze(accuracy(which_comp(2),1,:))-m,'LineWidth',2,'Color',margColours(1, :))

sigline   = nan(1,numel(time));
sigline(dpca_to_plot(2,:)==1) = -.075;

plot(time,sigline,'linewidth',4,'Color',margColours(1, :));
set(gca,'FontSize',20,'FontName','Arial')
xlabel('ripple time (sec)')
set(gca, 'XTickLabel', [])
ylabel('decoding perf')
set(gca,'TickDir','out')
box off;
set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left', ...
    'TickDir', 'out', 'LineWidth', 2);  

% dPC 4 (memory)
hold on
acc_resh    = permute(accuracyShuffle,[3,2,1]);
m           = squeeze(nanmean(acc_resh(:,:,which_comp(4))));
s           = squeeze(nanstd(acc_resh(:,:,which_comp(4))));


plot(time,squeeze(accuracy(which_comp(4),1,:))-m,'LineWidth',2,'color',margColours(2, :))
hold on

sigline   = nan(1,numel(time));
sigline(dpca_to_plot(4,:)==1) = -.06;

plot(time,sigline,'linewidth',4,'color',margColours(2, :));


vline(0)
set(gca,'FontSize',20,'FontName','Arial')
xlabel('ripple time (sec)')
ylabel('decoding perf')
set(gca, 'XTickLabel', [])
set(gca,'TickDir','out')
box off;
set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left', ...
    'TickDir', 'out', 'LineWidth', 2);  


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot scores %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Define component indices
comp_1 = 1;
comp_2 = 2;
comp_3 = 4;


component_ct1 = W(:, comp_1);  % First component capturing condition-time interaction
component_ct2 = W(:, comp_2);  % Second component capturing stimulus 
component_ct3 = W(:, comp_3);  % Third component capturing decision

time_before_idx = (time <= -0.4 & time >= -.8);


scores_ct1_before = nan(sum(time_before_idx), 4, 2);
scores_ct2_before = nan(sum(time_before_idx), 4, 2);
scores_ct3_before = nan(sum(time_before_idx), 4, 2);


for stim = 1:4
    for dec = 1:2
        % After the ripple
        data_after = squeeze(firingRatesAverage(:, stim, dec, time_before_idx));
        scores_ct1_before(:, stim, dec) = data_after' * component_ct1;
        scores_ct2_before(:, stim, dec) = data_after' * component_ct2;
        scores_ct3_before(:, stim, dec) = data_after' * component_ct3;
    end
end


colors = lines(4);

gcf = subplot('Position', [0.15 0.38 0.15 0.25]);
hold on;


for stim = 1:4
    % Correct decisions
    plot3(scores_ct1_before(:, stim, 1), scores_ct2_before(:, stim, 1), scores_ct3_before(:, stim, 1), ...
        'Color', colors(stim, :), 'LineWidth', 4, 'DisplayName', sprintf('Correct - Stimulus %d', stim));

    % Incorrect decisions
    plot3(scores_ct1_before(:, stim, 2), scores_ct2_before(:, stim, 2), scores_ct3_before(:, stim, 2), ...
        'Color', colors(stim, :), 'LineStyle', ':', 'LineWidth', 4, 'DisplayName', sprintf('Incorrect - Stimulus %d', stim));
end

xlabel(sprintf('dPC %d', comp_1));
ylabel(sprintf('dPC %d', comp_2));
zlabel(sprintf('dPC %d', comp_3));
% title('Pre ripple-events');
set(gca, 'FontSize', 16, 'LineWidth', 2, 'TickDir', 'out');
set(gca, 'FontSize', 20, 'FontName', 'Arial');
xlim([-200 200]);
ylim([-200 200]);
zlim([-200 200]);
grid on;
box on;
view(3)

hold off;



% after


component_ct1 = W(:, comp_1);  
component_ct2 = W(:, comp_2); 
component_ct3 = W(:, comp_3);


time_after_idx = (time >= 0.4 & time <= .8);  % Time points where any component is significant


scores_ct1_after = nan(sum(time_after_idx), 4, 2);  % Time points x Stimuli x Decisions
scores_ct2_after = nan(sum(time_after_idx), 4, 2);  % Time points x Stimuli x Decisions
scores_ct3_after = nan(sum(time_after_idx), 4, 2);  % Time points x Stimuli x Decisions


for stim = 1:4
    for dec = 1:2
        % After the ripple
        data_after = squeeze(firingRatesAverage(:, stim, dec, time_after_idx));
        scores_ct1_after(:, stim, dec) = data_after' * component_ct1;
        scores_ct2_after(:, stim, dec) = data_after' * component_ct2;
        scores_ct3_after(:, stim, dec) = data_after' * component_ct3;
    end
end


gcf = subplot('Position', [0.36 0.38 0.15 0.25]);
hold on;

for stim = 1:4
    % Correct decisions
    plot3(scores_ct1_after(:, stim, 1), scores_ct2_after(:, stim, 1), scores_ct3_after(:, stim, 1), ...
        'Color', colors(stim, :), 'LineWidth', 4, 'DisplayName', sprintf('Correct - Stimulus %d', stim));

    % Incorrect decisions
    plot3(scores_ct1_after(:, stim, 2), scores_ct2_after(:, stim, 2), scores_ct3_after(:, stim, 2), ...
        'Color', colors(stim, :), 'LineStyle', ':', 'LineWidth', 4, 'DisplayName', sprintf('Incorrect - Stimulus %d', stim));
end

xlabel(sprintf('dPC %d', comp_1));
ylabel(sprintf('dPC %d', comp_2));

% title('Post ripple-events');
set(gca, 'FontSize', 16, 'LineWidth', 2, 'TickDir', 'out');
set(gca, 'FontSize', 20, 'FontName', 'Arial');
xlim([-200 200]);
ylim([-200 200]);
zlim([-200 200]);
grid on;
box on;
view(3);
hold off;



% Calculate the Euclidean distances between points for each condition


distances_before = nan(4, 2,size(scores_ct1_after,1)-1); 
distances_after = nan(4, 2,size(scores_ct1_after,1)-1); 

for stim = 1:4
    for dec = 1:2
        % Before ripple
        points_before = [scores_ct1_before(:, stim, dec), scores_ct2_before(:, stim, dec), scores_ct3_before(:, stim, dec)];
%         distances_before(stim, dec) = mean(sqrt(sum(diff(points_before).^2, 2)));
        distances_before(stim, dec,:) = sqrt(sum(diff(points_before).^2, 2));

        % After ripple
        points_after = [scores_ct1_after(:, stim, dec), scores_ct2_after(:, stim, dec), scores_ct3_after(:, stim, dec)];
%         distances_after(stim, dec) = mean(sqrt(sum(diff(points_after).^2, 2)));
        distances_after(stim, dec,:) = sqrt(sum(diff(points_after).^2, 2));
    end
end


avg_distance_before = mean(distances_before, 'all');
avg_distance_after = mean(distances_after, 'all');

sem_distance_before = std(mean(mean((distances_before))))/sqrt(length((distances_before(:))));
sem_distance_after = std(mean(mean(distances_after)))/sqrt(length((distances_after(:))));

fprintf('Average Euclidean Distance (Before): %.2f\n', avg_distance_before);
fprintf('Average Euclidean Distance (After): %.2f\n', avg_distance_after);

[~,p_value_euc_dist,~,stats_euc] = ttest(distances_before(:),distances_after(:))


% kmeans and silhoutte score

features_after = [reshape(scores_ct1_after, [], 1), ...
    reshape(scores_ct2_after, [], 1), ...
    reshape(scores_ct3_after, [], 1)];

features_after = features_after(~any(isnan(features_after), 2), :);

features_before = [reshape(scores_ct1_before, [], 1), ...
    reshape(scores_ct2_before, [], 1), ...
    reshape(scores_ct3_before, [], 1)];

features_before = features_before(~any(isnan(features_before), 2), :);


num_clusters = 8;
[cluster_idx, cluster_centers] = kmeans(features_after, num_clusters, 'Replicates', 10);

% gcf             = subplot('Position',[.05 .10 .15 .25]);



max_clusters = 10;
silhouette_avg_after = zeros(max_clusters, 1);
silhouette_avg_before = zeros(max_clusters, 1);

for k = 1:max_clusters
    cluster_idx = kmeans(features_after, k, 'Replicates', 10);
    silhouette_vals = silhouette(features_after, cluster_idx);
    silhouette_avg_after(k) = mean(silhouette_vals);
end

% Find the peak silhouette score for AFTER
[val_after, peak_cluster_after] = max(silhouette_avg_after)

for k = 1:max_clusters
    cluster_idx = kmeans(features_before, k, 'Replicates', 10);
    silhouette_vals = silhouette(features_before, cluster_idx);
    silhouette_avg_before(k) = mean(silhouette_vals);
end

% Find the peak silhouette score for BEFORE
[val_before, peak_cluster_before] = max(silhouette_avg_before)



%% not in main figure

% Define component indices
comp_1 = 1;
comp_2 = 3;
comp_3 = 4;


component_ct1 = W(:, comp_1);  % First component capturing condition-time interaction
component_ct2 = W(:, comp_2);  % Second component capturing stimulus 
component_ct3 = W(:, comp_3);  % Third component capturing decision



time_before_idx = (time <= -0.4 & time >= -.8);

scores_ct1_before = nan(sum(time_before_idx), 4, 2);
scores_ct2_before = nan(sum(time_before_idx), 4, 2);
scores_ct3_before = nan(sum(time_before_idx), 4, 2);


for stim = 1:4
    for dec = 1:2
        % After the ripple
        data_after = squeeze(firingRatesAverage(:, stim, dec, time_before_idx));
        scores_ct1_before(:, stim, dec) = data_after' * component_ct1;
        scores_ct2_before(:, stim, dec) = data_after' * component_ct2;
        scores_ct3_before(:, stim, dec) = data_after' * component_ct3;
    end
end


colors = lines(4);



component_ct1 = W(:, comp_1);  
component_ct2 = W(:, comp_2); 
component_ct3 = W(:, comp_3);


time_after_idx = (time >= 0.4 & time <= .8);  % Time points where any component is significant


scores_ct1_after = nan(sum(time_after_idx), 4, 2);  % Time points x Stimuli x Decisions
scores_ct2_after = nan(sum(time_after_idx), 4, 2);  % Time points x Stimuli x Decisions
scores_ct3_after = nan(sum(time_after_idx), 4, 2);  % Time points x Stimuli x Decisions


for stim = 1:4
    for dec = 1:2
        % After the ripple
        data_after = squeeze(firingRatesAverage(:, stim, dec, time_after_idx));
        scores_ct1_after(:, stim, dec) = data_after' * component_ct1;
        scores_ct2_after(:, stim, dec) = data_after' * component_ct2;
        scores_ct3_after(:, stim, dec) = data_after' * component_ct3;
    end
end



% Calculate the Euclidean distances between points for each condition


distances_before = nan(4, 2,size(scores_ct1_after,1)-1); 
distances_after = nan(4, 2,size(scores_ct1_after,1)-1); 

for stim = 1:4
    for dec = 1:2
        % Before ripple
        points_before = [scores_ct1_before(:, stim, dec), scores_ct2_before(:, stim, dec), scores_ct3_before(:, stim, dec)];
%         distances_before(stim, dec) = mean(sqrt(sum(diff(points_before).^2, 2)));
        distances_before(stim, dec,:) = sqrt(sum(diff(points_before).^2, 2));

        % After ripple
        points_after = [scores_ct1_after(:, stim, dec), scores_ct2_after(:, stim, dec), scores_ct3_after(:, stim, dec)];
%         distances_after(stim, dec) = mean(sqrt(sum(diff(points_after).^2, 2)));
        distances_after(stim, dec,:) = sqrt(sum(diff(points_after).^2, 2));
    end
end


avg_distance_before = mean(distances_before, 'all');
avg_distance_after = mean(distances_after, 'all');

sem_distance_before = std(mean(mean((distances_before))))/sqrt(length((distances_before(:))));
sem_distance_after = std(mean(mean(distances_after)))/sqrt(length((distances_after(:))));

fprintf('Average Euclidean Distance (Before): %.2f\n', avg_distance_before);
fprintf('Average Euclidean Distance (After): %.2f\n', avg_distance_after);

[~,p_value_euc_dist,~,stats_euc] = ttest(distances_before(:),distances_after(:))


% kmeans and silhoutte score

features_after = [reshape(scores_ct1_after, [], 1), ...
    reshape(scores_ct2_after, [], 1), ...
    reshape(scores_ct3_after, [], 1)];

features_after = features_after(~any(isnan(features_after), 2), :);

features_before = [reshape(scores_ct1_before, [], 1), ...
    reshape(scores_ct2_before, [], 1), ...
    reshape(scores_ct3_before, [], 1)];

features_before = features_before(~any(isnan(features_before), 2), :);


num_clusters = 8;
[cluster_idx, cluster_centers] = kmeans(features_after, num_clusters, 'Replicates', 10);

% gcf             = subplot('Position',[.05 .10 .15 .25]);



max_clusters = 10;
silhouette_avg_after = zeros(max_clusters, 1);
silhouette_avg_before = zeros(max_clusters, 1);

for k = 1:max_clusters
    cluster_idx = kmeans(features_after, k, 'Replicates', 10);
    silhouette_vals = silhouette(features_after, cluster_idx);
    silhouette_avg_after(k) = mean(silhouette_vals);
end

% Find the peak silhouette score for AFTER
[val_after, peak_cluster_after] = max(silhouette_avg_after)

for k = 1:max_clusters
    cluster_idx = kmeans(features_before, k, 'Replicates', 10);
    silhouette_vals = silhouette(features_before, cluster_idx);
    silhouette_avg_before(k) = mean(silhouette_vals);
end

% Find the peak silhouette score for BEFORE
[val_before, peak_cluster_before] = max(silhouette_avg_before)




%% Do it per 200ms

% Define component indices
comp_1 = 1;
comp_2 = 2;
comp_3 = 4;

component_ct1 = W(:, comp_1);  % First component capturing condition-time interaction
component_ct2 = W(:, comp_2);  % Second component capturing stimulus 
component_ct3 = W(:, comp_3);  % Third component capturing decision

% Define time windows
window_size = 0.2; % 200 ms
window_edges = -1:window_size:1; % From -1000 ms to +1000 ms
num_windows = length(window_edges) - 1;

% Initialize figure

colors = lines(4);
distances = nan(4, 2, num_windows); % Store distances for each window

for w = 1:num_windows
    % Identify time indices within the current window
    time_idx = (time >= window_edges(w) & time < window_edges(w+1));
    
    scores_ct1 = nan(sum(time_idx), 4, 2);
    scores_ct2 = nan(sum(time_idx), 4, 2);
    scores_ct3 = nan(sum(time_idx), 4, 2);
  
    for stim = 1:4
        for dec = 1:2
            data_window = squeeze(firingRatesAverage(:, stim, dec, time_idx));
            scores_ct1(:, stim, dec) = data_window' * component_ct1;
            scores_ct2(:, stim, dec) = data_window' * component_ct2;
            scores_ct3(:, stim, dec) = data_window' * component_ct3;
            
            % Compute distances between consecutive points
            points_window = [scores_ct1(:, stim, dec), scores_ct2(:, stim, dec), scores_ct3(:, stim, dec)];
            distances(stim, dec, w) = mean(sqrt(sum(diff(points_window).^2, 2)), 'omitnan');
        end
    end
end

% Plot average and separate correct/incorrect distances over time

hold on;

gcf = subplot('Position', [0.13 0.07 0.4 0.2]);

times = window_edges(1:end-1) + window_size/2;
avg_distances = squeeze(mean(distances, [1 2], 'omitnan'));

plot(times, avg_distances, 'k', 'LineWidth', 2, 'DisplayName', 'Average');

xlabel('Time (s)');
ylabel('Cluster Distance');
title('Cluster Distance Over Time');
legend;
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'TickDir', 'out');

grid off;

vline(0)

