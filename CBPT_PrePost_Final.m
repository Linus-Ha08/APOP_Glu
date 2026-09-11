% Max-Statistic / Synchronized Sign-Flip CBPT for Pre/Post Aperiodic Exponent
% =========================================================================
% CORRECTED VERSION. Changes vs. the previous script:
%
%   [FIX 1 - CRITICAL] Permutation (null) cluster mass was indexed with
%       pos_mask(bins==b) / neg_mask(bins==b). bins is defined over the
%       SUBGRAPH nodes, so this maps subgraph-node positions onto channel
%       indices and collapses the null cluster mass toward ~0, making the
%       corrected p-values severely anti-conservative. It now uses
%       pos_idx(bins==b) / neg_idx(bins==b), exactly matching Section 5.
%
%   [FIX 2] Two-sided FWER. Positive and negative clusters are referred to a
%       SINGLE permutation distribution of the maximum absolute cluster mass
%       (taken across all drugs and both tails) and tested at full alpha.
%       This controls the two-sided error rate at alpha = 0.05 while
%       accounting for the dependence between the two tails.
%
%   [FIX 3] Permutation p-values use (1 + #>=obs) / (1 + num_perms)
%       (Phipson & Smyth, 2010) so they can never be exactly zero.
%
%   [MINOR] Removed a duplicated channel-lookup block in Section 3; preallocated
%       SubIDs with strings(); added an explicit tiledlayout for the topoplots.
%
%   [FIX 4] The cluster-DEFINING threshold is now its own parameter (cft_p).
%       It used to be read off `alpha`, so tightening the family-wise level to
%       0.01 silently redefined which channels form clusters - a different
%       test, not a stricter one. cft_p and alpha are now independent.
%
%   [FIX 5] rng seed set, so a rerun reproduces the same p-values. Without one,
%       a cluster sitting near alpha can flip between runs.
%
%   [FIX 6] The report carries the per-arm N. A surviving cluster from N = 3 and
%       one from N = 15 used to be indistinguishable in the CSV.
%
% INPUT: *_parameterised.mat written by APOP_Glu_Specparam.ipynb, read through
%   load_parameterised.m - recursive over the BIDS tree, maps eyestate EO/EC to
%   open/closed, drops sessions with no drug assignment, and errors if any
%   file's channel labels differ in content or order from the first file's
%   (the adjacency graph is built once and would otherwise be misaligned).
%
% SCOPE: FWER is controlled across all drug conditions and both directions
%   WITHIN one eyes condition. Eyes-open and eyes-closed are run separately
%   (set target_eyes) and are NOT mutually corrected.
%
% EFFECT SIZE: Cohen's d_z is computed on the electrodes the test selected
%   and is therefore a (selection-biased) descriptor of the identified
%   cluster, not an unbiased population estimate (Meyer et al., 2021).
%
% PERMUTATION RESOLUTION: signs are drawn once per participant and shared
%   across that participant's drug conditions, so the relevant permutation
%   space is 2^(number of unique participants), not 2^(N per drug). With
%   ~15 per drug but many unique participants pooled, 5000 Monte-Carlo
%   permutations give an adequate estimate.
% =========================================================================

clearvars; close all; clc;
eeglab_path = "D:\Linus\MATLAB_applications\eeglab2026.0.0";
addpath(eeglab_path)

% paths
inPath  = "D:\Linus\APOP_Glu\parameterised2_BIDS";   % = OUT_ROOT in APOP_Glu_Specparam.ipynb
outPath = "D:\Linus\APOP_Glu\Analysis\CBPT_PrePost_Final";
if ~exist(outPath, 'dir')
    mkdir(outPath);
end

%%% TARGET CONDITION  (run once with 'open', once with 'closed')
target_eyes = 'open';

cft_p = 0.05;   % cluster-DEFINING threshold (two-tailed p) - independent of alpha
alpha = 0.05;   % two-sided family-wise level
rng(42)         % reproducible permutations

%% 1. Load Data (one eyes condition)
[records, n_nodrug] = load_parameterised(inPath);
records = records(strcmpi([records.Eyes], target_eyes));
if isempty(records)
    error('No %s-eyes recordings found under %s', target_eyes, inPath);
end
fprintf("Data loaded. %d %s-eyes record(s); %d file(s) dropped for no drug assignment.\n", ...
    numel(records), target_eyes, n_nodrug);

my_labels = records(1).Chanlabels;
num_chans = numel(my_labels);

%% 2. Compare pre vs post and compute differences
all_uniqueIDs = [records.UniqueID];
all_drugs     = [records.Drug];
all_prepost   = [records.PrePost];

unique_drugs  = unique(all_drugs);   % load_parameterised already dropped unassigned drugs
unique_IDs    = unique(all_uniqueIDs);

diff_data = struct('UniqueID', {}, 'Drug', {}, 'DiffExp', {});
diff_idx  = 1;

for d = 1:length(unique_drugs)
    current_drug = unique_drugs(d);
    for s = 1:length(unique_IDs)
        current_id = unique_IDs(s);

        idx_pre  = find(all_uniqueIDs == current_id & all_drugs == current_drug & strcmpi(all_prepost, 'pre'));
        idx_post = find(all_uniqueIDs == current_id & all_drugs == current_drug & strcmpi(all_prepost, 'post'));

        if ~isempty(idx_pre) && ~isempty(idx_post)
            diff_exp = records(idx_post(1)).Exponent(:) - records(idx_pre(1)).Exponent(:);

            diff_data(diff_idx).UniqueID = current_id;
            diff_data(diff_idx).Drug     = current_drug;
            diff_data(diff_idx).DiffExp  = diff_exp'; % row vector (1 x nChan), Post - Pre
            diff_idx = diff_idx + 1;
        end
    end
end
fprintf("Differences calculated.\n");

%% 3. Setup Spatial Adjacency Matrix using EEGLAB native coordinates
eeglab_path = fileparts(which('eeglab'));
master_locs_file = fullfile(eeglab_path, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc');
if ~exist(master_locs_file, 'file')
    error('Could not find standard_1005.elc at: %s', master_locs_file);
end

eeglab nogui;

% Look up coordinates for the montage established in Section 1
chanlocs = struct('labels', my_labels);
chanlocs = pop_chanedit(chanlocs, 'lookup', master_locs_file);

% Spatial neighbourhood matrix from pure geometry
X = [chanlocs.X]';
Y = [chanlocs.Y]';
Z = [chanlocs.Z]';
coords = [X, Y, Z];

% Distance between every channel pair
dist_matrix = squareform(pdist(coords));

% Proximity threshold (units of the template). 40-45 mm ~ a 10-20
% neighbourhood; 55 mm is more permissive (larger clusters).
max_dist_threshold = 55;

% Adjacency: true if channels are close, but not 0 (self)
adj_mat = (dist_matrix <= max_dist_threshold) & (dist_matrix > 0);

% Neighbour counts -> REPORT THESE in the methods
neighbor_counts = sum(adj_mat, 2);
fprintf('Min neighbors: %d\nMax neighbors: %d\nMean neighbors: %.2f\n', ...
    min(neighbor_counts), max(neighbor_counts), mean(neighbor_counts));

% Graph object for the permutation engine
G_full = graph(adj_mat);
fprintf('Adjacency matrix built with %d channels.\n', num_chans);

%% 4. Max-Statistic Permutation Engine (synchronized sign-flip)
num_perms = 5000;

% Extract data matrices per drug for fast vectorized computation
drug_matrices = struct();
for d = 1:length(unique_drugs)
    idx = find(string({diff_data.Drug}) == unique_drugs(d));
    drug_matrices(d).Name = unique_drugs(d);
    drug_matrices(d).N    = length(idx);

    % D is Subjects x Channels
    D_mat   = zeros(length(idx), num_chans);
    sub_ids = strings(length(idx), 1);
    for i = 1:length(idx)
        D_mat(i, :) = diff_data(idx(i)).DiffExp;
        sub_ids(i)  = diff_data(idx(i)).UniqueID;
    end
    drug_matrices(d).D      = D_mat;
    drug_matrices(d).SubIDs = sub_ids;
    % Two-tailed critical t for CLUSTER FORMATION (cft_p, independent of alpha)
    drug_matrices(d).tcrit  = tinv(1 - (cft_p/2), length(idx) - 1);
end

fprintf('Starting permutations (%d iterations)...\n', num_perms);
global_max_abs = zeros(num_perms, 1);   % max |cluster mass| across drugs AND both tails

% Unique participants across the whole dataset -> one shared sign each
% (preserves dependence between drug conditions that share participants)
active_unique_IDs = unique(vertcat(drug_matrices.SubIDs));
num_active_subs   = length(active_unique_IDs);

wb = waitbar(0, 'Running permutations...');

for p = 1:num_perms
    if mod(p, 100) == 0; waitbar(p/num_perms, wb); end

    % Random signs for all subjects globally (+1 or -1), shared across drugs
    sub_signs = randi([0, 1], num_active_subs, 1) * 2 - 1;

    iter_max_abs = 0;

    for d = 1:length(drug_matrices)
        % Map global signs to the subjects in this drug condition
        [~, loc] = ismember(drug_matrices(d).SubIDs, active_unique_IDs);
        current_signs = sub_signs(loc);

        % Apply permutation (flip differences)
        perm_D = drug_matrices(d).D .* current_signs;

        % Fast one-sample (dependent) t-test per channel
        mean_D  = sum(perm_D, 1) / drug_matrices(d).N;
        std_D   = std(perm_D, 0, 1);
        t_stats = mean_D ./ (std_D / sqrt(drug_matrices(d).N));

        tcrit = drug_matrices(d).tcrit;

        % ---- Positive clusters (FIXED indexing) ----
        pos_mask = t_stats > tcrit;
        if any(pos_mask)
            pos_idx = find(pos_mask);
            G_sub   = subgraph(G_full, pos_idx);
            bins    = conncomp(G_sub);
            for b = 1:max(bins)
                cluster_chans = pos_idx(bins == b);
                cluster_mass  = sum(t_stats(cluster_chans));
                if abs(cluster_mass) > iter_max_abs
                    iter_max_abs = abs(cluster_mass);
                end
            end
        end

        % ---- Negative clusters (FIXED indexing) ----
        neg_mask = t_stats < -tcrit;
        if any(neg_mask)
            neg_idx = find(neg_mask);
            G_sub   = subgraph(G_full, neg_idx);
            bins    = conncomp(G_sub);
            for b = 1:max(bins)
                cluster_chans = neg_idx(bins == b);
                cluster_mass  = sum(t_stats(cluster_chans));
                if abs(cluster_mass) > iter_max_abs
                    iter_max_abs = abs(cluster_mass);
                end
            end
        end
    end

    % Largest absolute cluster mass across ALL drugs and both tails
    global_max_abs(p) = iter_max_abs;
end
close(wb);
fprintf('Permutations complete.\n');

%% 5. Evaluate Observed Data against the Max-Statistic Distribution
cbpt_report = table([], [], [], [], [], [], [], [], ...
    'VariableNames', {'Drug', 'Direction', 'N', 'ClusterMass', 'Omnitest_P_Value', 'NumChannels', 'Channels', 'Cohens_d'});
stats_results = struct();

for d = 1:length(drug_matrices)
    % Unpermuted (observed) data
    raw_D  = drug_matrices(d).D;                 % (Post - Pre) differences
    mean_D = sum(raw_D, 1) / drug_matrices(d).N;
    std_D  = std(raw_D, 0, 1);
    obs_t_stats = mean_D ./ (std_D / sqrt(drug_matrices(d).N));
    tcrit = drug_matrices(d).tcrit;

    sig_channels = [];

    % ---- Positive clusters ----
    pos_mask = obs_t_stats > tcrit;
    if any(pos_mask)
        pos_idx = find(pos_mask);
        G_sub   = subgraph(G_full, pos_idx);
        bins    = conncomp(G_sub);
        for b = 1:max(bins)
            cluster_chans = pos_idx(bins == b);
            cluster_mass  = sum(obs_t_stats(cluster_chans));

            % FWER-corrected p (Phipson & Smyth: +1 num & denom)
            p_val = (1 + sum(global_max_abs >= abs(cluster_mass))) / (1 + num_perms);

            % Two-sided control via single absolute-mass distribution
            if p_val <= alpha
                % Cohen's d_z (selection-biased descriptor of the cluster)
                subj_cluster_diffs = mean(raw_D(:, cluster_chans), 2);
                cohens_d = mean(subj_cluster_diffs) / std(subj_cluster_diffs);

                sig_channels = [sig_channels, cluster_chans];
                ch_names = strjoin(my_labels(cluster_chans), ', ');
                new_row = table(string(drug_matrices(d).Name), string("Increase (Post > Pre)"), ...
                    drug_matrices(d).N, cluster_mass, p_val, length(cluster_chans), string(ch_names), cohens_d, ...
                    'VariableNames', cbpt_report.Properties.VariableNames);
                cbpt_report = [cbpt_report; new_row];
            end
        end
    end

    % ---- Negative clusters ----
    neg_mask = obs_t_stats < -tcrit;
    if any(neg_mask)
        neg_idx = find(neg_mask);
        G_sub   = subgraph(G_full, neg_idx);
        bins    = conncomp(G_sub);
        for b = 1:max(bins)
            cluster_chans = neg_idx(bins == b);
            cluster_mass  = sum(obs_t_stats(cluster_chans));

            % FWER-corrected p (single absolute-mass distribution)
            p_val = (1 + sum(global_max_abs >= abs(cluster_mass))) / (1 + num_perms);

            % Two-sided control via single absolute-mass distribution
            if p_val <= alpha
                % Cohen's d_z (selection-biased descriptor of the cluster)
                subj_cluster_diffs = mean(raw_D(:, cluster_chans), 2);
                cohens_d = mean(subj_cluster_diffs) / std(subj_cluster_diffs);

                sig_channels = [sig_channels, cluster_chans];
                ch_names = strjoin(my_labels(cluster_chans), ', ');
                new_row = table(string(drug_matrices(d).Name), string("Decrease (Post < Pre)"), ...
                    drug_matrices(d).N, cluster_mass, p_val, length(cluster_chans), string(ch_names), cohens_d, ...
                    'VariableNames', cbpt_report.Properties.VariableNames);
                cbpt_report = [cbpt_report; new_row];
            end
        end
    end

    stats_results(d).Drug      = drug_matrices(d).Name;
    stats_results(d).t_stats   = obs_t_stats';
    stats_results(d).sig_chans = sig_channels;
end

% Save and display report
if ~isempty(cbpt_report)
    csv_name = fullfile(outPath, sprintf('Omnitest_CBPT_Report_Eyes_%s.csv', target_eyes));
    writetable(cbpt_report, csv_name);
    disp('--- FWER-CORRECTED SIGNIFICANT CLUSTERS (two-sided alpha = 0.05) ---');
    disp(cbpt_report);
else
    fprintf('No significant clusters surviving FWER correction.\n');
end

%% 6. Create Topoplots
chanlocs = struct('labels', my_labels);
chanlocs = pop_chanedit(chanlocs, 'lookup', master_locs_file);
for i = 1:length(chanlocs)
    if isfield(chanlocs, 'X') && ~isempty(chanlocs(i).X)
        old_X = chanlocs(i).X; old_Y = chanlocs(i).Y;
        chanlocs(i).X = old_Y; chanlocs(i).Y = -old_X;
    end
end
chanlocs = pop_chanedit(chanlocs, 'convert', {'cart2all'});

all_t_vals = [];
for d = 1:length(stats_results)
    all_t_vals = [all_t_vals; stats_results(d).t_stats];
end
max_t = max(abs(all_t_vals));
if isempty(max_t) || isnan(max_t) || max_t == 0; max_t = 3; end

figure('Color', 'w', 'Position', [100, 100, 1600, 800]);
tiledlayout('flow');
for d = 1:length(stats_results)
    nexttile;
    t_vals    = stats_results(d).t_stats;
    sig_chans = stats_results(d).sig_chans;

    if isempty(sig_chans)
        topoplot(t_vals, chanlocs, 'electrodes', 'off', 'maplimits', [-max_t max_t]);
    else
        topoplot(t_vals, chanlocs, 'electrodes', 'off', ...
            'emarker2', {sig_chans, '.', 'k', 15, 1}, 'maplimits', [-max_t max_t]);
    end
    title(sprintf('Drug: %s (t-values)', stats_results(d).Drug), 'Interpreter', 'none');
end
cb_t = colorbar; cb_t.Layout.Tile = 'east';
saveNameT = fullfile(outPath, sprintf('Topoplot_Tvals_Eyes_%s_Omnitest.png', target_eyes));
saveas(gcf, saveNameT);
fprintf('Topoplots saved to: %s\n', saveNameT);