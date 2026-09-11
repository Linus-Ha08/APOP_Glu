% Eyes-Open vs Eyes-Closed (EO - EC) CBPT on the aperiodic exponent
% =========================================================================
% Within-subject contrast throughout: EO - EC (exponent_open - exponent_closed),
% computed per (Subject, Drug, PrePost) session. Three analyses:
%
%   ANALYSIS 1  PRE, single CBPT, EO-EC pooled across subjects.
%               One value per subject = mean EO-EC over that subject's PRE
%               sessions (any drug). One-sample cluster test against 0.
%               -> "Is there an EO-EC effect at baseline?"
%
%   ANALYSIS 2  PRE, difference-of-differences, single CBPT.
%               Per subject: (EO-EC)_active,pre  -  (EO-EC)_placebo,pre, where
%               active = mean over the subject's active-drug PRE sessions and
%               placebo = the subject's placebo PRE session. One-sample
%               (paired) cluster test on subjects having BOTH.
%               -> "Do active and placebo differ at baseline?" (expect no).
%               This is the correct direct test; comparing a significant
%               active map to a non-significant placebo map is NOT a test of
%               their difference (Nieuwenhuis et al. 2011; Gelman & Stern 2006).
%
%   ANALYSIS 3  POST, EO-EC per drug, OMNIBUS max-statistic across drugs.
%               Per drug: one-sample EO-EC cluster test on the post sessions;
%               the permutation null pools the maximum absolute cluster mass 
%               ACROSS ALL drugs (synchronized sign-flip preserving cross-drug 
%               dependence).
%
% ENGINES (local functions):
%   - cbpt_core        : single-group absolute max-statistic (Analyses 1 & 2).
%   - run_omnibus_maxstat : pooled-across-drugs absolute max-statistic
%                           with synchronized sign-flip (Analysis 3).
%   Both: one-sample t per channel; clusters on a 55 mm adjacency graph;
%   cluster mass = sum of t; cluster mass indexed find(mask)(bins==b);
%   p = (1 + sum(null_abs >= abs(obs)))/(1 + num_perms) [Phipson & Smyth 2010].
%
% INFERENCE LEVELS:
%   - By using the absolute maximum cluster mass for the permutation null, 
%     the distribution inherently covers both tails. Significance for a 
%     two-sided test is therefore evaluated directly at alpha = 0.05.
%
% INPUT: *_parameterised.mat written by APOP_Glu_Specparam.ipynb, read through
%   load_parameterised.m - recursive over the BIDS tree, maps eyestate EO/EC to
%   open/closed, drops sessions with no drug assignment, and errors if any
%   file's channel labels differ in content or order from the first file's
%   (the adjacency graph is built once and would otherwise be misaligned).
%
% ASSUMPTIONS TO CHECK (top of script):
%   - active_drugs and placebo_label must match the Drug column of
%     APOP_Glu_Conditions.csv (compared lower-case; the script warns if absent).
%   - A subject's "active" value averages their active-drug sessions so a
%     person is counted once, not 3-4x.
%
% REPORTING: every candidate cluster is written out with its per-arm N and a
%   Significant flag, so near-misses and small-N arms stay visible. The p-value
%   column is P_FWER: it is already corrected by the max-statistic null. Earlier
%   versions named it RawP, which read as an uncorrected p-value. It never was.
%
% REPRODUCIBILITY: rng seeded, so a rerun reproduces the same p-values.
%
% EFFECT SIZE: Cohen's d_z on the test-selected electrodes is a
%   selection-biased descriptor of the cluster, not an unbiased population
%   estimate (Meyer et al., 2021).
% =========================================================================
clearvars; close all; clc;
eeglab_path = "D:\Linus\MATLAB_applications\eeglab2026.0.0";
addpath(eeglab_path)

% paths
inPath  = "D:\Linus\APOP_Glu\parameterised_BIDS";   % = OUT_ROOT in APOP_Glu_Specparam.ipynb
outPath = "D:\Linus\APOP_Glu\CBPT_eyes";
if ~exist(outPath, 'dir'); mkdir(outPath); end

% --- Conditions (must match the Drug column of APOP_Glu_Conditions.csv) ---
active_drugs  = lower(["dextromethorphan","nimodipine","perampanel"]);   % APOP_Glu panel
% active_drugs  = lower(["alprazolam","baclofen","diazepam","zolpidem"]); % APOP_GABA panel
placebo_label = "placebo";

% --- Thresholds ---
cft_p      = 0.05;      % cluster-DEFINING threshold (two-tailed p)
alpha      = 0.05;      % two-sided FWER (tested against absolute max distribution)
num_perms  = 5000;
min_N      = 3;         % minimum subjects required to run a test
rng(42)                 % reproducible permutations

%% 1. Load ALL data (both eyes)
[records, n_nodrug] = load_parameterised(inPath);
fprintf('Loaded %d records (both eyes conditions); %d file(s) dropped for no drug assignment.\n', ...
    numel(records), n_nodrug);

my_labels = records(1).Chanlabels;
num_chans = numel(my_labels);

%% 2. Build EO - EC per (Subject, Drug, PrePost)
ru = [records.UniqueID]; rd = [records.Drug]; rp = [records.PrePost]; reye = [records.Eyes];
keys  = ru + "|" + rd + "|" + rp;
ukeys = unique(keys);

eoec = struct('UniqueID', {}, 'Drug', {}, 'PrePost', {}, 'EOEC', {});
n_skip = 0;
for i = 1:numel(ukeys)
    idxs = find(keys == ukeys(i));
    io = idxs(reye(idxs) == "open");
    ic = idxs(reye(idxs) == "closed");
    
    if isempty(io) || isempty(ic)
        n_skip = n_skip + 1; continue;          % need both eyes to form EO-EC
    end
    
    eo = records(io(1)).Exponent(:);
    ec = records(ic(1)).Exponent(:);
    
    eoec(end+1).UniqueID = records(io(1)).UniqueID;
    eoec(end).Drug    = records(io(1)).Drug;
    eoec(end).PrePost = records(io(1)).PrePost;
    eoec(end).EOEC    = (eo - ec)';             % 1 x nChan
end

if isempty(eoec); error('No EO-EC sessions could be formed (missing open/closed pairs).'); end
fprintf('Formed %d EO-EC sessions; skipped %d session(s) missing an eyes condition.\n', numel(eoec), n_skip);

% Vectorized handles aligned with the rows of EO
uu = [eoec.UniqueID];  dd = [eoec.Drug];  pp = [eoec.PrePost];
EO = vertcat(eoec.EOEC);            % nSessions x nChan

% Sanity: are the requested conditions present?
present = unique(dd);
for a = active_drugs
    if ~any(present == a); warning('Active drug "%s" not found in data.', a); end
end
if ~any(present == placebo_label); warning('Placebo label "%s" not found in data.', placebo_label); end

%% 3. Spatial adjacency graph (55 mm on the standard 10-05 template)
eeglab_root = fileparts(which('eeglab'));
master_locs_file = fullfile(eeglab_root, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc');
if ~exist(master_locs_file, 'file'); error('standard_1005.elc not found at %s', master_locs_file); end

eeglab nogui;
chanlocs = struct('labels', my_labels);
chanlocs = pop_chanedit(chanlocs, 'lookup', master_locs_file);

coords   = [[chanlocs.X]', [chanlocs.Y]', [chanlocs.Z]'];
dist_matrix = squareform(pdist(coords));
max_dist_threshold = 55;
adj_mat = (dist_matrix <= max_dist_threshold) & (dist_matrix > 0);
neighbor_counts = sum(adj_mat, 2);

fprintf('Adjacency: %d channels; neighbours min %d / mean %.2f / max %d.\n', ...
    num_chans, min(neighbor_counts), mean(neighbor_counts), max(neighbor_counts));

G_full = graph(adj_mat);

% Channel locations for topoplots (axis convention used previously)
chanlocs_plot = struct('labels', my_labels);
chanlocs_plot = pop_chanedit(chanlocs_plot, 'lookup', master_locs_file);
for i = 1:numel(chanlocs_plot)
    if isfield(chanlocs_plot,'X') && ~isempty(chanlocs_plot(i).X)
        ox = chanlocs_plot(i).X; oy = chanlocs_plot(i).Y;
        chanlocs_plot(i).X = oy; chanlocs_plot(i).Y = -ox;
    end
end
chanlocs_plot = pop_chanedit(chanlocs_plot, 'convert', {'cart2all'});

subj_ids = unique(uu);

%% ========================================================================
%% ANALYSIS 1: PRE, single EO-EC CBPT, pooled across subjects
%% ========================================================================
fprintf('\n=== ANALYSIS 1: EO-EC at baseline (PRE), pooled ===\n');
D1 = zeros(0, num_chans);
for i = 1:numel(subj_ids)
    sel = (uu == subj_ids(i)) & (pp == "pre");
    if any(sel); D1(end+1,:) = mean(EO(sel,:), 1); end
end
fprintf('N subjects (PRE): %d\n', size(D1,1));

if size(D1,1) >= min_N
    [obs_t1, cl1] = cbpt_core(D1, G_full, num_perms, cft_p, "EO > EC", "EO < EC");
    % Compare against alpha directly because we use an absolute max distribution
    sig1 = arrayfun(@(c) c.P_FWER <= alpha, cl1);
    sigch1 = [];  for k = 1:numel(cl1); if sig1(k); sigch1 = [sigch1, cl1(k).Chans]; end; end
    
    rep1 = build_report_single(cl1, sig1, my_labels);
    if ~isempty(rep1)
        writetable(rep1, fullfile(outPath, 'A1_PRE_EOEC_pooled.csv'));
        disp(rep1);
    else
        fprintf('No candidate clusters formed.\n');
    end
    
    draw_topos({obs_t1}, "PRE EO-EC (pooled)", {sigch1}, chanlocs_plot, ...
        fullfile(outPath, 'A1_PRE_EOEC_pooled_topo.png'));
else
    warning('Analysis 1 skipped: fewer than %d subjects.', min_N);
end

%% ========================================================================
%% ANALYSIS 2: PRE, difference-of-differences (active - placebo), single CBPT
%% ========================================================================
fprintf('\n=== ANALYSIS 2: (EO-EC)_active - (EO-EC)_placebo at baseline (PRE) ===\n');
D2 = zeros(0, num_chans);
for i = 1:numel(subj_ids)
    sa = (uu == subj_ids(i)) & (pp == "pre") & ismember(dd, active_drugs);
    sp = (uu == subj_ids(i)) & (pp == "pre") & (dd == placebo_label);
    if any(sa) && any(sp)
        D2(end+1,:) = mean(EO(sa,:), 1) - mean(EO(sp,:), 1);
    end
end
fprintf('N subjects with both active & placebo PRE: %d\n', size(D2,1));

if size(D2,1) >= min_N
    [obs_t2, cl2] = cbpt_core(D2, G_full, num_perms, cft_p, "Active > Placebo", "Active < Placebo");
    % Compare against alpha directly because we use an absolute max distribution
    sig2 = arrayfun(@(c) c.P_FWER <= alpha, cl2);
    sigch2 = [];  for k = 1:numel(cl2); if sig2(k); sigch2 = [sigch2, cl2(k).Chans]; end; end
    
    rep2 = build_report_single(cl2, sig2, my_labels);
    if ~isempty(rep2)
        writetable(rep2, fullfile(outPath, 'A2_PRE_active_minus_placebo.csv'));
        disp(rep2);
    else
        fprintf('No candidate clusters formed (no baseline active-vs-placebo difference detected).\n');
    end
    
    draw_topos({obs_t2}, "PRE (Active - Placebo) EO-EC", {sigch2}, chanlocs_plot, ...
        fullfile(outPath, 'A2_PRE_active_minus_placebo_topo.png'));
else
    warning('Analysis 2 skipped: fewer than %d subjects with both.', min_N);
end

%% ========================================================================
%% ANALYSIS 3: POST, EO-EC per drug, OMNIBUS max-statistic across drugs
%% ========================================================================
fprintf('\n=== ANALYSIS 3: EO-EC (POST) per drug, omnibus max-statistic ===\n');
drugs_post = unique(dd);

% Build per-drug POST matrices
drug_matrices = struct('Name', {}, 'D', {}, 'SubIDs', {}, 'N', {}, 'tcrit', {});
for di = 1:numel(drugs_post)
    drug = drugs_post(di);
    sel  = (dd == drug) & (pp == "post");
    if nnz(sel) < min_N
        fprintf('  %-16s skipped (N = %d < %d)\n', drug, nnz(sel), min_N); continue;
    end
    
    drug_matrices(end+1).Name = drug;
    drug_matrices(end).D      = EO(sel,:);
    drug_matrices(end).SubIDs = uu(sel)';            % column of subject IDs
    drug_matrices(end).N      = nnz(sel);
    drug_matrices(end).tcrit  = tinv(1 - cft_p/2, nnz(sel) - 1);
    fprintf('  %-16s N = %d\n', drug, nnz(sel));
end

if isempty(drug_matrices)
    warning('Analysis 3 skipped: no drug had >= %d post subjects.', min_N);
else
    [post_t, post_names, cl3] = run_omnibus_maxstat(drug_matrices, G_full, num_perms);
    
    is_sig3 = arrayfun(@(c) c.P_FWER <= alpha, cl3);
    
    if isempty(cl3)
        fprintf('No candidate clusters formed in any post drug condition.\n');
    else
        rep3 = build_report_omnibus(cl3, is_sig3, my_labels);
        writetable(rep3, fullfile(outPath, 'A3_POST_perdrug_omnibus_maxstat.csv'));
        fprintf('%d candidate cluster(s); %d survive omnibus FWER (alpha = %.3f).\n', ...
            numel(cl3), sum(is_sig3), alpha);
        disp(rep3);
    end
    
    % significant channels per drug for plotting
    sig_cell3 = repmat({[]}, numel(post_names), 1);
    for k = 1:numel(cl3)
        if is_sig3(k)
            j = cl3(k).DrugIdx;
            sig_cell3{j} = [sig_cell3{j}, cl3(k).Chans];
        end
    end
    
    % Ensure post_names is treated as strings, concatenate, then convert to a cell array
    titles3 = cellstr("POST EO-EC: " + string(post_names)); 
    
    draw_topos(post_t, titles3, sig_cell3, chanlocs_plot, ...
        fullfile(outPath, 'A3_POST_perdrug_topo.png'));
end

fprintf('\nDone. Outputs in: %s\n', outPath);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function [obs_t, clusters] = cbpt_core(D, G_full, num_perms, cft_p, pos_label, neg_label)
% Single-group one-sample (sign-flip) cluster-permutation test on a
% Subjects x Channels difference matrix D. Single absolute max-statistic null.
    [N, ~] = size(D);
    tcrit = tinv(1 - cft_p/2, N - 1);
    
    obs_t = (sum(D,1) ./ N) ./ (std(D,0,1) ./ sqrt(N));     % 1 x nchan
    
    null_abs = zeros(num_perms,1);
    
    for p = 1:num_perms
        s  = randi([0,1], N, 1) * 2 - 1;
        pD = D .* s;
        t  = (sum(pD,1) ./ N) ./ (std(pD,0,1) ./ sqrt(N));
        
        iter_max_abs = 0; 
        
        pm = t > tcrit;
        if any(pm)
            pidx = find(pm); bins = conncomp(subgraph(G_full, pidx));
            for b = 1:max(bins)
                cm = sum(t(pidx(bins==b))); 
                if abs(cm) > iter_max_abs
                    iter_max_abs = abs(cm); 
                end
            end
        end
        
        nm = t < -tcrit;
        if any(nm)
            nidx = find(nm); bins = conncomp(subgraph(G_full, nidx));
            for b = 1:max(bins)
                cm = sum(t(nidx(bins==b))); 
                if abs(cm) > iter_max_abs
                    iter_max_abs = abs(cm); 
                end
            end
        end
        
        null_abs(p) = iter_max_abs;
    end
    
    clusters = struct('Direction', {}, 'N', {}, 'ClusterMass', {}, 'Chans', {}, 'P_FWER', {}, 'Cohens_d', {});
    
    pm = obs_t > tcrit;
    if any(pm)
        pidx = find(pm); bins = conncomp(subgraph(G_full, pidx));
        for b = 1:max(bins)
            cc = pidx(bins==b); cm = sum(obs_t(cc));
            pv = (1 + sum(null_abs >= abs(cm))) / (1 + num_perms);
            scd = mean(D(:,cc), 2);
            clusters(end+1) = struct('Direction', pos_label, 'N', N, 'ClusterMass', cm, ...
                'Chans', cc, 'P_FWER', pv, 'Cohens_d', mean(scd)/std(scd));
        end
    end
    
    nm = obs_t < -tcrit;
    if any(nm)
        nidx = find(nm); bins = conncomp(subgraph(G_full, nidx));
        for b = 1:max(bins)
            cc = nidx(bins==b); cm = sum(obs_t(cc));
            pv = (1 + sum(null_abs >= abs(cm))) / (1 + num_perms);
            scd = mean(D(:,cc), 2);
            clusters(end+1) = struct('Direction', neg_label, 'N', N, 'ClusterMass', cm, ...
                'Chans', cc, 'P_FWER', pv, 'Cohens_d', mean(scd)/std(scd));
        end
    end
end

function [post_t, post_names, clusters] = run_omnibus_maxstat(drug_matrices, G_full, num_perms)
% Omnibus max-statistic across drugs, single tail (absolute max).
% Synchronized sign-flip: one shared sign per participant across their drugs,
% preserving cross-drug dependence. The null pools the absolute maximum 
% cluster mass ACROSS ALL drugs and both tails. Each observed cluster is
% referred to that pooled absolute extremum distribution. tcrit is read per drug.
    nd = numel(drug_matrices);
    active_unique_IDs = unique(vertcat(drug_matrices.SubIDs));
    nsub = numel(active_unique_IDs);
    
    global_max_abs = zeros(num_perms, 1);
    
    for p = 1:num_perms
        s_all = randi([0,1], nsub, 1) * 2 - 1;
        iter_max_abs = 0;
        
        for d = 1:nd
            [~, loc] = ismember(drug_matrices(d).SubIDs, active_unique_IDs);
            s  = s_all(loc);
            Nn = drug_matrices(d).N;
            pD = drug_matrices(d).D .* s;
            t  = (sum(pD,1) ./ Nn) ./ (std(pD,0,1) ./ sqrt(Nn));
            tcrit = drug_matrices(d).tcrit;
            
            pm = t > tcrit;
            if any(pm)
                pidx = find(pm); bins = conncomp(subgraph(G_full, pidx));
                for b = 1:max(bins)
                    cm = sum(t(pidx(bins==b))); 
                    if abs(cm) > iter_max_abs
                        iter_max_abs = abs(cm); 
                    end
                end
            end
            
            nm = t < -tcrit;
            if any(nm)
                nidx = find(nm); bins = conncomp(subgraph(G_full, nidx));
                for b = 1:max(bins)
                    cm = sum(t(nidx(bins==b))); 
                    if abs(cm) > iter_max_abs
                        iter_max_abs = abs(cm); 
                    end
                end
            end
        end
        global_max_abs(p) = iter_max_abs;
    end
    
    post_t = cell(nd,1); post_names = strings(nd,1);
    clusters = struct('Drug', {}, 'Direction', {}, 'N', {}, 'ClusterMass', {}, 'Chans', {}, 'P_FWER', {}, 'Cohens_d', {}, 'DrugIdx', {});
    
    for d = 1:nd
        D = drug_matrices(d).D; Nn = drug_matrices(d).N; tcrit = drug_matrices(d).tcrit;
        obs_t = (sum(D,1) ./ Nn) ./ (std(D,0,1) ./ sqrt(Nn));
        post_t{d} = obs_t; post_names(d) = drug_matrices(d).Name;
        
        pm = obs_t > tcrit;
        if any(pm)
            pidx = find(pm); bins = conncomp(subgraph(G_full, pidx));
            for b = 1:max(bins)
                cc = pidx(bins==b); cm = sum(obs_t(cc));
                pv = (1 + sum(global_max_abs >= abs(cm))) / (1 + num_perms);
                scd = mean(D(:,cc), 2);
                clusters(end+1) = struct('Drug', drug_matrices(d).Name, 'Direction', "EO > EC", 'N', Nn, ...
                    'ClusterMass', cm, 'Chans', cc, 'P_FWER', pv, 'Cohens_d', mean(scd)/std(scd), 'DrugIdx', d);
            end
        end
        
        nm = obs_t < -tcrit;
        if any(nm)
            nidx = find(nm); bins = conncomp(subgraph(G_full, nidx));
            for b = 1:max(bins)
                cc = nidx(bins==b); cm = sum(obs_t(cc));
                pv = (1 + sum(global_max_abs >= abs(cm))) / (1 + num_perms);
                scd = mean(D(:,cc), 2);
                clusters(end+1) = struct('Drug', drug_matrices(d).Name, 'Direction', "EO < EC", 'N', Nn, ...
                    'ClusterMass', cm, 'Chans', cc, 'P_FWER', pv, 'Cohens_d', mean(scd)/std(scd), 'DrugIdx', d);
            end
        end
    end
end

function T = build_report_single(clusters, sigflags, my_labels)
% P_FWER is the max-statistic corrected p-value, not a raw per-cluster p.
    m = numel(clusters);
    if m == 0, T = table(); return; end
    Direction = strings(m,1); N = zeros(m,1); ClusterMass = zeros(m,1); P_FWER = zeros(m,1);
    Significant = false(m,1); NumChannels = zeros(m,1); Channels = strings(m,1); Cohens_d = zeros(m,1);
    for k = 1:m
        Direction(k)   = clusters(k).Direction;
        N(k)           = clusters(k).N;
        ClusterMass(k) = clusters(k).ClusterMass;
        P_FWER(k)      = clusters(k).P_FWER;
        Significant(k) = sigflags(k);
        NumChannels(k) = numel(clusters(k).Chans);
        Channels(k)    = strjoin(my_labels(clusters(k).Chans), ', ');
        Cohens_d(k)    = clusters(k).Cohens_d;
    end
    T = table(Direction, N, ClusterMass, P_FWER, Significant, NumChannels, Channels, Cohens_d);
    T = sortrows(T, 'P_FWER');
end

function T = build_report_omnibus(clusters, is_sig, my_labels)
    m = numel(clusters);
    if m == 0, T = table(); return; end
    Drug = strings(m,1); Direction = strings(m,1); N = zeros(m,1); ClusterMass = zeros(m,1);
    P_FWER = zeros(m,1); Significant = false(m,1);
    NumChannels = zeros(m,1); Channels = strings(m,1); Cohens_d = zeros(m,1);
    for k = 1:m
        Drug(k)        = clusters(k).Drug;
        Direction(k)   = clusters(k).Direction;
        N(k)           = clusters(k).N;
        ClusterMass(k) = clusters(k).ClusterMass;
        P_FWER(k)      = clusters(k).P_FWER;
        Significant(k) = is_sig(k);
        NumChannels(k) = numel(clusters(k).Chans);
        Channels(k)    = strjoin(my_labels(clusters(k).Chans), ', ');
        Cohens_d(k)    = clusters(k).Cohens_d;
    end
    T = table(Drug, Direction, N, ClusterMass, P_FWER, Significant, NumChannels, Channels, Cohens_d);
    T = sortrows(T, 'P_FWER');
end

function draw_topos(t_maps, titles, sig_cell, chanlocs_plot, savepath)
    n = numel(t_maps);
    allt = [];  for i = 1:n; allt = [allt, t_maps{i}]; end
    max_t = max(abs(allt));
    if isempty(max_t) || isnan(max_t) || max_t == 0; max_t = 3; end
    figure('Color', 'w', 'Position', [100, 100, max(500, 380*min(n,4)), 380*ceil(n/4)]);
    tiledlayout('flow');
    for i = 1:n
        nexttile;
        sc = sig_cell{i};
        if isempty(sc)
            topoplot(t_maps{i}, chanlocs_plot, 'electrodes', 'off', 'maplimits', [-max_t max_t]);
        else
            topoplot(t_maps{i}, chanlocs_plot, 'electrodes', 'off', ...
                'emarker2', {sc, '.', 'k', 15, 1}, 'maplimits', [-max_t max_t]);
        end
        if isstring(titles) || ischar(titles)
            title(titles, 'Interpreter', 'none');
        else
            title(titles(i), 'Interpreter', 'none');
        end
    end
    cb = colorbar; cb.Layout.Tile = 'east';
    saveas(gcf, savepath);
end