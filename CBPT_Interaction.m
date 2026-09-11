% Eyes x Phase INTERACTION on the aperiodic exponent
%   Delta(EO-EC) = (EO-EC)_post - (EO-EC)_pre, per drug
% =========================================================================
% Within-subject double difference per (Subject, Drug):
%   INT = [exponent_open_post - exponent_closed_post]
%       - [exponent_open_pre  - exponent_closed_pre ]
% i.e. the drug-induced CHANGE in the eyes-open-vs-closed contrast.
% One-sample cluster test of INT against 0, per drug, with the SAME omnibus
% max-statistic procedure as Analysis 3 of the EO-EC script:
%   - clusters on a 55 mm adjacency graph; cluster mass = sum of t;
%   - synchronized sign-flip permutation (one shared sign per participant
%     across their drugs), 5000 iterations;
%   - SINGLE absolute null distribution POOLED across all drugs and both 
%     tails (omnibus FWER across drugs + space + both tails);
%   - evaluated directly at alpha = 0.05;
%   - p = (1 + sum(null_abs >= abs(obs)))/(1 + num_perms) [Phipson & Smyth 2010];
%   - cluster mass indexed find(mask)(bins==b) [never mask(bins==b)].
%
% A subject contributes to a drug only if all four sessions exist for that
% drug (open/closed x pre/post). Cohen's d_z is a selection-biased descriptor
% of the cluster, not an unbiased population estimate (Meyer et al., 2021).
%
% EXCLUSIONS: nimodipine and placebo are dropped before the family is formed
%   (see exclude_drugs in Section 1). Dropping placebo removes the only control
%   for session-order / time-of-day / arousal drift in Delta(EO-EC), so a
%   surviving cluster cannot be separated from a pure pre-to-post session
%   effect. Dropping any arm also shrinks the omnibus family, which lowers the
%   permutation maximum and therefore lowers the p-values of the arms that
%   remain. The exclusion is not neutral - report it in the methods.
%
% INPUT: *_parameterised.mat written by APOP_Glu_Specparam.ipynb, read through
%   load_parameterised.m - recursive over the BIDS tree, maps eyestate EO/EC to
%   open/closed, drops sessions with no drug assignment, and errors if any
%   file's channel labels differ in content or order from the first file's
%   (the adjacency graph is built once and would otherwise be misaligned).
%
% REPORTING: every candidate cluster is written out with its per-arm N and a
%   Significant flag. The p-value column is P_FWER - already corrected by the
%   max-statistic null.
% =========================================================================
clearvars; close all; clc;
eeglab_path = "D:\Linus\MATLAB_applications\eeglab2026.0.0";
addpath(eeglab_path)

% paths
inPath  = "D:\Linus\APOP_Glu\parameterised2_BIDS";   % = OUT_ROOT in APOP_Glu_Specparam.ipynb
outPath = "D:\Linus\APOP_Glu\Analysis4\CBPT_Interaction_QC";
if ~exist(outPath, 'dir'); mkdir(outPath); end

% --- Thresholds ---
cft_p      = 0.05;      % cluster-DEFINING threshold (two-tailed p)
alpha      = 0.05;      % two-sided FWER (tested directly against absolute max dist)
num_perms  = 5000;
min_N      = 3;         % minimum subjects required to keep a drug in the family
rng(42)

%% 1. Load ALL data (both eyes)
% Drugs held out of this analysis. See EXCLUSIONS in the header: this changes
% both what the test can claim and the p-values of the arms that remain.
exclude_drugs = ["nimodipine", "placebo"];

[records, n_nodrug] = load_parameterised(inPath);
n_loaded = numel(records);
records  = records(~ismember([records.Drug], lower(exclude_drugs)));
if isempty(records)
    error('Every record was excluded - check exclude_drugs against the data.');
end
fprintf(['Loaded %d record(s) (both eyes); %d file(s) dropped for no drug ' ...
         'assignment; %d excluded (%s).\n'], numel(records), n_nodrug, ...
         n_loaded - numel(records), strjoin(cellstr(exclude_drugs), ', '));

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
        n_skip = n_skip + 1; continue;
    end
    eo = records(io(1)).Exponent(:);
    ec = records(ic(1)).Exponent(:);
    eoec(end+1).UniqueID = records(io(1)).UniqueID;
    eoec(end).Drug    = records(io(1)).Drug;
    eoec(end).PrePost = records(io(1)).PrePost;
    eoec(end).EOEC    = (eo - ec)';            % 1 x nChan
end

if isempty(eoec); error('No EO-EC sessions could be formed (missing open/closed pairs).'); end
fprintf('Formed %d EO-EC sessions; skipped %d session(s) missing an eyes condition.\n', numel(eoec), n_skip);

uu = [eoec.UniqueID];  dd = [eoec.Drug];  pp = [eoec.PrePost];
EO = vertcat(eoec.EOEC);            % nSessions x nChan

%% 3. Build the interaction: INT = (EO-EC)_post - (EO-EC)_pre per (Subject, Drug)
key2  = uu + "|" + dd;             % subject|drug
ukey2 = unique(key2);

int_rec = struct('UniqueID', {}, 'Drug', {}, 'INT', {});
n_incomplete = 0;
for i = 1:numel(ukey2)
    rows  = find(key2 == ukey2(i));
    rpre  = rows(pp(rows) == "pre");
    rpost = rows(pp(rows) == "post");
    if isempty(rpre) || isempty(rpost)
        n_incomplete = n_incomplete + 1; continue;   % need pre AND post EO-EC
    end
    int_rec(end+1).UniqueID = uu(rpre(1));
    int_rec(end).Drug = dd(rpre(1));
    int_rec(end).INT  = EO(rpost(1),:) - EO(rpre(1),:);   % (EO-EC)post - (EO-EC)pre
end

if isempty(int_rec); error('No (Subject,Drug) had both pre and post EO-EC.'); end
fprintf('Built %d interaction maps; %d (subject,drug) lacked a complete pre+post pair.\n', ...
    numel(int_rec), n_incomplete);

iu = [int_rec.UniqueID];  idr = [int_rec.Drug];  II = vertcat(int_rec.INT);

%% 4. Spatial adjacency graph (55 mm on the standard 10-05 template)
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

chanlocs_plot = struct('labels', my_labels);
chanlocs_plot = pop_chanedit(chanlocs_plot, 'lookup', master_locs_file);
for i = 1:numel(chanlocs_plot)
    if isfield(chanlocs_plot,'X') && ~isempty(chanlocs_plot(i).X)
        ox = chanlocs_plot(i).X; oy = chanlocs_plot(i).Y;
        chanlocs_plot(i).X = oy; chanlocs_plot(i).Y = -ox;
    end
end
chanlocs_plot = pop_chanedit(chanlocs_plot, 'convert', {'cart2all'});

%% 5. Per-drug interaction CBPT, OMNIBUS max-statistic across drugs
fprintf('\n=== INTERACTION: Delta(EO-EC) = (EO-EC)post - (EO-EC)pre, per drug ===\n');
drugs_all = unique(idr);
drug_matrices = struct('Name', {}, 'D', {}, 'SubIDs', {}, 'N', {}, 'tcrit', {});

for di = 1:numel(drugs_all)
    drug = drugs_all(di);
    sel  = (idr == drug);
    if nnz(sel) < min_N
        fprintf('  %-16s skipped (N = %d < %d)\n', drug, nnz(sel), min_N); continue;
    end
    drug_matrices(end+1).Name = drug;
    drug_matrices(end).D      = II(sel,:);
    drug_matrices(end).SubIDs = iu(sel)';
    drug_matrices(end).N      = nnz(sel);
    drug_matrices(end).tcrit  = tinv(1 - cft_p/2, nnz(sel) - 1);
    fprintf('  %-16s N = %d\n', drug, nnz(sel));
end

if isempty(drug_matrices)
    error('No drug had >= %d complete subjects for the interaction.', min_N);
end

% Function call corrected to pass the labels properly
[obs_maps, drug_names, clusters] = run_omnibus_maxstat(drug_matrices, G_full, num_perms, ...
    "Post > Pre (EO-EC increased)", "Post < Pre (EO-EC decreased)");

% Evaluated against full alpha directly because absolute max covers both tails
is_sig = arrayfun(@(c) c.P_FWER <= alpha, clusters);

if isempty(clusters)
    fprintf('No candidate clusters formed in any drug.\n');
else
    rep = build_report_omnibus(clusters, is_sig, my_labels);
    writetable(rep, fullfile(outPath, 'Interaction_EOEC_post_minus_pre_omnibus.csv'));
    fprintf('%d candidate cluster(s); %d survive omnibus FWER (alpha = %.3f).\n', ...
        numel(clusters), sum(is_sig), alpha);
    disp(rep);
end

% significant channels per drug for plotting
sig_cell = repmat({[]}, numel(drug_names), 1);
for k = 1:numel(clusters)
    if is_sig(k)
        j = clusters(k).DrugIdx;
        sig_cell{j} = [sig_cell{j}, clusters(k).Chans];
    end
end

titles = cellstr("Delta(EO-EC): " + string(drug_names));
draw_topos(obs_maps, titles, sig_cell, chanlocs_plot, ...
    fullfile(outPath, 'Interaction_EOEC_post_minus_pre_topo.png'));

fprintf('\nDone. Outputs in: %s\n', outPath);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

% Added pos_label and neg_label variables to function definition
function [post_t, post_names, clusters] = run_omnibus_maxstat(drug_matrices, G_full, num_perms, pos_label, neg_label)
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
                % Updated to use dynamic pos_label
                clusters(end+1) = struct('Drug', drug_matrices(d).Name, 'Direction', pos_label, 'N', Nn, ...
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
                % Updated to use dynamic neg_label
                clusters(end+1) = struct('Drug', drug_matrices(d).Name, 'Direction', neg_label, 'N', Nn, ...
                    'ClusterMass', cm, 'Chans', cc, 'P_FWER', pv, 'Cohens_d', mean(scd)/std(scd), 'DrugIdx', d);
            end
        end
    end
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