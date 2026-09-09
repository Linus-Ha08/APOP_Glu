% Final APOP_Glu preprocessing script 09/09/2026
% Windows Version

% to check stil
% line noise clean 
% ICLabel thresholds

% takes raw data and performs the following preprocessing steps:
% 1. Load raw BrainVision data (BIDS layout from bidsify_apop_glu.py)
% 2. Identify EOG channels by label (aux channels '31'/'32'), not by index
% 3. Resample to 1000 Hz (raw data must already be >= 1000 Hz)
% 4. Filter: 80 Hz lowpass, 0.5 Hz highpass
% 5. Line noise removal: Zapline-plus + CleanLine (50 Hz)
% 6. clean_rawdata: reject channels/segments, ASR-correct bursts (its own
%    highpass is disabled -- the explicit 0.5 Hz filter above is the only
%    highpass applied, so the data isn't filtered twice)
% 7. ICA (extended infomax) + ICLabel, remove flagged components
% 8. Interpolate EEG channels back to the full pre-clean_rawdata set
% 9. Re-reference to average (FCz recovered); EOG channels are excluded
%    from the reference computation
%
% Per-file design variables:
%   PrePost  -- 'pre' / 'post' (acq- entity in the filename)
%   EyeState -- 'EO' / 'EC'    (task-rest entity in the filename)
%   Condition, Drug -- looked up by Subject + Day in APOP_Glu_Conditions.csv;
%                      anything without a unique match is flagged 'undefined'
% All four are written to the log file and to EEG.etc of each saved dataset.
%
% returns continuous, cleaned data

%% Initialize
clear; clc;
eeglabRoot  = "D:\Linus\MATLAB_applications\eeglab2026.0.0";
ICLabelRoot = "D:\Linus\MATLAB_applications\eeglab2026.0.0\plugins\ICLabel";
rawpath     = "D:\Linus\APOP_Glu\raw_BIDS";
CondFile  = "D:\Linus\APOP_Glu\scripts\APOP_Glu_Conditions.csv";
outpath     = "D:\Linus\Local\APOP\preprocessed_BIDS";

% labels of EOG chans
eogRawLabels = {'31', '32'};
eogNewLabels = {'EOG_L', 'EOG_R'};

if ~exist(outpath, 'dir')
    mkdir(outpath);
end

% logfile
LogFile = fullfile(outpath, 'APOP_Glu_preprocessing_log.csv');

%% Log setup
% one struct per file, turned into a table row with struct2table.
logFieldDefaults = struct( ...
    'Subject_ID',           "", ...
    'Raw_Filename',         "", ...
    'Day',                  NaN, ...
    'PrePost',              "", ...
    'EyeState',             "", ...
    'Condition',            "undefined", ...
    'Drug',                 "undefined", ...
    'Status',               "", ...
    'Error_Message',        "", ...
    'Loaded_Srate',         NaN, ...
    'Loaded_Points',        NaN, ...
    'Loaded_Seconds',       NaN, ...
    'Loaded_Channels',      NaN, ...
    'EEG_Channels',         NaN, ...
    'EOG_Channels',         NaN, ...
    'Resampled_Srate',      NaN, ...
    'Removed_Points',       NaN, ...
    'Removed_Seconds',      NaN, ...
    'Removed_Channels_N',   NaN, ...
    'Removed_Channels',     "", ...
    'Total_Removed_ICs',    NaN, ...
    'Eye_Rejected_ICs',     NaN, ...
    'Muscle_Rejected_ICs',  NaN, ...
    'Heart_Rejected_ICs',   NaN, ...
    'Line_Rejected_ICs',    NaN, ...
    'Channel_Rejected_ICs', NaN);

if exist(LogFile, 'file')
    logData = readtable(LogFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    logData = alignLogTable(logData, logFieldDefaults, LogFile);
else
    logData = struct2table(logFieldDefaults);
    logData(1,:) = [];
end

%% Condition / Drug lookup table
if exist(CondFile, 'file') ~= 2
    error('Condition file not found: %s', CondFile);
end

condTbl   = readtable(CondFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
condNames = strtrim(erase(string(condTbl.Properties.VariableNames), char(65279)));  % trim + strip BOM

iSubj = find(strcmpi(condNames, "Subject"),   1);
iDay  = find(strcmpi(condNames, "Day"),       1);
iCond = find(strcmpi(condNames, "Condition"), 1);
iDrug = find(strcmpi(condNames, "Drug"),      1);
if isempty(iSubj) || isempty(iDay) || isempty(iCond) || isempty(iDrug)
    error('%s must contain Subject, Day, Condition and Drug columns (found: %s).', ...
        CondFile, strjoin(condNames, ', '));
end

condSubject = strtrim(toString(condTbl{:, iSubj}));
condDay     = toNumeric(condTbl{:, iDay});
condCond    = strtrim(toString(condTbl{:, iCond}));
condDrug    = strtrim(toString(condTbl{:, iDrug}));

% empty / missing cells in the table itself count as undefined
condCond(ismissing(condCond) | condCond == "") = "undefined";
condDrug(ismissing(condDrug) | condDrug == "") = "undefined";

fprintf('Loaded %d Condition/Drug entries from %s\n', height(condTbl), CondFile);

addpath(eeglabRoot);
addpath(ICLabelRoot);
eeglab nogui;

%% Loop through raw data files in BIDS format
files = dir(fullfile(rawpath, 'sub-*', 'ses-day*', '*', '*', '*.vhdr'));
files = files(~[files.isdir]);

fprintf('Found %d files in %s\n', numel(files), rawpath);

if isempty(files)
    error('No files found in %s. Exiting.', rawpath);
end

%% Process each file
for iFile = 1:numel(files)
    rawFile   = files(iFile).name;
    rawFolder = files(iFile).folder;
    fprintf('\nProcessing file %d/%d: %s\n', iFile, numel(files), rawFile);

    % Parse BIDS entities from the filename (bidsify_apop_glu.py's FILENAME_TEMPLATE)
    tok = regexp(rawFile, ['^sub-(?<subject>[^_]+)_ses-day(?<day>\d+)_' ...
        'task-rest(?<eyestate>EO|EC)_acq-(?<prepost>pre|post)_eeg\.vhdr$'], 'names', 'once');
    if isempty(tok)
        error('Filename does not match the expected BIDS pattern: %s', rawFile);
    end

    [~, base_name, ~] = fileparts(rawFile);
    savedSetName = [base_name '_preprocessed.set'];
    savedSetPath = fullfile(outpath, savedSetName);

    if exist(savedSetPath, 'file') == 2
        fprintf('File %s already processed. Skipping...\n', savedSetName);
        continue;
    end

    logEntry = logFieldDefaults;
    logEntry.Subject_ID   = string(tok.subject);
    logEntry.Raw_Filename = string(rawFile);
    logEntry.Day          = str2double(tok.day);
    logEntry.PrePost      = string(tok.prepost);
    logEntry.EyeState     = string(tok.eyestate);

    %% Condition / Drug for this subject and day
    % defaults are already "undefined"; only a unique match overwrites them
    isMatch = strcmpi(condSubject, strtrim(logEntry.Subject_ID)) & condDay == logEntry.Day;
    nMatch  = nnz(isMatch);
    if nMatch == 1
        logEntry.Condition = condCond(isMatch);
        logEntry.Drug      = condDrug(isMatch);
    elseif nMatch == 0
        warning('APOP:NoCondition', ...
            'No Condition/Drug entry for subject %s day %d in %s -- flagged as undefined.', ...
            logEntry.Subject_ID, logEntry.Day, CondFile);
    else
        warning('APOP:AmbiguousCondition', ...
            '%d Condition/Drug entries for subject %s day %d in %s -- flagged as undefined.', ...
            nMatch, logEntry.Subject_ID, logEntry.Day, CondFile);
    end
    fprintf('  Subject %s | Day %d | PrePost %s | EyeState %s | Condition %s | Drug %s\n', ...
        logEntry.Subject_ID, logEntry.Day, logEntry.PrePost, logEntry.EyeState, ...
        logEntry.Condition, logEntry.Drug);

    try
        %% Load raw data
        EEG = pop_loadbv(rawFolder, rawFile);
        EEG.setname = char(logEntry.Subject_ID);

        logEntry.Loaded_Srate    = EEG.srate;
        logEntry.Loaded_Points   = EEG.pnts;
        logEntry.Loaded_Seconds  = EEG.xmax;
        logEntry.Loaded_Channels = EEG.nbchan;

        %% chanlocs & EOG chans 
        EEG = pop_chanedit(EEG, 'lookup', fullfile(eeglabRoot, 'plugins', 'dipfit', ...
            'standard_BEM', 'elec', 'standard_1005.elc'));

        % The aux EOG channels are found by their raw labels ('31'/'32'),
        % whatever position they sit in, and renamed to EOG_L / EOG_R.
        % Done after the template lookup so they keep no scalp coordinates.
        for iEOG = 1:numel(eogRawLabels)
            eogIdx = find(strcmpi(strtrim({EEG.chanlocs.labels}), eogRawLabels{iEOG}));
            if numel(eogIdx) ~= 1
                error('Expected exactly one channel labelled ''%s'', found %d in %s.', ...
                    eogRawLabels{iEOG}, numel(eogIdx), rawFile);
            end
            EEG.chanlocs(eogIdx).labels = eogNewLabels{iEOG};
            EEG.chanlocs(eogIdx).type   = 'EOG';
        end
        EEG = eeg_checkset(EEG);

        eogChanIdx = find(strcmpi({EEG.chanlocs.type}, 'EOG'));
        logEntry.EEG_Channels = EEG.nbchan - numel(eogChanIdx);
        logEntry.EOG_Channels = numel(eogChanIdx);

        %% resample 
        if EEG.srate < 1000
            error('Sampling rate is below 1000 Hz. Please check the data.');
        end
        if EEG.srate > 1000
            EEG = pop_resample(EEG, 1000);
        end
        logEntry.Resampled_Srate = EEG.srate;

        %% filter
        % lowpass at 80 Hz, highpass at 0.5 Hz (clean_rawdata's own highpass
        % is disabled below, so the data is only highpassed once, here)
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 80);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);

        % line noise removal
        EEG = pop_zapline_plus(EEG, 'noisefreqs', 50, 'plotResults', 0);
        EEG = pop_cleanline(EEG, 'linefreqs', 50, 'plotfigures', 0);

        %% clean_rawdata
        % record pre cleaning
        original_chanlocs    = EEG.chanlocs;
        before_clean_points  = EEG.pnts;
        before_clean_seconds = EEG.pnts / EEG.srate;
        chans_before_clean   = {EEG.chanlocs.labels};

        % clean data 
        EEG = pop_clean_rawdata(EEG, ...
            'FlatlineCriterion', 5, ...     % reject chans if flat > 5 s
            'ChannelCriterion', 0.7, ...    % bad if chan correlated less than this to own reconstruction based on other chans; default = 0.8
            'LineNoiseCriterion', 4, ...    % more line noise relative to its signal than this value in sd
            'Highpass', 'off', ...          % highpass already applied explicitly above
            'BurstCriterion', 20, ...       % mark activity bursts that are 20 SDs above the mean (for later rejection); default = 5
            'WindowCriterion', 0.25, ...    % after ASR criterion: reject windows where more than 25% of channels show bursts
            'BurstRejection', 'off', ...    % on = reject marked bursts, off = ASR will attempt to clean them
            'Distance', 'Euclidian', ...
            'WindowCriterionTolerances', '[-Inf 7]');
        EEG = eeg_checkset(EEG);

        % record post cleaning
        after_clean_points  = EEG.pnts;
        after_clean_seconds = EEG.pnts / EEG.srate;
        chans_after_clean   = {EEG.chanlocs.labels};

        removed_chans     = setdiff(chans_before_clean, chans_after_clean);
        removed_chans_str = strjoin(removed_chans, ', ');
        if isempty(removed_chans_str)
            removed_chans_str = 'None';
        end

        logEntry.Removed_Points     = before_clean_points - after_clean_points;
        logEntry.Removed_Seconds    = before_clean_seconds - after_clean_seconds;
        logEntry.Removed_Channels_N = numel(removed_chans);
        logEntry.Removed_Channels   = string(removed_chans_str);

        %% ICA & ICLabel 
        EEG = pop_runica(EEG, 'extended', 1);
        EEG = pop_iclabel(EEG, 'default');

        % ICLabel class order: 1 brain, 2 muscle, 3 eye, 4 heart,
        %                      5 line noise, 6 channel noise, 7 other
        % reject if > 0.5 for muscle, eye, heart, channel noise, or > 0.8 for line noise
        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;

        EEG = pop_icflag(EEG, [NaN NaN; ...
                    0.50 1;                 % muscle
                    0.50 1;                 % eye
                    0.50 1;                 % heart
                    0.80 1;                 % line noise
                    0.50 1;                 % channel noise
                    NaN NaN]);

        rej_mask = logical(EEG.reject.gcompreject(:));
        rejICs   = find(rej_mask);

        muscle_rejICs  = find(rej_mask & ic_probs(:,2) >= 0.50);
        eye_rejICs     = find(rej_mask & ic_probs(:,3) >= 0.50);
        heart_rejICs   = find(rej_mask & ic_probs(:,4) >= 0.50);
        line_rejICs    = find(rej_mask & ic_probs(:,5) >= 0.80);
        channel_rejICs = find(rej_mask & ic_probs(:,6) >= 0.50);

        logEntry.Total_Removed_ICs    = numel(rejICs);
        logEntry.Eye_Rejected_ICs     = numel(eye_rejICs);
        logEntry.Muscle_Rejected_ICs  = numel(muscle_rejICs);
        logEntry.Heart_Rejected_ICs   = numel(heart_rejICs);
        logEntry.Line_Rejected_ICs    = numel(line_rejICs);
        logEntry.Channel_Rejected_ICs = numel(channel_rejICs);

        % reject ICs
        if ~isempty(rejICs)
            EEG = pop_subcomp(EEG, rejICs, 0);
        end

        % wipe ICA metadata
        EEG.icaweights  = [];
        EEG.icasphere   = [];
        EEG.icawinv     = [];
        EEG.icachansind = [];
        EEG.icaact      = [];
        EEG.reject.gcompreject = [];
        if isfield(EEG.etc, 'ic_classification')
            EEG.etc.ic_classification = [];
        end

        %% interpolate chans
        % Only  non-EOG channels are interpolated
        original_eeg_chanlocs = original_chanlocs(~strcmpi({original_chanlocs.type}, 'EOG'));
        EEG = pop_interp(EEG, original_eeg_chanlocs, 'spherical');
        EEG = eeg_checkset(EEG);

        %% re-reference to avg. and recover active reference (FCz)
        % EOG channels are excluded from the reference computation
        eogChanIdx = find(strcmpi({EEG.chanlocs.type}, 'EOG'));

        tmpl = readlocs(fullfile(eeglabRoot, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc'));
        fcz_loc = tmpl(find(strcmpi({tmpl.labels}, 'FCz'), 1));
        if ~isfield(fcz_loc, 'type'), fcz_loc.type = ''; end
        if ~isfield(fcz_loc, 'ref'), fcz_loc.ref = ''; end
        if ~isfield(fcz_loc, 'urchan'), fcz_loc.urchan = []; end

        EEG = pop_reref(EEG, [], 'refloc', fcz_loc, 'exclude', eogChanIdx);

        %% save metadata 
        EEG.etc.subject_id          = char(logEntry.Subject_ID);
        EEG.etc.raw_filename        = char(logEntry.Raw_Filename);
        EEG.etc.day                 = logEntry.Day;
        EEG.etc.prepost             = char(logEntry.PrePost);    % 'pre' / 'post'
        EEG.etc.eye_state           = char(logEntry.EyeState);   % 'EO' / 'EC'
        EEG.etc.condition           = char(logEntry.Condition);  % drug condition 1-4, or 'undefined'
        EEG.etc.drug                = char(logEntry.Drug);       % drug name, or 'undefined'
        EEG.etc.original_chanlocs   = original_chanlocs;
        EEG.etc.removed_points      = logEntry.Removed_Points;
        EEG.etc.removed_seconds     = logEntry.Removed_Seconds;
        EEG.etc.removed_channels_n  = logEntry.Removed_Channels_N;
        EEG.etc.removed_channels    = char(logEntry.Removed_Channels);
        EEG.etc.rejected_ICs        = rejICs;
        EEG.etc.rejected_ICs_counts = struct( ...
            'total',   logEntry.Total_Removed_ICs, ...
            'eye',     logEntry.Eye_Rejected_ICs, ...
            'muscle',  logEntry.Muscle_Rejected_ICs, ...
            'heart',   logEntry.Heart_Rejected_ICs, ...
            'line',    logEntry.Line_Rejected_ICs, ...
            'channel', logEntry.Channel_Rejected_ICs);

        %% save preprocessed dataset
        EEG = pop_saveset(EEG, 'filename', savedSetName, 'filepath', outpath);

        logEntry.Status        = string('Success');
        logEntry.Error_Message = string('');
        fprintf('Finished processing file %d/%d: %s\n', iFile, numel(files), rawFile);

    catch ME
        logEntry.Status        = string('Error');
        logEntry.Error_Message = string(ME.message);
        fprintf(2, 'ERROR in file %s\n%s\n', rawFile, ME.message);
    end

    %% update logfile 
    logData = [logData; struct2table(logEntry)];
    writetable(logData, LogFile);
    fprintf('Log file updated: %s\n', LogFile);
end

disp('Finished preprocessing.');

%% Local helper functions
function s = toString(x)
% column of any table type -> string column
    if isstring(x)
        s = x;
    elseif isnumeric(x)
        s = strings(numel(x), 1);
        s(~isnan(x)) = string(x(~isnan(x)));
    else
        s = string(x);
    end
    s = s(:);
end

function v = toNumeric(x)
% column of any table type -> double column
    if isnumeric(x)
        v = double(x);
    else
        v = str2double(string(x));
    end
    v = v(:);
end

function T = alignLogTable(T, defaults, LogFile)
% Bring an existing log file onto the current schema: add missing columns,
% drop unknown ones (e.g. the old Phase/Condition naming), fix column order
% and column types so struct2table rows can be appended.
    % Migrate the previous naming so old rows keep their meaning:
    % Phase -> PrePost, and the old Condition (EO/EC) -> EyeState. The new
    % Condition/Drug columns are then filled with 'undefined' for those rows.
    if ismember('Phase', T.Properties.VariableNames)
        if ~ismember('PrePost', T.Properties.VariableNames)
            T = renamevars(T, 'Phase', 'PrePost');
        end
        if ismember('Condition', T.Properties.VariableNames) && ...
                ~ismember('EyeState', T.Properties.VariableNames)
            T = renamevars(T, 'Condition', 'EyeState');
        end
    end

    expected = fieldnames(defaults);
    extra    = setdiff(T.Properties.VariableNames, expected, 'stable');
    if ~isempty(extra)
        warning('APOP:LogSchema', ...
            ['%s has columns that are not part of the current schema (%s). ' ...
             'They are dropped and the file is rewritten without them -- ' ...
             'keep a copy if you still need them.'], LogFile, strjoin(extra, ', '));
    end

    missing = setdiff(expected, T.Properties.VariableNames, 'stable');
    for k = 1:numel(missing)
        name = missing{k};
        if isstring(defaults.(name))
            T.(name) = repmat(defaults.(name), height(T), 1);
        else
            T.(name) = nan(height(T), 1);
        end
    end

    T = T(:, expected);

    for k = 1:numel(expected)
        name = expected{k};
        if isstring(defaults.(name))
            T.(name) = toString(T.(name));
        else
            T.(name) = toNumeric(T.(name));
        end
    end
end
