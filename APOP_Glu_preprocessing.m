% Final APOP_Glu preprocessing script 09/09/2026
% Windows Version

% takes raw data and performs the following preprocessing steps:
% 1. Load raw BrainVision data (BIDS layout from bidsify_apop_glu.py)
% 2. Identify EOG channels
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
% returns continuous, cleaned data

%% Initialize
clear; clc;
eeglabRoot  = "D:\Linus\MATLAB_applications\eeglab2026.0.0";
ICLabelRoot = "D:\Linus\MATLAB_applications\eeglab2026.0.0\plugins\ICLabel";
rawpath     = "D:\Linus\APOP_Glu\raw_BIDS";
outpath     = "D:\Linus\Local\APOP\PreprocessedData";

if ~exist(outpath, 'dir')
    mkdir(outpath);
end

% logfile
LogFile = fullfile(outpath, 'APOP_Glu_preprocessing_log.csv');

%% Log setup -- one struct per file, turned into a table row with struct2table.
% Text fields use string(...) rather than plain char so that rows with
% different string lengths (e.g. different subject codes, error messages)
% concatenate cleanly into the growing table.
logFieldDefaults = struct( ...
    'Subject_ID',           "", ...
    'Raw_Filename',         "", ...
    'Day',                  NaN, ...
    'Phase',                "", ...
    'Condition',            "", ...
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
    logData = readtable(LogFile);
else
    logData = struct2table(logFieldDefaults);
    logData(1,:) = [];
end

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
        'task-rest(?<cond>EO|EC)_acq-(?<phase>pre|post)_eeg\.vhdr$'], 'names', 'once');
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
    logEntry.Phase        = string(tok.phase);
    logEntry.Condition    = string(tok.cond);

    try
        %% ---------------- LOAD ----------------
        EEG = pop_loadbv(rawFolder, rawFile);
        EEG.setname = char(logEntry.Subject_ID);

        logEntry.Loaded_Srate    = EEG.srate;
        logEntry.Loaded_Points   = EEG.pnts;
        logEntry.Loaded_Seconds  = EEG.xmax;
        logEntry.Loaded_Channels = EEG.nbchan;

        %% ---------------- CHANNEL LOCATIONS ----------------
        EEG = pop_chanedit(EEG, 'lookup', fullfile(eeglabRoot, 'plugins', 'dipfit', ...
            'standard_BEM', 'elec', 'standard_1005.elc'));

        % Positions 31/32 are the aux EOG channels in this montage
        % (verify this holds for every subject before trusting it).
        EEG.chanlocs(31).labels = 'EOG_L';
        EEG.chanlocs(31).type   = 'EOG';
        EEG.chanlocs(32).labels = 'EOG_R';
        EEG.chanlocs(32).type   = 'EOG';
        EEG = eeg_checkset(EEG);

        eogChanIdx = find(strcmpi({EEG.chanlocs.type}, 'EOG'));
        logEntry.EEG_Channels = EEG.nbchan - numel(eogChanIdx);
        logEntry.EOG_Channels = numel(eogChanIdx);

        %% ---------------- RESAMPLE ----------------
        if EEG.srate < 1000
            error('Sampling rate is below 1000 Hz. Please check the data.');
        end
        if EEG.srate > 1000
            EEG = pop_resample(EEG, 1000);
        end
        logEntry.Resampled_Srate = EEG.srate;

        %% ---------------- FILTER ----------------
        % lowpass at 80 Hz, highpass at 0.5 Hz (clean_rawdata's own highpass
        % is disabled below, so the data is only highpassed once, here)
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 80);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);

        % line noise removal
        EEG = pop_zapline_plus(EEG, 'noisefreqs', 50, 'plotResults', 0);
        EEG = pop_cleanline(EEG, 'linefreqs', 50, 'plotfigures', 0);

        %% ---------------- BEFORE CLEAN_RAWDATA ----------------
        original_chanlocs    = EEG.chanlocs;
        before_clean_points  = EEG.pnts;
        before_clean_seconds = EEG.pnts / EEG.srate;
        chans_before_clean   = {EEG.chanlocs.labels};

        %% ---------------- CLEAN_RAWDATA ----------------
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

        %% ---------------- AFTER CLEAN_RAWDATA ----------------
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

        %% ---------------- ICA + ICLabel ----------------
        EEG = pop_runica(EEG, 'extended', 1);
        EEG = pop_iclabel(EEG, 'default');

        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;

        % mark if:  eye > 50%,
        %           muscle > 50%,
        %           heart > 50%,
        %           line noise > 80%,
        %           channel noise > 50%
        EEG = pop_icflag(EEG, [NaN NaN; ...
                    0.50 1;                 % eye
                    0.50 1;                 % muscle
                    0.50 1;                 % heart
                    0.80 1;                 % line noise
                    0.50 1;                 % channel noise
                    NaN NaN]);

        rej_mask = EEG.reject.gcompreject;
        rejICs   = find(rej_mask);

        eye_rejICs     = find(rej_mask & ic_probs(:,3) >= 0.50);
        muscle_rejICs  = find(rej_mask & ic_probs(:,2) >= 0.50);
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

        %% ---------------- Interpolate EEG channels ----------------
        % Only the non-EOG channels are restored to their pre-clean_rawdata
        % set -- spherical interpolation from scalp neighbours isn't valid
        % for a dropped EOG channel, so an EOG channel that clean_rawdata
        % removed is simply left removed (visible via Removed_Channels above).
        original_eeg_chanlocs = original_chanlocs(~strcmpi({original_chanlocs.type}, 'EOG'));
        EEG = pop_interp(EEG, original_eeg_chanlocs, 'spherical');
        EEG = eeg_checkset(EEG);

        %% ---------------- Re-reference to average ----------------
        % EOG channels are excluded from the reference computation and are
        % left un-re-referenced.
        eogChanIdx = find(strcmpi({EEG.chanlocs.type}, 'EOG'));

        tmpl = readlocs(fullfile(eeglabRoot, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc'));
        fcz_loc = tmpl(find(strcmpi({tmpl.labels}, 'FCz'), 1));
        if ~isfield(fcz_loc, 'type'), fcz_loc.type = ''; end
        if ~isfield(fcz_loc, 'ref'), fcz_loc.ref = ''; end
        if ~isfield(fcz_loc, 'urchan'), fcz_loc.urchan = []; end

        EEG = pop_reref(EEG, [], 'refloc', fcz_loc, 'exclude', eogChanIdx);

        %% ---------------- SAVE METADATA ----------------
        EEG.etc.subject_id          = char(logEntry.Subject_ID);
        EEG.etc.raw_filename        = char(logEntry.Raw_Filename);
        EEG.etc.day                 = logEntry.Day;
        EEG.etc.phase               = char(logEntry.Phase);
        EEG.etc.condition           = char(logEntry.Condition);
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

        %% ---------------- SAVE SET ----------------
        EEG = pop_saveset(EEG, 'filename', savedSetName, 'filepath', outpath);

        logEntry.Status        = string('Success');
        logEntry.Error_Message = string('');
        fprintf('Finished processing file %d/%d: %s\n', iFile, numel(files), rawFile);

    catch ME
        logEntry.Status        = string('Error');
        logEntry.Error_Message = string(ME.message);
        fprintf(2, 'ERROR in file %s\n%s\n', rawFile, ME.message);
    end

    %% ---------------- UPDATE LOG FILE (always -- success or error) ----------------
    logData = [logData; struct2table(logEntry)];
    writetable(logData, LogFile);
    fprintf('Log file updated: %s\n', LogFile);
end

disp('Finished preprocessing.');
