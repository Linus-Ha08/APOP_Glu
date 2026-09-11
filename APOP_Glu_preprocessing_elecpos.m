% Final APOP_Glu preprocessing script 11/09/2026
% Windows Version

% takes raw data and performs the following preprocessing steps:
%  1. Load raw BrainVision data (BIDS layout from bidsify_apop_glu.py)
%  2. Identify EOG channels by label (aux channels '31'/'32'), not by index
%  3. Resample to 1000 Hz (raw data must already be >= 1000 Hz)
%  4. Filter: 80 Hz lowpass, 0.5 Hz highpass (EOG channels included)
%  5. Split the EOG channels off. 
%  6. Line noise removal: Zapline-plus + CleanLine (50 Hz)
%  7. clean_rawdata
%  8. ICA (extended infomax) + ICLabel, remove flagged components.
%  9. EOG QC: correlations between the EOG traces and the ICs / the EEG
%     before and after component removal
% 10. Interpolate EEG channels back to the full pre-clean_rawdata set
% 11. Re-reference to average (FCz recovered)
% 12. Re-append the EOG channels to the saved dataset 
%
% Per-file design variables:
%   PrePost  -- 'pre' / 'post' (acq- entity in the filename)
%   EyeState -- 'EO' / 'EC'    (task-rest entity in the filename)
%   Condition, Drug -- looked up by Subject + Day in APOP_Glu_Conditions.csv;
%                      anything without a unique match is flagged 'undefined'
% All four are written to the log file and to EEG.etc of each saved dataset.
%
% EOG QC columns in the log:
%   EOG_Corr_Before / EOG_Corr_After -- largest |r| between any EEG channel
%       and any EOG channel, before and after IC removal. After should be
%       clearly lower than before.
%   Max_IC_EOG_Corr    -- largest |r| between ANY IC and an EOG channel
%   Max_EyeIC_EOG_Corr -- largest |r| among the ICs ICLabel called 'eye'
%       If Max_IC_EOG_Corr is high but Max_EyeIC_EOG_Corr is low or NaN,
%       ICLabel missed an ocular component -- inspect that file.
%   ICA_Rank -- rank of the data entering ICA. If it is below EEG_Channels
%       (ASR can reduce rank) the decomposition is over-parameterised and
%       you may want 'pca', ICA_Rank in pop_runica.
%
% Recordings shorter than minRecordingSeconds are skipped before any
% processing and logged with Status 'Skipped' and the measured length.
%
% The output mirrors the input BIDS layout inside outpath:
%   preprocessed_BIDS/sub-<CODE>/ses-day<N>/<pre|post>/<EO|EC>/
%       sub-<CODE>_ses-day<N>_task-rest<EO|EC>_acq-<pre|post>_eeg_preprocessed.set
%   plus the matching .fdt. The subfolder is taken from the input folder, and is
%   cross-checked against the filename entities before the file is processed.
%
% Re-running is safe. A recording is skipped when both its .set and .fdt are
% already in the output tree, or when the log already records it as 'Skipped'
% (too short -- those never produce a .set). Everything else, including files
% that errored, is retried, and its old log row is replaced, not duplicated.
%
% ELECTRODE POSITIONS (this is what differs from APOP_Glu_preprocessing.m)
% Instead of the idealised standard_1005 template coordinates, every channel
% gets the MEASURED position from Paolo's pp_Paolo .mat for that subject and
% day (rsopen_avgref.elec: elecpos + label + unit). The positions are written
% into EEG.chanlocs BEFORE clean_rawdata, so the real geometry drives
%   - the RANSAC channel rejection in clean_rawdata
%   - the ICLabel topography feature (its main cue for eye/muscle ICs)
%   - the spherical interpolation of rejected channels
%   - the recovered FCz reference position, when the file contains FCz
% The template lookup still runs first, only to give every channel the full
% set of chanlocs fields; the measured coordinates then overwrite it.
%
% The .mat is found by globbing pp_Paolo/Subject_<CODE>_*/<CODE>_*day<N>.mat,
% because the date in the filename is not predictable. Exactly one match is
% required. It is loaded once per subject/day and cached.
%
% Two guards, both of which fail the file rather than letting it through:
%   - orientation: EEGLAB (X nose, Y left ear) and FieldTrip (X right ear,
%     Y nose) disagree about the axes. coordSystem picks the conversion, and
%     verifyOrientation then checks it against the 10-20 naming rules (odd
%     numbers left, even right; F* anterior, P*/O* posterior). A wrong
%     convention otherwise runs to completion and merely rotates every
%     topography, which is invisible in the output.
%   - completeness: if any EEG channel has no measured position the file is
%     failed by default (requireAllElecPositions), because a montage mixing
%     measured and template coordinates is worse than either alone.
%
% returns continuous, cleaned data

%% Initialize
clear; clc;
eeglabRoot  = 'D:\Linus\MATLAB_applications\eeglab2026.0.0';
ICLabelRoot = 'D:\Linus\MATLAB_applications\eeglab2026.0.0\plugins\ICLabel';
rawpath     = 'D:\Linus\APOP_Glu\raw_BIDS';

% Separate output tree, so this run does not overwrite or get confused with
% the template-coordinate output of APOP_Glu_preprocessing.m.
outpath     = 'D:\Linus\APOP_Glu\preprocessed_elecpos_BIDS';

% measured electrode positions
elecRoot    = 'D:\Linus\APOP_Glu\pp_Paolo';
elecVar     = 'rsopen_avgref';   % struct inside each .mat that holds .elec
coordSystem = 'ras';             % 'ras' = FieldTrip (X right, Y nose, Z up)
                                 % 'eeglab' = already X nose, Y left ear, Z up
                                 % verifyOrientation checks this for you
requireAllElecPositions = true;   % false = warn instead of failing the file

% Recordings shorter than this are skipped, not preprocessed. Below ~32 s
% zapline-plus errors out in findpeaks (segmentLength 1 s vs minChunkLength 30),
% and ASR/ICA have too little data to be trustworthy anyway.
minRecordingSeconds = 60;

% labels of EOG chans
eogRawLabels = {'31', '32'};
eogNewLabels = {'EOG_L', 'EOG_R'};

% Keep the EOG channels in the saved dataset (as untouched reference traces)?
% The QC columns are logged either way.
keepEOGInOutput = true;

if ~exist(outpath, 'dir')
    mkdir(outpath);
end

% logfile
LogFile = fullfile(outpath, 'APOP_Glu_preprocessing_log.csv');

% condition/drug table (lives next to this script)
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;   % running by-section: fall back to the current folder
end
CondFile = fullfile(scriptDir, 'APOP_Glu_Conditions.csv');

%% Log setup
% one struct per file, turned into a table row with struct2table.
% The text defaults below MUST be string scalars (""), not char ('').
%  - '' is 0x0 char, i.e. ZERO rows, so struct2table cannot reconcile it with
%    the 1x1 NaN fields: "fields have different numbers of rows".
%  - alignLogTable dispatches on isstring(defaults.(name)) to decide which
%    columns are text. With char defaults every text column is routed through
%    toNumeric instead, str2double("AK") gives NaN, and writetable blanks the
%    whole column -- which silently wipes Subject_ID/Raw_Filename/Status/etc
%    of every previously logged row each time the script is restarted.
% Single quotes are for PATHS only (see the note at the top), never here.
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
    'Elec_File',            "", ...
    'Elec_Matched',         NaN, ...
    'Elec_FCz',             "", ...
    'Resampled_Srate',      NaN, ...
    'Removed_Points',       NaN, ...
    'Removed_Seconds',      NaN, ...
    'Removed_Channels_N',   NaN, ...
    'Removed_Channels',     "", ...
    'ICA_Channels',         NaN, ...
    'ICA_Rank',             NaN, ...
    'Total_Removed_ICs',    NaN, ...
    'Eye_Rejected_ICs',     NaN, ...
    'Muscle_Rejected_ICs',  NaN, ...
    'Heart_Rejected_ICs',   NaN, ...
    'Line_Rejected_ICs',    NaN, ...
    'Channel_Rejected_ICs', NaN, ...
    'EOG_Corr_Before',      NaN, ...
    'EOG_Corr_After',       NaN, ...
    'Max_IC_EOG_Corr',      NaN, ...
    'Max_EyeIC_EOG_Corr',   NaN);

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

iSubj = find(strcmpi(condNames, 'Subject'),   1);
iDay  = find(strcmpi(condNames, 'Day'),       1);
iCond = find(strcmpi(condNames, 'Condition'), 1);
iDrug = find(strcmpi(condNames, 'Drug'),      1);
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

% 'eeglab nogui' adds only the TOP-LEVEL plugin folders. The eegplugin_*
% functions -- the ones that add a plugin's own subfolders via genpath -- are
% called by eeglab.m only when the GUI is built (eeglab.m line ~1143,
% 'if ~nouiflag'). CleanLine is the one installed plugin that keeps runtime
% code in subfolders: hlp_varargin2struct and arg_define in
% external/bcilab_partial, the multitaper functions in
% external/chronux_2_modified. Without this line pop_cleanline dies with
% 'Undefined function 'hlp_varargin2struct' for input arguments of type 'cell''.
% genpath skips private/ folders, so nothing here shadows EEGLAB.
addpath(genpath(fullfile(eeglabRoot, 'plugins', 'Cleanline2.1')));

%% Measured electrode positions
if exist(elecRoot, 'dir') ~= 7
    error('APOP:MissingElecRoot', 'Electrode root folder not found: %s', elecRoot);
end
% one .mat serves all four recordings of a subject/day -- load it once
elecCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

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

    % Mirror the input BIDS tree inside outpath:
    %   preprocessed_BIDS/sub-<CODE>/ses-day<N>/<pre|post>/<EO|EC>/
    % The subfolder is derived from the INPUT FOLDER rather than rebuilt from the
    % filename entities, so the output stays a faithful mirror of raw_BIDS. The
    % two are cross-checked at the top of the try block below.
    relFolder    = relativeFolder(rawFolder, rawpath);
    savedSetDir  = fullfile(outpath, relFolder);
    savedSetPath = fullfile(savedSetDir, savedSetName);

    %% resume: has this recording already been preprocessed?
    % pop_saveset writes a .set (header) and a .fdt (the data). BOTH must be
    % present to count as done: checking only the .set would permanently skip a
    % file whose run was killed between the two writes, leaving a header that
    % points at data which never arrived. The files on disk are the truth here,
    % not the log -- the log can be deleted without losing finished work.
    savedFdtPath = fullfile(savedSetDir, [base_name '_preprocessed.fdt']);
    if exist(savedSetPath, 'file') == 2 && exist(savedFdtPath, 'file') == 2
        fprintf('  Already preprocessed, skipping: %s\n', ...
            fullfile(relFolder, savedSetName));
        continue;
    elseif exist(savedSetPath, 'file') == 2 || exist(savedFdtPath, 'file') == 2
        warning('APOP:IncompleteSave', ...
            ['%s has only one of its .set/.fdt pair -- a previous run was '  ...
             'interrupted mid-save. Reprocessing and overwriting it.'], base_name);
    end

    % A recording judged too short never produces a .set, so without this it
    % would be reloaded and re-logged on every restart.
    if height(logData) > 0 && any(logData.Raw_Filename == string(rawFile) & ...
                                  logData.Status == "Skipped")
        fprintf('  Previously skipped (shorter than %g s), not retrying: %s\n', ...
            minRecordingSeconds, rawFile);
        continue;
    end

    logEntry = logFieldDefaults;
    logEntry.Subject_ID   = string(tok.subject);
    logEntry.Raw_Filename = string(rawFile);
    logEntry.Day          = str2double(tok.day);
    logEntry.PrePost      = string(tok.prepost);
    logEntry.EyeState     = string(tok.eyestate);

    %% Condition / Drug for this subject and day
    % defaults are already 'undefined'; only a unique match overwrites them
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
        %% folder / filename consistency
        % The output folder comes from the input FOLDER; Subject/Day/PrePost/
        % EyeState in the log come from the FILENAME. If the two disagree the
        % saved dataset would contradict its own log row, so fail the file here
        % rather than filing it somewhere ambiguous. This is what catches a
        % recording that was moved by hand without also being renamed.
        folderTok = regexp(relFolder, ...
            ['^sub-(?<subject>[^\\/]+)[\\/]ses-day(?<day>\d+)[\\/]' ...
             '(?<prepost>pre|post)[\\/](?<eyestate>EO|EC)$'], 'names', 'once');
        if isempty(folderTok)
            error('APOP:BadFolderLayout', ...
                ['Input folder ''%s'' does not match the expected ' ...
                 'sub-<ID>/ses-day<N>/<pre|post>/<EO|EC> layout.'], relFolder);
        end
        if ~strcmp(folderTok.subject, tok.subject) || ...
                str2double(folderTok.day) ~= str2double(tok.day) || ...
                ~strcmp(folderTok.prepost, tok.prepost) || ...
                ~strcmp(folderTok.eyestate, tok.eyestate)
            error('APOP:PathEntityMismatch', ...
                ['Folder says sub-%s / day %s / %s / %s but the filename says ' ...
                 'sub-%s / day %s / %s / %s. Rename the recording in raw_data and ' ...
                 're-run bidsify_apop_glu.py before preprocessing it.'], ...
                folderTok.subject, folderTok.day, folderTok.prepost, folderTok.eyestate, ...
                tok.subject, tok.day, tok.prepost, tok.eyestate);
        end

        %% Load raw data
        EEG = pop_loadbv(rawFolder, rawFile);
        EEG.setname = char(logEntry.Subject_ID);

        logEntry.Loaded_Srate    = EEG.srate;
        logEntry.Loaded_Points   = EEG.pnts;
        logEntry.Loaded_Seconds  = EEG.xmax;
        logEntry.Loaded_Channels = EEG.nbchan;

        %% minimum recording length
        % Raised as an error with a dedicated identifier so the catch block below
        % logs it as 'Skipped' rather than 'Error' -- the Loaded_* columns are
        % already filled in, so the log still records how long the file actually is.
        if EEG.xmax < minRecordingSeconds
            error('APOP:TooShort', ...
                'Recording is %.1f s, below the %g s minimum -- skipped, not preprocessed.', ...
                EEG.xmax, minRecordingSeconds);
        end

        %% chanlocs & EOG chans 
        % The template lookup runs only to give every channel the full set of
        % chanlocs fields. The measured coordinates overwrite X/Y/Z a few lines
        % below, before anything spatial happens.
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
        if isempty(eogChanIdx)
            error('No EOG channels found in %s.', rawFile);
        end
        logEntry.EEG_Channels = EEG.nbchan - numel(eogChanIdx);
        logEntry.EOG_Channels = numel(eogChanIdx);

        %% measured electrode positions -> EEG.chanlocs
        % Done here, after the EOG channels have been renamed (so they cannot
        % accidentally match an entry in the electrode file) and before the
        % resampling, so every later step sees the real geometry.
        elecKey = sprintf('%s_day%d', tok.subject, str2double(tok.day));
        if isKey(elecCache, elecKey)
            elecInfo = elecCache(elecKey);
        else
            elecInfo = loadElecPositions(elecRoot, tok.subject, ...
                str2double(tok.day), elecVar, coordSystem);
            elecCache(elecKey) = elecInfo;
        end
        logEntry.Elec_File = string(elecInfo.file);

        [EEG, nElecMatched] = applyElecPositions(EEG, elecInfo, requireAllElecPositions);
        EEG = eeg_checkset(EEG);
        logEntry.Elec_Matched = nElecMatched;
        fprintf('  Measured positions written to %d/%d channels from %s\n', ...
            nElecMatched, EEG.nbchan, elecInfo.file);

        %% resample 
        if EEG.srate < 1000
            error('Sampling rate is below 1000 Hz. Please check the data.');
        end
        if EEG.srate > 1000
            EEG = pop_resample(EEG, 1000);
        end
        logEntry.Resampled_Srate = EEG.srate;

        %% filter
        EEG = pop_eegfiltnew(EEG, 'hicutoff', 80);
        EEG = pop_eegfiltnew(EEG, 'locutoff', 0.5);

        %% split off the EOG channels
        % The EOG data is kept as a plain matrix and re-appended at the end.
        eogData     = EEG.data(eogChanIdx, :);
        eogChanlocs = EEG.chanlocs(eogChanIdx);
        eogLabels   = {eogChanlocs.labels};
        EEG         = pop_select(EEG, 'nochannel', eogChanIdx);
        EEG         = eeg_checkset(EEG);
        fprintf('  Split off %d EOG channel(s) (%s); %d EEG channels continue.\n', ...
            numel(eogChanIdx), strjoin(eogLabels, ', '), EEG.nbchan);

        %% line noise removal (EEG channels only)
        EEG = pop_zapline_plus(EEG, 'noisefreqs', 50, 'plotResults', 0);
        EEG = pop_cleanline(EEG, 'linefreqs', 50, 'plotfigures', 0);

        %% clean_rawdata
        % record pre cleaning
        original_chanlocs    = EEG.chanlocs;
        before_clean_points  = EEG.pnts;
        before_clean_seconds = EEG.pnts / EEG.srate;
        chans_before_clean   = {EEG.chanlocs.labels};

        % clean data
        % NOTE: settings    (1) Burst rejection off; BurstCriterion 20
        %                   (2) Burst rejection on ; BurstCriterion 30
        EEG = pop_clean_rawdata(EEG, ...
            'FlatlineCriterion', 5, ...     % reject chans if flat > 5 s
            'ChannelCriterion', 0.8, ...    % bad if chan correlated less than this to own reconstruction based on other chans; default = 0.8
            'LineNoiseCriterion', 4, ...    % more line noise relative to its signal than this value in sd
            'Highpass', 'off', ...          % highpass already applied explicitly above
            'BurstCriterion', 30, ...       % mark activity bursts that are 20 SDs above the mean (for later rejection); default = 5
            'WindowCriterion', 0.25, ...    % after ASR criterion: reject windows where more than 25% of channels show bursts
            'BurstRejection', 'on', ...    % on = reject marked bursts, off = ASR will attempt to clean them
            'Distance', 'Euclidian', ...
            'WindowCriterionTolerances', [-Inf 7]);   % numeric -- see note above
        asrLastErr = lasterr;   % capture before eeg_checkset can overwrite it
        EEG = eeg_checkset(EEG);

        %% ASR silent-failure guard
        % clean_artifacts wraps clean_asr in 'try ... catch, lasterr; return; end'
        % on its 'euclidian' branch (clean_artifacts lines 266-276), which is the
        % branch this script's 'Distance','Euclidian' selects. If ASR throws, the
        % error is printed to the command window and clean_artifacts RETURNS EARLY
        % with EEG unchanged since before ASR -- bad channels dropped, but no burst
        % correction and no window rejection -- and raises nothing. The file would
        % otherwise be logged as Success, indistinguishable from a properly cleaned
        % one, which is exactly how differential cleaning gets into a dataset.
        %
        % Detection: line 274 is the ONLY early return in clean_artifacts.
        % pop_clean_rawdata rmfields any pre-existing EEG.etc.clean_sample_mask
        % before calling it, and clean_windows -- which always runs here, since
        % WindowCriterion and WindowCriterionTolerances are both numeric -- always
        % writes the field back (clean_windows lines 156-165). Nothing after that
        % point removes it. So a missing field means, and only means, that ASR
        % failed. Note the script's EOG alignment check below does NOT catch this:
        % with no samples removed, size(eogData,2) == EEG.pnts still holds.
        %
        % Do not 'fix' this by switching Distance to 'riemannian' to get a real
        % error -- that needs clean_rawdata/manopt, which 'eeglab nogui' does not
        % put on the path either.
        if ~isfield(EEG.etc, 'clean_sample_mask') || isempty(EEG.etc.clean_sample_mask)
            error('APOP:ASRFailed', ...
                ['clean_rawdata returned no EEG.etc.clean_sample_mask: clean_asr ' ...
                 'threw and clean_artifacts returned early, so this file received ' ...
                 'neither burst correction nor window rejection. Last MATLAB error ' ...
                 'recorded at that point: %s'], asrLastErr);
        end

        % keep the EOG traces sample-aligned with the EEG.
        if isfield(EEG.etc, 'clean_sample_mask') && ~isempty(EEG.etc.clean_sample_mask)
            sampleMask = logical(EEG.etc.clean_sample_mask(:))';
            if numel(sampleMask) ~= size(eogData, 2)
                error(['clean_sample_mask has %d entries but the EOG data has %d samples ' ...
                       '-- cannot realign the EOG channels.'], ...
                       numel(sampleMask), size(eogData, 2));
            end
            eogData = eogData(:, sampleMask);
        end
        if size(eogData, 2) ~= EEG.pnts
            error('EOG data has %d samples but the EEG has %d after cleaning -- alignment lost.', ...
                size(eogData, 2), EEG.pnts);
        end

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
        % ASR can leave the data rank deficient; logged so it is visible.
        logEntry.ICA_Channels = EEG.nbchan;
        logEntry.ICA_Rank     = rank(double(EEG.data));
        if logEntry.ICA_Rank < EEG.nbchan
            warning('APOP:RankDeficient', ...
                ['%s: rank %d < %d channels entering ICA. Consider ' ...
                 '''pca'', %d in pop_runica.'], ...
                rawFile, logEntry.ICA_Rank, EEG.nbchan, logEntry.ICA_Rank);
        end

        rng(42, 'twister');   % deterministic ICA -- reset before every decomposition
        EEG = pop_runica(EEG, 'extended', 1);
        EEG = pop_iclabel(EEG, 'default');

        % ICLabel class order: 1 brain, 2 muscle, 3 eye, 4 heart,
        %                      5 line noise, 6 channel noise, 7 other
        % reject if any prob > 0.6
        ic_probs = EEG.etc.ic_classification.ICLabel.classifications;

        EEG = pop_icflag(EEG, [NaN NaN; ...
                    0.60 1;                 % muscle
                    0.60 1;                 % eye
                    0.60 1;                 % heart
                    0.60 1;                 % line noise
                    0.60 1;                 % channel noise
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

        % EOG QC (before component removal)
        % ICLabel does not look at the EOG channels, so they give an
        % independent read on whether the ocular components were caught.
        icaact  = EEG.icaweights * EEG.icasphere * double(EEG.data(EEG.icachansind, :));
        R_ic    = corrMatrix(icaact, eogData);          % nIC x nEOG
        maxPerIC = max(abs(R_ic), [], 2);
        clear icaact R_ic;

        logEntry.Max_IC_EOG_Corr = max(maxPerIC);
        if ~isempty(eye_rejICs)
            logEntry.Max_EyeIC_EOG_Corr = max(maxPerIC(eye_rejICs));
        end
        logEntry.EOG_Corr_Before = max(abs(corrMatrix(EEG.data, eogData)), [], 'all');

        % reject ICs
        if ~isempty(rejICs)
            EEG = pop_subcomp(EEG, rejICs, 0);
        end

        % EOG QC (after component removal)
        logEntry.EOG_Corr_After = max(abs(corrMatrix(EEG.data, eogData)), [], 'all');
        fprintf(['  EOG QC: EEG-EOG |r| %.2f -> %.2f | strongest IC-EOG |r| %.2f ' ...
                 '(eye ICs: %.2f)\n'], ...
            logEntry.EOG_Corr_Before, logEntry.EOG_Corr_After, ...
            logEntry.Max_IC_EOG_Corr, logEntry.Max_EyeIC_EOG_Corr);
        if logEntry.Max_IC_EOG_Corr >= 0.5 && ...
                (isnan(logEntry.Max_EyeIC_EOG_Corr) || logEntry.Max_EyeIC_EOG_Corr < 0.5)
            warning('APOP:MissedEyeIC', ...
                ['%s: an IC correlates |r| = %.2f with the EOG but ICLabel flagged ' ...
                 'no eye component that strong -- inspect this file.'], ...
                rawFile, logEntry.Max_IC_EOG_Corr);
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
        % EEG channels only
        EEG = pop_interp(EEG, original_chanlocs, 'spherical');
        EEG = eeg_checkset(EEG);

        %% re-reference to avg. and recover active reference (FCz)
        % Only EEG channels
        % The template entry is used only for its field structure; the
        % coordinates are replaced by the measured FCz when the electrode
        % file has one, so the recovered reference lands in the same montage
        % as every other channel.
        tmpl = readlocs(fullfile(eeglabRoot, 'plugins', 'dipfit', 'standard_BEM', 'elec', 'standard_1005.elc'));
        fcz_loc = tmpl(find(strcmpi({tmpl.labels}, 'FCz'), 1));
        if ~isfield(fcz_loc, 'type'), fcz_loc.type = ''; end
        if ~isfield(fcz_loc, 'ref'), fcz_loc.ref = ''; end
        if ~isfield(fcz_loc, 'urchan'), fcz_loc.urchan = []; end

        mFcz = find(strcmpi(elecInfo.labels, 'FCz'), 1);
        if ~isempty(mFcz)
            fcz_loc.X = elecInfo.xyz(mFcz, 1);
            fcz_loc.Y = elecInfo.xyz(mFcz, 2);
            fcz_loc.Z = elecInfo.xyz(mFcz, 3);
            fcz_loc   = recomputeAngles(fcz_loc);
            logEntry.Elec_FCz = "measured";
        else
            warning('APOP:TemplateFCz', ...
                ['%s: the electrode file has no FCz, so the recovered reference ' ...
                 'channel keeps its template position while every other channel ' ...
                 'is measured.'], rawFile);
            logEntry.Elec_FCz = "template";
        end

        EEG = pop_reref(EEG, [], 'refloc', fcz_loc);

        %% re-append the EOG channels
        % Filtered and cut to the same samples
        if keepEOGInOutput
            EEG = appendChannels(EEG, eogData, eogChanlocs);
            EEG = eeg_checkset(EEG);
        end

        %% save metadata 
        EEG.etc.subject_id          = char(logEntry.Subject_ID);
        EEG.etc.raw_filename        = char(logEntry.Raw_Filename);
        EEG.etc.day                 = logEntry.Day;
        EEG.etc.prepost             = char(logEntry.PrePost);    % 'pre' / 'post'
        EEG.etc.eye_state           = char(logEntry.EyeState);   % 'EO' / 'EC'
        EEG.etc.condition           = char(logEntry.Condition);  % drug condition 1-4, or 'undefined'
        EEG.etc.drug                = char(logEntry.Drug);       % drug name, or 'undefined'
        EEG.etc.original_chanlocs   = original_chanlocs;         % EEG channels before clean_rawdata
        EEG.etc.eog_channels        = eogLabels;
        EEG.etc.eog_processing      = ['filtered 0.5-80 Hz and cut to the retained ' ...
                                       'samples only; excluded from zapline, clean_rawdata/ASR, ' ...
                                       'ICA, interpolation and the average reference'];
        EEG.etc.chanlocs_source     = sprintf( ...
            ['measured elecpos from %s (%s.elec); coordSystem=%s; unit=%s; ' ...
             '%d channels matched; FCz=%s'], ...
            elecInfo.file, elecVar, coordSystem, elecInfo.unit, ...
            logEntry.Elec_Matched, logEntry.Elec_FCz);
        EEG.etc.removed_points      = logEntry.Removed_Points;
        EEG.etc.removed_seconds     = logEntry.Removed_Seconds;
        EEG.etc.removed_channels_n  = logEntry.Removed_Channels_N;
        EEG.etc.removed_channels    = char(logEntry.Removed_Channels);
        EEG.etc.ica_channels        = logEntry.ICA_Channels;
        EEG.etc.ica_rank            = logEntry.ICA_Rank;
        EEG.etc.rejected_ICs        = rejICs;
        EEG.etc.rejected_ICs_counts = struct( ...
            'total',   logEntry.Total_Removed_ICs, ...
            'eye',     logEntry.Eye_Rejected_ICs, ...
            'muscle',  logEntry.Muscle_Rejected_ICs, ...
            'heart',   logEntry.Heart_Rejected_ICs, ...
            'line',    logEntry.Line_Rejected_ICs, ...
            'channel', logEntry.Channel_Rejected_ICs);
        EEG.etc.eog_qc              = struct( ...
            'eeg_eog_corr_before', logEntry.EOG_Corr_Before, ...
            'eeg_eog_corr_after',  logEntry.EOG_Corr_After, ...
            'max_ic_eog_corr',     logEntry.Max_IC_EOG_Corr, ...
            'max_eye_ic_eog_corr', logEntry.Max_EyeIC_EOG_Corr);

        %% save preprocessed dataset
        if exist(savedSetDir, 'dir') ~= 7
            mkdir(savedSetDir);
        end
        EEG = pop_saveset(EEG, 'filename', savedSetName, 'filepath', savedSetDir);

        logEntry.Status        = string('Success');
        logEntry.Error_Message = string('');
        fprintf('Finished processing file %d/%d: %s\n', iFile, numel(files), rawFile);
        fprintf('  Saved to %s\n', fullfile(relFolder, savedSetName));

    catch ME
        if strcmp(ME.identifier, 'APOP:TooShort')
            logEntry.Status        = string('Skipped');
            logEntry.Error_Message = string(ME.message);
            fprintf(2, 'SKIPPED file %s\n%s\n', rawFile, ME.message);
        else
            logEntry.Status        = string('Error');
            logEntry.Error_Message = string(ME.message);
            fprintf(2, 'ERROR in file %s\n%s\n', rawFile, ME.message);
        end
    end

    %% update logfile 
    % Replace any earlier row for this recording instead of appending a second
    % one, so a retried file (previous run errored, or was interrupted) does not
    % leave a stale row behind. Raw_Filename is unique across the whole tree,
    % so it is a safe key.
    if height(logData) > 0
        logData(logData.Raw_Filename == string(rawFile), :) = [];
    end
    logData = [logData; struct2table(logEntry)];
    writetable(logData, LogFile);
    fprintf('Log file updated: %s\n', LogFile);
end

disp('Finished preprocessing.');

%% Local helper functions
function elecInfo = loadElecPositions(elecRoot, subject, day, elecVar, coordSystem)
% Find and read the measured electrode positions for one subject and day.
% The date in the filename is not predictable, so the file is globbed as
%   <elecRoot>/Subject_<CODE>_*/<CODE>_*day<N>.mat
% and exactly one match is required -- zero or several is a data-organisation
% problem that should stop the file, not be guessed at.
    pat  = fullfile(elecRoot, ['Subject_' subject '_*'], ...
                    [subject '_*day' num2str(day) '.mat']);
    hits = dir(pat);
    hits = hits(~[hits.isdir]);
    if isempty(hits)
        error('APOP:NoElecFile', 'No electrode file matches %s', pat);
    elseif numel(hits) > 1
        error('APOP:AmbiguousElecFile', ...
            '%d electrode files match %s (%s) -- cannot choose between them.', ...
            numel(hits), pat, strjoin({hits.name}, ', '));
    end
    f = fullfile(hits(1).folder, hits(1).name);

    S = load(f, elecVar);
    if ~isfield(S, elecVar)
        error('APOP:MissingElecVar', '%s contains no variable ''%s''.', f, elecVar);
    end
    if ~isfield(S.(elecVar), 'elec')
        error('APOP:MissingElec', '%s.%s has no .elec field (found: %s).', ...
            f, elecVar, strjoin(fieldnames(S.(elecVar))', ', '));
    end
    srcElec = S.(elecVar).elec;
    if ~isfield(srcElec, 'elecpos') || ~isfield(srcElec, 'label')
        error('APOP:MissingPosOrLabel', ...
            '%s.elec needs both elecpos and label (found: %s).', ...
            elecVar, strjoin(fieldnames(srcElec)', ', '));
    end

    pos    = double(srcElec.elecpos);
    labels = cellstr(string(srcElec.label(:)));
    if size(pos, 1) ~= numel(labels)
        error('APOP:ElecCountMismatch', '%d positions but %d labels in %s.', ...
            size(pos, 1), numel(labels), f);
    end

    % scale to mm
    unitStr = 'unspecified';
    if isfield(srcElec, 'unit') && ~isempty(srcElec.unit)
        unitStr = char(string(srcElec.unit));
    end
    switch lower(unitStr)
        case 'mm'   % already
        case 'cm',  pos = pos * 10;
        case 'dm',  pos = pos * 100;
        case 'm',   pos = pos * 1000;
        otherwise
            warning('APOP:UnknownUnit', ...
                'Unit ''%s'' in %s not recognised -- positions used as they are.', ...
                unitStr, f);
    end

    % axis convention -> EEGLAB (X nose, Y left ear, Z vertex)
    switch lower(coordSystem)
        case 'ras',    xyz = [ pos(:,2), -pos(:,1), pos(:,3) ];
        case 'eeglab', xyz = pos;
        otherwise
            error('APOP:BadCoordSystem', 'coordSystem must be ''ras'' or ''eeglab''.');
    end

    verifyOrientation(labels, xyz, coordSystem, f);

    elecInfo = struct('labels', {labels}, 'xyz', xyz, 'file', f, 'unit', unitStr);
    fprintf('  Loaded %d measured positions from %s (unit ''%s'')\n', ...
        numel(labels), f, unitStr);
end

function verifyOrientation(labels, xyz, coordSystem, elecFile)
% Catch a wrong axis convention automatically instead of hoping someone reads
% a printout. A flipped convention does not error anywhere downstream -- it
% just rotates every topography, which is invisible in the saved data.
%
% Uses the 10-20 naming rules, so it needs no hardcoded electrode list:
%   odd trailing number = LEFT hemisphere, even = RIGHT
%   Fp*/AF*/F* are anterior, P*/PO*/O* posterior
% In EEGLAB coordinates X points at the nose and Y at the LEFT ear, so the
% left group must have the larger mean Y and the anterior group the larger
% mean X.
    n       = numel(labels);
    isLeft  = false(n, 1);
    isRight = false(n, 1);
    for i = 1:n
        t = regexp(labels{i}, '(\d+)$', 'tokens', 'once');
        if isempty(t), continue; end
        if mod(str2double(t{1}), 2) == 1
            isLeft(i) = true;
        else
            isRight(i) = true;
        end
    end
    isAnt  = ~cellfun(@isempty, regexpi(labels, '^(Fp|AF|F)\d', 'once'));
    isPost = ~cellfun(@isempty, regexpi(labels, '^(PO|P|O)\d', 'once'));
    isAnt  = isAnt(:);
    isPost = isPost(:);

    % A magnitude threshold, not just a sign test. On a left/right symmetric
    % montage a wrong convention gives a separation of almost exactly zero, so
    % which side of a sign test it lands on is decided by rounding noise -- a
    % coin flip. A correct convention separates the groups by most of the head.
    scale  = max(max(xyz, [], 1) - min(xyz, [], 1));
    minSep = 0.25 * scale;
    hint   = sprintf(['Set coordSystem to the other value (currently ''%s'') and ' ...
                      're-run. Electrode file: %s'], coordSystem, elecFile);

    if nnz(isLeft) >= 2 && nnz(isRight) >= 2
        sepLR = mean(xyz(isLeft, 2)) - mean(xyz(isRight, 2));
        if sepLR <= minSep
            error('APOP:BadOrientation', ...
                ['Left/right axis is wrong: odd-numbered (left) electrodes sit only ' ...
                 '%.1f mm further along EEGLAB-Y than even-numbered (right) ones ' ...
                 '(need > %.1f mm on a %.0f mm head). %s'], sepLR, minSep, scale, hint);
        end
    end
    if nnz(isAnt) >= 2 && nnz(isPost) >= 2
        sepAP = mean(xyz(isAnt, 1)) - mean(xyz(isPost, 1));
        if sepAP <= minSep
            error('APOP:BadOrientation', ...
                ['Front/back axis is wrong: frontal electrodes sit only %.1f mm ' ...
                 'further along EEGLAB-X than posterior ones (need > %.1f mm on a ' ...
                 '%.0f mm head). %s'], sepAP, minSep, scale, hint);
        end
    end

    [~, iFront] = max(xyz(:,1));
    [~, iLeft]  = max(xyz(:,2));
    [~, iTop]   = max(xyz(:,3));
    fprintf(['  Orientation ok (coordSystem ''%s''): most anterior %s, ' ...
             'most leftward %s, most superior %s\n'], ...
        coordSystem, labels{iFront}, labels{iLeft}, labels{iTop});
end

function [EEG, nMatched] = applyElecPositions(EEG, elecInfo, requireAll)
% Overwrite X/Y/Z in EEG.chanlocs with the measured positions, matching by
% label, then refresh the derived angles.
    isEog   = strcmpi({EEG.chanlocs.type}, 'EOG');
    matched = false(1, EEG.nbchan);
    for c = 1:EEG.nbchan
        m = find(strcmpi(elecInfo.labels, EEG.chanlocs(c).labels), 1);
        if isempty(m), continue; end
        EEG.chanlocs(c).X = elecInfo.xyz(m, 1);
        EEG.chanlocs(c).Y = elecInfo.xyz(m, 2);
        EEG.chanlocs(c).Z = elecInfo.xyz(m, 3);
        matched(c) = true;
    end
    nMatched = nnz(matched);
    if nMatched == 0
        error('APOP:NoElecMatch', ...
            'No label in %s matches any channel in this recording.', elecInfo.file);
    end

    % EOG channels are expected to be unmatched and must stay without coordinates
    missedEeg = {EEG.chanlocs(~matched & ~isEog).labels};
    if ~isempty(missedEeg)
        msg = sprintf(['%d EEG channel(s) have no measured position (%s) in %s. ' ...
            'Mixing measured and template coordinates in one montage corrupts the ' ...
            'RANSAC channel rejection, the ICLabel topography and the ' ...
            'interpolation.'], numel(missedEeg), strjoin(missedEeg, ', '), elecInfo.file);
        if requireAll
            error('APOP:IncompleteElec', '%s', msg);
        else
            warning('APOP:IncompleteElec', ...
                '%s Keeping the template coordinates for those channels.', msg);
        end
    end

    EEG.chanlocs = recomputeAngles(EEG.chanlocs);
end

function chanlocs = recomputeAngles(chanlocs)
% Refresh theta/radius/sph_* from X/Y/Z. topoplot, pop_interp and ICLabel read
% the angles, not the cartesian coordinates, so overwriting X/Y/Z alone would
% change nothing visible.
%
% Only channels that actually have coordinates are handed to convertlocs: it
% builds [chanlocs.X] internally, and an empty entry (the EOG channels) would
% silently shift every later channel's coordinates by one.
    hasXYZ = find(~cellfun(@isempty, {chanlocs.X}));
    if isempty(hasXYZ), return; end
    conv = convertlocs(chanlocs(hasXYZ), 'cart2all');
    for f = {'theta', 'radius', 'sph_theta', 'sph_phi', 'sph_radius'}
        if ~isfield(conv, f{1}), continue; end
        for j = 1:numel(hasXYZ)
            chanlocs(hasXYZ(j)).(f{1}) = conv(j).(f{1});
        end
    end
end

function rel = relativeFolder(folder, root)
% Path of FOLDER relative to ROOT, e.g. 'sub-AK\ses-day1\pre\EO'.
% Used to mirror the input BIDS tree in the output folder. Case-insensitive
% because dir() may not return the root spelled exactly as it was given.
    folder = regexprep(char(folder), '[\\/]+$', '');   % drop trailing separators
    root   = regexprep(char(root),   '[\\/]+$', '');
    if numel(folder) < numel(root) || ~strncmpi(folder, root, numel(root))
        error('APOP:PathOutsideRoot', ...
            'Folder ''%s'' is not below the raw data root ''%s''.', folder, root);
    end
    rel = regexprep(folder(numel(root)+1:end), '^[\\/]+', '');
end

function R = corrMatrix(A, B)
% Pearson correlation between every row of A and every row of B.
% A: nA x time, B: nB x time -> R: nA x nB. Toolbox-free.
    A = double(A);
    B = double(B);
    A = A - mean(A, 2);
    B = B - mean(B, 2);
    A = A ./ (sqrt(sum(A.^2, 2)) + eps);
    B = B ./ (sqrt(sum(B.^2, 2)) + eps);
    R = A * B';
end

function EEG = appendChannels(EEG, newData, newChanlocs)
% Append extra channels (data + chanlocs) to the end of an EEG struct,
% harmonising the chanlocs fields first so the structs can be concatenated.
    if size(newData, 2) ~= size(EEG.data, 2)
        error('appendChannels: %d samples to append but the dataset has %d.', ...
            size(newData, 2), size(EEG.data, 2));
    end

    fEEG = fieldnames(EEG.chanlocs);
    fNew = fieldnames(newChanlocs);
    for k = 1:numel(fEEG)
        if ~isfield(newChanlocs, fEEG{k})
            [newChanlocs.(fEEG{k})] = deal([]);
        end
    end
    for k = 1:numel(fNew)
        if ~isfield(EEG.chanlocs, fNew{k})
            [EEG.chanlocs.(fNew{k})] = deal([]);
        end
    end
    newChanlocs = orderfields(newChanlocs, fieldnames(EEG.chanlocs));

    EEG.data     = [EEG.data; cast(newData, 'like', EEG.data)];
    EEG.chanlocs = [EEG.chanlocs(:)', newChanlocs(:)'];
    EEG.nbchan   = size(EEG.data, 1);
end

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
