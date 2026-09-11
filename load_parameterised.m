function [records, n_nodrug] = load_parameterised(inPath, studyLabel)
%LOAD_PARAMETERISED  Read the APOP *_parameterised.mat tree into a record array.
%
%   records = LOAD_PARAMETERISED(inPath) walks inPath RECURSIVELY for
%   *_parameterised.mat files (as written by APOP_Glu_Specparam.ipynb, which
%   nests them as sub-*/ses-day*/{pre,post}/{EO,EC}/) and returns a struct
%   array carrying the fields the CBPT_* scripts expect:
%
%       UniqueID    subject identifier (string)
%       Drug        drug name, lower-case (string)
%       PrePost     "pre" | "post"
%       Eyes        "open" | "closed"
%       Exponent    1 x nChan aperiodic exponent per channel
%       Chanlabels  1 x nChan cell of channel labels
%
%   records = LOAD_PARAMETERISED(inPath, studyLabel) prefixes UniqueID with
%   studyLabel, so subject codes from two studies cannot collide when pooled.
%   Omit it for a single-study analysis.
%
%   [records, n_nodrug] = LOAD_PARAMETERISED(...) also returns how many files
%   were dropped for having no drug assignment.
%
%   NAME MAPPING (the .mat files are written by scipy.io.savemat from Python,
%   so every variable is lower-case and MATLAB field access is case-sensitive):
%       eyestate "EO"/"EC"  ->  Eyes "open"/"closed"
%       subject             ->  UniqueID   (there is no Study field)
%       drug, prepost       ->  Drug, PrePost
%       exponent_ch         ->  Exponent
%       chanlabels          ->  Chanlabels
%
%   Sessions whose drug is "undefined" (the notebook's value when a
%   (subject, day) pair is missing from APOP_Glu_Conditions.csv) or "unknown"
%   are dropped, so they cannot enter an analysis as a drug arm of their own.
%
%   Channel labels are checked against the first file. A file whose labels
%   differ in CONTENT OR ORDER raises an error: the CBPT adjacency graph is
%   built once from the first file's montage, so a reordered file would be
%   silently misaligned onto the wrong electrodes.

    if nargin < 2 || strlength(string(studyLabel)) == 0
        prefix = "";
    else
        prefix = string(studyLabel) + "_";
    end

    files = dir(fullfile(inPath, '**', '*_parameterised.mat'));
    if isempty(files)
        error('load_parameterised:noFiles', ...
              'No *_parameterised.mat files found under %s', inPath);
    end

    records = struct('Eyes', {}, 'UniqueID', {}, 'Drug', {}, 'PrePost', {}, ...
                     'Exponent', {}, 'Chanlabels', {});
    ref_labels = {};
    n_nodrug   = 0;

    for f = files'
        data = load(fullfile(f.folder, f.name));

        % --- channel labels must match the first file exactly, order included ---
        labels = cellstr(data.chanlabels);
        labels = labels(:)';
        if isempty(ref_labels)
            ref_labels = labels;
        elseif ~isequal(labels, ref_labels)
            error('load_parameterised:chanMismatch', ...
                 ['Channel labels differ (content or order) in %s.\n' ...
                  'The adjacency graph assumes one common montage across files, ' ...
                  'so this must be resolved rather than ignored.'], f.name);
        end

        % --- drop sessions with no drug assignment ---
        drug = lower(strtrim(string(data.drug)));
        if ismember(drug, ["undefined", "unknown", ""])
            n_nodrug = n_nodrug + 1;
            continue
        end

        % --- eyestate is stored EO/EC; the scripts speak open/closed ---
        switch upper(strtrim(string(data.eyestate)))
            case "EO", eyes = "open";
            case "EC", eyes = "closed";
            otherwise
                error('load_parameterised:eyestate', ...
                      'Unrecognised eyestate "%s" in %s', string(data.eyestate), f.name);
        end

        records(end+1).Eyes     = eyes;                       %#ok<AGROW>
        records(end).UniqueID   = prefix + strtrim(string(data.subject));
        records(end).Drug       = drug;
        records(end).PrePost    = lower(strtrim(string(data.prepost)));
        records(end).Exponent   = double(data.exponent_ch(:))';
        records(end).Chanlabels = ref_labels;
    end

    if isempty(records)
        error('load_parameterised:noUsableFiles', ...
              ['%d file(s) found under %s but none had a drug assignment. ' ...
               'Check that APOP_Glu_Conditions.csv covers every (Subject, Day).'], ...
              numel(files), inPath);
    end
end
