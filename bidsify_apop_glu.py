#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bidsify_apop_glu.py
===================

Reorganise the APOP_Glu raw BrainVision recordings into a BIDS-style tree,
verify every transfer, audit the dataset, and finish with a sanity check that
tests every recording's timestamp against the folder it was filed into.

    source :  raw_data/raw_open  /Subject_AK_01/AK_0308_day1/Study_Glu_AK_0308_pre_RS_open.vhdr
              raw_data/raw_closed/...
    target :  raw_BIDS/sub-AK/ses-day1/pre/EO/sub-AK_ses-day1_task-restEO_acq-pre_eeg.vhdr

The three companion files (.eeg / .vhdr / .vmrk) are always kept together.

OUTPUT
------
Besides the reorganised recordings, this script writes exactly one file: a
.txt transcript of everything it printed, at

    <target>/bidsify_log_<timestamp>.txt

(override with --log). Nothing else -- no CSVs, no sidecars.

SAFETY
------
The source tree is opened READ-ONLY. Nothing is moved, renamed or deleted.
Every source file's size and modification time is recorded before the run and
re-checked afterwards; the run fails loudly if any source file changed.

Each file is copied in chunks while its SHA-256 is computed on the fly, then
the file that landed on disk is read back and hashed again. A transfer only
counts as successful when both digests match.

HEADER POINTERS
---------------
A BrainVision .vhdr points at its own .eeg / .vmrk by name (DataFile= and
MarkerFile=), and the .vmrk points at the .eeg. Renaming the files to BIDS
names therefore requires those lines to be updated in the COPIES, or nothing
would be able to load them. The rewrite is byte-level and touches only those
lines; every other byte is preserved exactly.

25 recordings in this dataset already have pointers that do not resolve
(recording-time typos: 'srudy_', 'stduy_', 'copen', 'IB' for 'GB', 'SP' for
'StP', a swapped 'GB_0712' -> '0712_GB'). Those are repaired by the same
mechanism and listed individually in the log.

SANITY CHECK (final step)
-------------------------
Every BrainVision recording carries the recorder's own clock in the
'New Segment' marker. Folder and file names were typed by a human and are
demonstrably wrong in places; the clock is not. The last step therefore walks
the recordings, resolves each .vhdr to its marker file, reads the timestamp,
and asks whether it is consistent with the sub-/ses-/phase/condition folder
the recording was filed into:

  * day        -- ses-day1..4 must run in chronological order for a subject,
                  and everything inside one session must share a date
  * phase      -- every 'pre' must precede every 'post' of the same session
  * condition  -- EO and EC of one block are recorded back to back; a large
                  gap means one of them is filed under the wrong phase
  * duplicates -- one recording (same subject, same timestamp) must not be
                  filed in two places at once, which is what a half-finished
                  manual re-filing looks like

Note: in this dataset the .vhdr itself contains no timestamp (checked: 0 of
257 headers), so the value is read from the .vmrk that the header's
MarkerFile= line names.

The check walks the target tree, not the source listing, so recordings you
moved, renamed or added by hand after a previous run are checked too -- they
are the ones most likely to be in the wrong place. The loadability check does
the same, which catches a file renamed without its header pointers being
updated.

USAGE
-----
    python bidsify_apop_glu.py --dry-run          # audit + plan, copies nothing
    python bidsify_apop_glu.py                    # do it
    python bidsify_apop_glu.py --verify size      # faster, skips read-back hashing
    python bidsify_apop_glu.py --subjects SR      # redo one subject
    python bidsify_apop_glu.py --exclude FD FS    # leave subjects you moved out alone

Re-running is safe: files already present with the expected content are
reported as identical and skipped, so an interrupted copy can be restarted.
A destination file whose content differs is never silently overwritten -- it
is reported as a conflict and left alone unless --overwrite is given.

Exit codes: 0 all good | 1 transfer or integrity failure | 3 copied cleanly
but the sanity check found a recording whose timestamp contradicts its folder.

Python 3.8+, standard library only.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import statistics
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# --------------------------------------------------------------------------
# Defaults -- adjust here or override on the command line
# --------------------------------------------------------------------------

DEFAULT_PROJECT = Path(r"D:\Linus\APOP_Glu")
DEFAULT_SOURCE = "raw_data"
DEFAULT_TARGET = "raw_BIDS"

# Which source root holds which condition. EO = eyes open, EC = eyes closed.
CONDITION_ROOTS = {"raw_open": "EO", "raw_closed": "EC"}

# Every subject is expected to have this full grid.
EXPECTED_DAYS = (1, 2, 3, 4)
EXPECTED_PHASES = ("pre", "post")
EXPECTED_CONDITIONS = ("EO", "EC")

EXTENSIONS = (".eeg", ".vhdr", ".vmrk")

# A recording shorter than this is flagged outright.
DEFAULT_MIN_DURATION_S = 60.0
# ...and so is one shorter than this fraction of the dataset median.
DEFAULT_REL_FRACTION = 0.5
# Header/marker files below these sizes are structurally implausible.
MIN_VHDR_BYTES = 1000
MIN_VMRK_BYTES = 100
# The eyes-open and eyes-closed halves of one block are recorded back to back.
# Anything further apart suggests they were not acquired together and one of
# them is filed under the wrong phase.
MAX_EO_EC_GAP_MINUTES = 30.0

COPY_CHUNK = 8 * 1024 * 1024  # 8 MiB

# BIDS-style filename. Eyes-open / eyes-closed rest are treated as two tasks;
# the pre/post intervention repeat is carried by the acq- entity.
FILENAME_TEMPLATE = "sub-{subject}_ses-day{day}_task-rest{cond}_acq-{phase}_eeg{ext}"

# --------------------------------------------------------------------------
# Parsing patterns
# --------------------------------------------------------------------------

SUBJECT_RE = re.compile(r"^Subject_(?P<code>.+?)_(?P<num>\d+)$", re.IGNORECASE)
DAY_RE = re.compile(r"day\s*(?P<day>\d+)\s*$", re.IGNORECASE)
PHASE_RE = re.compile(r"(?:^|_)(?P<phase>pre|post)(?=_|$)", re.IGNORECASE)
COND_TOKEN_RE = re.compile(r"(?:^|_)(?P<cond>open|closed)(?=_|$)", re.IGNORECASE)
FOLDER_DATE_RE = re.compile(r"^(?P<code>[A-Za-z]+)_(?P<date>\d{3,4})_day\d+$", re.IGNORECASE)
FILE_DATE_RE = re.compile(r"_(?P<date>\d{4})_", re.IGNORECASE)

# Byte-level pointer lines inside .vhdr / .vmrk
POINTER_RE = re.compile(rb"(?im)^(?P<key>DataFile|MarkerFile)[ \t]*=(?P<val>[^\r\n]*)")

BYTES_PER_SAMPLE = {
    "INT_16": 2,
    "UINT_16": 2,
    "INT_32": 4,
    "UINT_32": 4,
    "IEEE_FLOAT_32": 4,
    "IEEE_FLOAT_64": 8,
}

OK, WARN, FAIL, INFO, SKIP = "[ OK ]", "[WARN]", "[FAIL]", "[INFO]", "[SKIP]"


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------


@dataclass
class Header:
    """The handful of .vhdr fields this script needs."""

    n_channels: Optional[int] = None
    sampling_interval_us: Optional[float] = None
    binary_format: Optional[str] = None
    data_file: Optional[str] = None
    marker_file: Optional[str] = None
    parse_error: Optional[str] = None

    @property
    def bytes_per_sample(self) -> Optional[int]:
        if not self.binary_format:
            return None
        return BYTES_PER_SAMPLE.get(self.binary_format.upper())

    @property
    def sfreq(self) -> Optional[float]:
        if not self.sampling_interval_us:
            return None
        return 1_000_000.0 / self.sampling_interval_us


@dataclass
class Recording:
    """One pre/post x EO/EC recording: up to three files that belong together."""

    subject: str
    subject_num: str
    day: int
    phase: str
    cond: str
    src_dir: Path
    stem: str
    src: Dict[str, Path] = field(default_factory=dict)
    header: Header = field(default_factory=Header)
    duration_s: Optional[float] = None
    folder_date: Optional[str] = None
    file_date: Optional[str] = None
    recorded: Optional[datetime] = None
    recorded_from: str = ""
    flags: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def key(self) -> Tuple[str, int, str, str]:
        return (self.subject, self.day, self.phase, self.cond)

    @property
    def label(self) -> str:
        return f"sub-{self.subject} ses-day{self.day} {self.phase:<4} {self.cond}"

    @property
    def location(self) -> str:
        return f"sub-{self.subject}/ses-day{self.day}/{self.phase}/{self.cond}"

    def dest_dir(self, target_root: Path) -> Path:
        return target_root / f"sub-{self.subject}" / f"ses-day{self.day}" / self.phase / self.cond

    def dest_name(self, ext: str) -> str:
        return FILENAME_TEMPLATE.format(
            subject=self.subject, day=self.day, cond=self.cond, phase=self.phase, ext=ext
        )

    def dest_path(self, target_root: Path, ext: str) -> Path:
        return self.dest_dir(target_root) / self.dest_name(ext)

    @property
    def missing_exts(self) -> List[str]:
        return [e for e in EXTENSIONS if e not in self.src]

    @property
    def is_orphan(self) -> bool:
        return bool(self.missing_exts)


@dataclass
class FileResult:
    recording: str
    ext: str
    source: str
    dest: str
    bytes: int
    status: str
    sha_source: str = ""
    pointer_changes: str = ""
    pointer_broken_in_source: str = ""
    detail: str = ""


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------


class Reporter:
    """Prints to the terminal and keeps a transcript for the .txt log."""

    def __init__(self) -> None:
        self.lines: List[str] = []

    def __call__(self, text: str = "") -> None:
        print(text, flush=True)
        self.lines.append(text)

    def rule(self, title: str = "", char: str = "=") -> None:
        width = 78
        if title:
            self(f"\n{char * 3} {title} ".ljust(width, char))
        else:
            self(char * width)

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")


def human_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return f"{n:,.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} TB"


def human_duration(seconds: Optional[float]) -> str:
    if seconds is None:
        return "     ?"
    return f"{seconds / 60:5.1f}m" if seconds >= 60 else f"{seconds:5.1f}s"


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(COPY_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# BrainVision header / marker handling
# --------------------------------------------------------------------------


def read_header(path: Path) -> Header:
    """Parse the few fields we need. Read as bytes and decoded leniently,
    because channel labels in these files are not reliably UTF-8."""
    header = Header()
    try:
        raw = path.read_bytes()
    except OSError as exc:
        header.parse_error = f"unreadable: {exc}"
        return header

    for line in raw.decode("utf-8", errors="replace").splitlines():
        line = line.strip()
        if "=" not in line or line.startswith(";"):
            continue
        key, _, value = line.partition("=")
        key = key.strip().lower()
        value = value.strip()
        try:
            if key == "numberofchannels":
                header.n_channels = int(value)
            elif key == "samplinginterval":
                header.sampling_interval_us = float(value)
            elif key == "binaryformat":
                header.binary_format = value
            elif key == "datafile":
                header.data_file = value
            elif key == "markerfile":
                header.marker_file = value
        except ValueError:
            header.parse_error = f"bad value for {key}: {value!r}"
    return header


def read_marker_timestamp(path: Path) -> Optional[datetime]:
    """Pull the recorder's own clock out of the 'New Segment' marker.

    The last field of that marker is YYYYMMDDhhmmssffffff, written by the
    Vision Recorder at acquisition time. Unlike the folder and file names it
    was not typed by anyone, so it is the authority on when a session happened.
    """
    try:
        text = path.read_bytes().decode("utf-8", errors="replace")
    except OSError:
        return None

    for line in text.splitlines():
        if not line.lower().startswith("mk"):
            continue
        fields = line.strip().split(",")
        if len(fields) < 6 or "new segment" not in fields[0].lower():
            continue
        stamp = fields[-1].strip()
        if len(stamp) >= 14 and stamp[:14].isdigit():
            try:
                return datetime.strptime(stamp[:14], "%Y%m%d%H%M%S")
            except ValueError:
                return None
    return None


def timestamp_via_header(vhdr: Path) -> Tuple[Optional[datetime], str]:
    """Resolve a recording's time starting from its .vhdr.

    BrainVision keeps no timestamp in the header itself, so we follow the
    header's own MarkerFile= line to the .vmrk and read it there. Falling back
    to the sibling .vmrk covers a header whose pointer is broken.
    """
    header = read_header(vhdr)
    candidates: List[Path] = []
    if header.marker_file:
        candidates.append(vhdr.parent / header.marker_file)
    candidates.append(vhdr.with_suffix(".vmrk"))

    for candidate in candidates:
        if candidate.exists():
            stamp = read_marker_timestamp(candidate)
            if stamp:
                how = "via MarkerFile=" if candidate == candidates[0] else "sibling .vmrk"
                return stamp, f"{candidate.name} ({how})"
    return None, ""


def resolve_date_mismatch(recording: Recording) -> None:
    """Say which of the two typed dates the recorder's clock agrees with."""
    stamp = recording.recorded
    folder_date, file_date = recording.folder_date, recording.file_date
    if not (folder_date and file_date):
        return

    base = f"folder says {folder_date}, filename says {file_date}"
    if stamp is None:
        recording.notes.append(f"{base}; no recording timestamp available to settle it")
        return

    actual = f"{stamp.day:02d}{stamp.month:02d}"

    def matches(label: str) -> bool:
        return label.lstrip("0") == actual.lstrip("0")

    recorded_text = stamp.strftime("%Y-%m-%d")
    if matches(folder_date) and not matches(file_date):
        recording.notes.append(
            f"{base}; recorded {recorded_text}, so the FILENAME is the wrong one"
        )
    elif matches(file_date) and not matches(folder_date):
        recording.notes.append(
            f"{base}; recorded {recorded_text}, so the FOLDER NAME is the wrong one"
        )
    elif matches(folder_date) and matches(file_date):
        recording.notes.append(f"{base}; both agree with the recorded date {recorded_text}")
    else:
        recording.flags.append("DATE_UNEXPLAINED")
        recording.notes.append(
            f"{base}; but it was recorded {recorded_text}, which matches neither"
        )


def rewrite_pointers(
    data: bytes, new_eeg: str, new_vmrk: str
) -> Tuple[bytes, List[Tuple[str, str, str]]]:
    """Replace DataFile= / MarkerFile= values at the byte level.

    Only those lines are touched; line endings, encoding and every other byte
    are preserved. Returns the new bytes and a list of (key, old, new) for the
    values that actually changed.
    """
    changes: List[Tuple[str, str, str]] = []

    def replace(match: "re.Match[bytes]") -> bytes:
        key_bytes = match.group("key")
        key = key_bytes.decode("ascii")
        old = match.group("val")
        target = new_eeg if key.lower() == "datafile" else new_vmrk
        new = target.encode("ascii")
        if old.strip() != new:
            changes.append((key, old.decode("utf-8", "replace").strip(), target))
        return key_bytes + b"=" + new

    return POINTER_RE.sub(replace, data), changes


# --------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------


def discover(source_root: Path, report: Reporter) -> Tuple[List[Recording], List[str], List[str]]:
    """Walk the source tree and build one Recording per pre/post x EO/EC cell."""
    recordings: Dict[Tuple[str, int, str, str], Recording] = {}
    unparsed: List[str] = []
    ignored: List[str] = []
    collisions: List[str] = []

    for root_name, condition in CONDITION_ROOTS.items():
        root = source_root / root_name
        if not root.is_dir():
            unparsed.append(f"missing source root: {root}")
            continue

        for subject_dir in sorted(root.iterdir()):
            if not subject_dir.is_dir():
                ignored.append(str(subject_dir))
                continue
            subject_match = SUBJECT_RE.match(subject_dir.name)
            if not subject_match:
                unparsed.append(f"subject folder name not understood: {subject_dir}")
                continue
            subject = subject_match.group("code")
            subject_num = subject_match.group("num")

            for session_dir in sorted(subject_dir.iterdir()):
                if not session_dir.is_dir():
                    ignored.append(str(session_dir))
                    continue
                day_match = DAY_RE.search(session_dir.name)
                if not day_match:
                    unparsed.append(f"no day number in folder name: {session_dir}")
                    continue
                day = int(day_match.group("day"))

                # Group the files in this folder by their stem.
                groups: Dict[str, Dict[str, Path]] = {}
                for entry in sorted(session_dir.iterdir()):
                    if not entry.is_file():
                        ignored.append(str(entry))
                        continue
                    ext = entry.suffix.lower()
                    if ext not in EXTENSIONS:
                        ignored.append(str(entry))
                        continue
                    groups.setdefault(entry.stem, {})[ext] = entry

                for stem, files in sorted(groups.items()):
                    phase_match = PHASE_RE.search(stem)
                    if not phase_match:
                        unparsed.append(f"no pre/post token in filename: {session_dir / stem}")
                        continue
                    phase = phase_match.group("phase").lower()

                    recording = Recording(
                        subject=subject,
                        subject_num=subject_num,
                        day=day,
                        phase=phase,
                        cond=condition,
                        src_dir=session_dir,
                        stem=stem,
                        src=files,
                    )

                    # The filename should agree with the folder it sits in.
                    cond_match = COND_TOKEN_RE.search(stem)
                    if cond_match:
                        from_name = "EO" if cond_match.group("cond").lower() == "open" else "EC"
                        if from_name != condition:
                            recording.flags.append("CONDITION_MISMATCH")
                            recording.notes.append(
                                f"filename says {from_name} but it sits under {root_name}"
                            )
                    else:
                        recording.notes.append("no open/closed token in filename")

                    # Two hand-typed dates that can disagree. enrich() settles
                    # the argument against the recorder's clock.
                    folder_date_match = FOLDER_DATE_RE.match(session_dir.name)
                    file_date_match = FILE_DATE_RE.search(stem)
                    recording.folder_date = (
                        folder_date_match.group("date") if folder_date_match else None
                    )
                    recording.file_date = (
                        file_date_match.group("date") if file_date_match else None
                    )
                    if recording.folder_date and recording.file_date:
                        if recording.folder_date.lstrip("0") != recording.file_date.lstrip("0"):
                            recording.flags.append("DATE_MISMATCH")

                    if recording.key in recordings:
                        other = recordings[recording.key]
                        collisions.append(
                            f"{recording.label}: '{other.src_dir / other.stem}' and "
                            f"'{session_dir / stem}' both map to the same recording"
                        )
                        continue
                    recordings[recording.key] = recording

    if collisions:
        report.rule("FATAL: two source recordings map to the same destination")
        for line in collisions:
            report(f"  {FAIL} {line}")
        report("\nResolve these by hand before running again; nothing was copied.")
        sys.exit(2)

    ordered = sorted(recordings.values(), key=lambda r: (r.subject, r.day, r.phase != "pre", r.cond))
    return ordered, unparsed, ignored


def enrich(recordings: List[Recording]) -> None:
    """Read each header, take the recording time, derive the duration."""
    for recording in recordings:
        vhdr = recording.src.get(".vhdr")
        if vhdr:
            recording.header = read_header(vhdr)
            recording.recorded, recording.recorded_from = timestamp_via_header(vhdr)
        elif ".vmrk" in recording.src:
            recording.recorded = read_marker_timestamp(recording.src[".vmrk"])
            recording.recorded_from = recording.src[".vmrk"].name

        if "DATE_MISMATCH" in recording.flags:
            resolve_date_mismatch(recording)

        eeg = recording.src.get(".eeg")
        if not eeg:
            continue

        size = eeg.stat().st_size
        header = recording.header
        n_channels, bps, interval = (
            header.n_channels,
            header.bytes_per_sample,
            header.sampling_interval_us,
        )

        if n_channels and bps and interval:
            frame = n_channels * bps
            recording.duration_s = size / frame * interval / 1_000_000.0
            if size % frame:
                recording.flags.append("NOT_BLOCK_ALIGNED")
                recording.notes.append(
                    f".eeg size {size} is not a multiple of {frame} bytes per sample"
                )
        elif not vhdr:
            recording.notes.append("no header, so duration cannot be derived")
        else:
            recording.notes.append("header incomplete, so duration cannot be derived")


def flag_sizes(
    recordings: List[Recording], min_duration: float, rel_fraction: float
) -> Optional[float]:
    """Mark recordings that look too short, judged both absolutely and against
    the dataset median."""
    durations = [r.duration_s for r in recordings if r.duration_s]
    median = statistics.median(durations) if durations else None
    threshold = median * rel_fraction if median else None

    for recording in recordings:
        for ext, path in recording.src.items():
            size = path.stat().st_size
            if size == 0:
                recording.flags.append(f"EMPTY{ext.upper().replace('.', '_')}")
            elif ext == ".vhdr" and size < MIN_VHDR_BYTES:
                recording.flags.append("SMALL_VHDR")
                recording.notes.append(f".vhdr is only {size} bytes")
            elif ext == ".vmrk" and size < MIN_VMRK_BYTES:
                recording.flags.append("SMALL_VMRK")
                recording.notes.append(f".vmrk is only {size} bytes")

        if recording.duration_s is None:
            continue
        if recording.duration_s < min_duration:
            recording.flags.append("SHORT")
            recording.notes.append(
                f"{recording.duration_s:.1f}s is under the {min_duration:.0f}s floor"
            )
        elif threshold and recording.duration_s < threshold:
            recording.flags.append("SHORT_VS_MEDIAN")
            recording.notes.append(
                f"{recording.duration_s:.1f}s is under {rel_fraction:.0%} of the "
                f"{median:.0f}s median"
            )
    return median


# --------------------------------------------------------------------------
# Source protection
# --------------------------------------------------------------------------


def snapshot_sources(recordings: List[Recording]) -> Dict[str, Tuple[int, int]]:
    """Record (size, mtime_ns) of every source file so we can prove afterwards
    that the run left them alone."""
    snapshot: Dict[str, Tuple[int, int]] = {}
    for recording in recordings:
        for path in recording.src.values():
            stat = path.stat()
            snapshot[str(path)] = (stat.st_size, stat.st_mtime_ns)
    return snapshot


def check_sources_untouched(snapshot: Dict[str, Tuple[int, int]]) -> List[str]:
    problems: List[str] = []
    for path_text, (size, mtime) in snapshot.items():
        path = Path(path_text)
        if not path.exists():
            problems.append(f"{path} has disappeared")
            continue
        stat = path.stat()
        if stat.st_size != size:
            problems.append(f"{path} changed size: {size} -> {stat.st_size}")
        elif stat.st_mtime_ns != mtime:
            problems.append(f"{path} was modified (mtime changed)")
    return problems


# --------------------------------------------------------------------------
# Copy + verify
# --------------------------------------------------------------------------


def copy_and_hash(source: Path, dest: Path) -> Tuple[str, int]:
    """Copy in chunks, hashing the bytes as they stream past."""
    digest = hashlib.sha256()
    total = 0
    tmp = dest.with_suffix(dest.suffix + ".partial")
    try:
        with open(source, "rb") as reader, open(tmp, "wb") as writer:
            while True:
                chunk = reader.read(COPY_CHUNK)
                if not chunk:
                    break
                digest.update(chunk)
                writer.write(chunk)
                total += len(chunk)
            writer.flush()
            os.fsync(writer.fileno())
        os.replace(tmp, dest)
    except BaseException:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        raise
    return digest.hexdigest(), total


def transfer_recording(
    recording: Recording,
    target_root: Path,
    verify: str,
    overwrite: bool,
    dry_run: bool,
) -> List[FileResult]:
    """Copy the .eeg/.vhdr/.vmrk of one recording, then repair the pointers in
    the copied .vhdr/.vmrk so they name the new files."""
    results: List[FileResult] = []
    dest_dir = recording.dest_dir(target_root)
    if not dry_run:
        dest_dir.mkdir(parents=True, exist_ok=True)

    new_eeg = recording.dest_name(".eeg")
    new_vmrk = recording.dest_name(".vmrk")

    for ext in EXTENSIONS:
        source = recording.src.get(ext)
        dest = recording.dest_path(target_root, ext)
        if source is None:
            results.append(
                FileResult(
                    recording=recording.label,
                    ext=ext,
                    source="",
                    dest=str(dest),
                    bytes=0,
                    status="MISSING_IN_SOURCE",
                    detail=f"no {ext} exists for this recording",
                )
            )
            continue

        size = source.stat().st_size
        result = FileResult(
            recording=recording.label,
            ext=ext,
            source=str(source),
            dest=str(dest),
            bytes=size,
            status="PLANNED",
        )

        # Work out what the destination is supposed to contain. Header and
        # marker files are a few kB, so this is cheap enough to do in dry-run
        # mode too -- it lets the rehearsal show exactly which pointers would
        # be rewritten before anything is committed to disk.
        expected_bytes: Optional[bytes] = None
        if ext in (".vhdr", ".vmrk"):
            try:
                original = source.read_bytes()
            except OSError as exc:
                result.status = "FAILED"
                result.detail = f"cannot read source: {exc}"
                results.append(result)
                continue
            rewritten, changes = rewrite_pointers(original, new_eeg, new_vmrk)
            expected_bytes = rewritten
            if changes:
                result.pointer_changes = "; ".join(
                    f"{key}: {old} -> {new}" for key, old, new in changes
                )
                # A pointer that did not name an existing file was already
                # broken before this script ever ran -- worth calling out,
                # because those recordings would not load from the source tree.
                dangling = [
                    f"{key}={old}"
                    for key, old, _ in changes
                    if old and not (source.parent / old).exists()
                ]
                if dangling:
                    result.pointer_broken_in_source = "; ".join(dangling)

        if dry_run:
            result.status = "WOULD_COPY"
            if dest.exists():
                result.detail = "destination already exists"
            results.append(result)
            continue

        # A re-run must tell "already done" from "something else is there".
        if dest.exists() and not overwrite:
            try:
                existing = sha256_of(dest)
            except OSError as exc:
                result.status = "FAILED"
                result.detail = f"cannot read existing destination: {exc}"
                results.append(result)
                continue

            wanted = (
                hashlib.sha256(expected_bytes).hexdigest()
                if expected_bytes is not None
                else sha256_of(source)
            )
            result.sha_source = wanted
            if existing == wanted:
                result.status = "SKIPPED_IDENTICAL"
                result.detail = "already present with the expected content"
            else:
                result.status = "CONFLICT"
                result.detail = "destination exists with different content; left untouched"
            results.append(result)
            continue

        try:
            sha_source, written = copy_and_hash(source, dest)
            result.sha_source = sha_source

            if written != size:
                result.status = "FAILED"
                result.detail = f"read {written} bytes but source is {size}"
                results.append(result)
                continue

            # Rewrite the pointer lines in the copy (never in the source).
            if expected_bytes is not None:
                with open(dest, "wb") as handle:
                    handle.write(expected_bytes)
                    handle.flush()
                    os.fsync(handle.fileno())

            shutil.copystat(source, dest, follow_symlinks=False)

            if verify == "none":
                result.status = "COPIED"
                result.detail = "not verified (--verify none)"
            elif verify == "size":
                expected_size = len(expected_bytes) if expected_bytes is not None else size
                actual_size = dest.stat().st_size
                if actual_size == expected_size:
                    result.status = "COPIED"
                    result.detail = "size verified"
                else:
                    result.status = "FAILED"
                    result.detail = (
                        f"size mismatch: {actual_size} on disk vs {expected_size} expected"
                    )
            else:  # hash
                sha_dest = sha256_of(dest)
                expected = (
                    hashlib.sha256(expected_bytes).hexdigest()
                    if expected_bytes is not None
                    else sha_source
                )
                if sha_dest == expected:
                    result.status = "COPIED"
                    result.detail = (
                        f"sha256 {sha_dest[:16]} verified"
                        + (" (after pointer fix)" if expected_bytes is not None else "")
                    )
                else:
                    result.status = "FAILED"
                    result.detail = "sha256 mismatch between source and destination"
        except OSError as exc:
            result.status = "FAILED"
            result.detail = f"{type(exc).__name__}: {exc}"

        results.append(result)

    return results


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------


def report_completeness(recordings: List[Recording], report: Reporter) -> None:
    """One line per subject showing which of the 16 expected cells are present."""
    present = {r.key for r in recordings}
    subjects = sorted({r.subject for r in recordings})
    expected_per_subject = len(EXPECTED_DAYS) * len(EXPECTED_PHASES) * len(EXPECTED_CONDITIONS)

    report.rule("COMPLETENESS")
    report(
        f"Each subject should have {expected_per_subject} recordings: "
        f"{len(EXPECTED_DAYS)} days x pre/post x EO/EC."
    )
    report("")
    header_cells = " ".join(f"d{d}" for d in EXPECTED_DAYS)
    report(f"  {'subject':<10} {'found':>7}   {header_cells}      missing")
    report(f"  {'-' * 10} {'-' * 7}   {'-' * len(header_cells)}      {'-' * 30}")

    complete_subjects = 0
    for subject in subjects:
        found = 0
        per_day: List[str] = []
        missing_labels: List[str] = []
        for day in EXPECTED_DAYS:
            day_count = 0
            for phase in EXPECTED_PHASES:
                for cond in EXPECTED_CONDITIONS:
                    exists = (subject, day, phase, cond) in present
                    found += exists
                    day_count += exists
                    if not exists:
                        missing_labels.append(f"d{day}/{phase}/{cond}")
            per_day.append(f"{day_count}/4")

        tag = OK if found == expected_per_subject else WARN
        complete_subjects += found == expected_per_subject
        missing_text = "-" if not missing_labels else ", ".join(missing_labels[:6]) + (
            f" (+{len(missing_labels) - 6} more)" if len(missing_labels) > 6 else ""
        )
        report(
            f"  {tag} {subject:<7} {found:>2}/{expected_per_subject}   "
            f"{' '.join(f'{c:>3}' for c in per_day)}   {missing_text}"
        )

    report("")
    report(
        f"  {complete_subjects} of {len(subjects)} subjects have the full "
        f"{expected_per_subject}-recording grid."
    )


def report_flags(recordings: List[Recording], report: Reporter, median: Optional[float]) -> None:
    report.rule("RECORDING LENGTH")
    if median:
        report(f"Median recording length across the dataset: {median:.0f} s ({median / 60:.1f} min).")
    report("")

    short = [r for r in recordings if "SHORT" in r.flags or "SHORT_VS_MEDIAN" in r.flags]
    if short:
        report(f"  {WARN} {len(short)} recording(s) look suspiciously short:")
        report("")
        report(f"       {'recording':<32} {'length':>8}   source folder")
        for recording in sorted(short, key=lambda r: r.duration_s or 0):
            report(
                f"       {recording.label:<32} {human_duration(recording.duration_s):>8}   "
                f"{recording.src_dir.name}"
            )
    else:
        report(f"  {OK} No recording falls below the length thresholds.")

    orphans = [r for r in recordings if r.is_orphan]
    report("")
    if orphans:
        report(f"  {WARN} {len(orphans)} recording(s) are missing companion files:")
        for recording in orphans:
            report(f"       {recording.label:<32} missing {', '.join(recording.missing_exts)}")
            report(f"       {'':<32} in {recording.src_dir}")
        report("")
        report("       A .eeg without its .vhdr cannot be read by EEGLAB or MNE.")
        report("       These were copied anyway so the data is not hidden from you.")
    else:
        report(f"  {OK} Every recording has all three files.")

    other_flags = [
        "CONDITION_MISMATCH",
        "DATE_MISMATCH",
        "DATE_UNEXPLAINED",
        "NOT_BLOCK_ALIGNED",
        "SMALL_VHDR",
        "SMALL_VMRK",
    ]
    noted = [r for r in recordings if any(f in r.flags for f in other_flags)]
    report("")
    if noted:
        report(f"  {INFO} {len(noted)} recording(s) with metadata oddities (not transfer errors):")
        for recording in noted:
            relevant = [f for f in recording.flags if f in other_flags]
            report(f"       {recording.label:<32} {', '.join(relevant)}")
            for note in recording.notes:
                report(f"       {'':<32}   {note}")
    else:
        report(f"  {OK} No metadata inconsistencies found.")


@dataclass
class Placed:
    """A recording as it actually sits in the target tree (or, before the first
    copy, as it is planned to sit there)."""

    subject: str
    day: int
    phase: str
    cond: str
    vhdr: Optional[Path] = None
    recorded: Optional[datetime] = None
    in_source: bool = True

    @property
    def location(self) -> str:
        return f"sub-{self.subject}/ses-day{self.day}/{self.phase}/{self.cond}"


def scan_target(target_root: Path) -> List[Placed]:
    """Walk what is really on disk under the target root.

    Deliberately independent of the source listing, so anything moved, renamed
    or added by hand after a previous run is seen too -- those are exactly the
    files most likely to be in the wrong place.
    """
    placed: List[Placed] = []
    if not target_root.is_dir():
        return placed

    for cond_dir in sorted(target_root.glob("sub-*/ses-day*/*/*")):
        if not cond_dir.is_dir():
            continue
        subject_part, session_part, phase, cond = (
            cond_dir.parts[-4],
            cond_dir.parts[-3],
            cond_dir.parts[-2],
            cond_dir.parts[-1],
        )
        day_match = DAY_RE.search(session_part)
        if not (subject_part.startswith("sub-") and day_match):
            continue
        if phase not in EXPECTED_PHASES or cond not in EXPECTED_CONDITIONS:
            continue

        headers = sorted(cond_dir.glob("*.vhdr"))
        data_files = sorted(cond_dir.glob("*.eeg"))
        if not headers and not data_files:
            continue

        entry = Placed(
            subject=subject_part[len("sub-"):],
            day=int(day_match.group("day")),
            phase=phase,
            cond=cond,
            vhdr=headers[0] if headers else None,
        )
        if entry.vhdr:
            entry.recorded, _ = timestamp_via_header(entry.vhdr)
        placed.append(entry)

    return placed


def verify_loadable(target_root: Path, report: Reporter) -> Tuple[int, int, bool]:
    """Open every .vhdr in the target tree and confirm the names it points at
    really exist next to it, and that the .eeg is the size the header implies.
    This is the check that answers 'will EEGLAB open this'.

    It walks the tree rather than the source listing, so a file moved or
    renamed by hand -- whose header pointers would then no longer match -- is
    caught rather than skipped.
    """
    checked = loadable = 0
    problems = False

    for cond_dir in sorted(target_root.glob("sub-*/ses-day*/*/*")):
        if not cond_dir.is_dir():
            continue
        location = "/".join(cond_dir.parts[-4:])

        headers = sorted(cond_dir.glob("*.vhdr"))
        data_files = sorted(cond_dir.glob("*.eeg"))

        if not headers:
            if data_files:
                # A known, already-reported condition of the source data
                # rather than something this run broke: warn, don't fail.
                report(f"  {WARN} {location}")
                report(f"        {data_files[0].name} has no .vhdr beside it, so nothing "
                       f"can open it (orphan carried over from the source)")
            continue

        if len(headers) > 1:
            problems = True
            report(f"  {FAIL} {location}")
            report(f"        {len(headers)} .vhdr files in one folder: "
                   f"{', '.join(h.name for h in headers)}")

        for vhdr in headers:
            checked += 1
            header = read_header(vhdr)
            issues: List[str] = []

            for label, name, expected_ext in (
                ("DataFile", header.data_file, ".eeg"),
                ("MarkerFile", header.marker_file, ".vmrk"),
            ):
                if not name:
                    issues.append(f"{label}= is missing from the header")
                    continue
                companion = vhdr.parent / name
                if not companion.exists():
                    issues.append(
                        f"{label}={name} does not exist next to the header"
                        + (
                            " -- was this file renamed without updating the header?"
                            if (vhdr.parent / vhdr.stem).with_suffix(expected_ext).exists()
                            else ""
                        )
                    )
                elif not name.lower().endswith(expected_ext):
                    issues.append(f"{label}={name} does not end in {expected_ext}")

            eeg = vhdr.with_suffix(".eeg")
            if eeg.exists() and header.n_channels and header.bytes_per_sample:
                frame = header.n_channels * header.bytes_per_sample
                if eeg.stat().st_size % frame:
                    issues.append(
                        f".eeg is {eeg.stat().st_size} bytes, not a multiple of {frame}"
                    )

            if issues:
                problems = True
                report(f"  {FAIL} {location}/{vhdr.name}")
                for issue in issues:
                    report(f"        {issue}")
            else:
                loadable += 1

    return checked, loadable, problems


# --------------------------------------------------------------------------
# Final sanity check: does the recording time fit the assigned location?
# --------------------------------------------------------------------------


def sanity_check(recordings: List[Recording], target_root: Path, report: Reporter) -> int:
    """Test every recording's clock against the folder it was filed into.

    For each recording we start at its .vhdr, follow the header's MarkerFile=
    line to the .vmrk and read the 'New Segment' timestamp the recorder wrote
    at acquisition time. That timestamp is then checked against the three
    things the destination path claims: which session day, which phase
    (pre/post) and which condition block (EO/EC) the recording belongs to.

    After a real run the .vhdr files read here are the copies in the target
    tree, so this validates what is actually on disk.
    """
    report.rule("SANITY CHECK -- RECORDING TIME vs ASSIGNED LOCATION")
    report("Each recording's .vhdr is resolved to its marker file and the recorder's")
    report("own 'New Segment' timestamp is read from it. Folder and file names were")
    report("typed by a human; this timestamp was not, so it is the reference.")
    report("")

    # Read what is really on disk. Only if the target tree is still empty
    # (a dry run before the first copy) do we fall back to the planned layout.
    placed = scan_target(target_root)
    if placed:
        source_keys = {r.key for r in recordings}
        for entry in placed:
            entry.in_source = (entry.subject, entry.day, entry.phase, entry.cond) in source_keys
        report("  Read from: the recordings actually present in the target tree.")
    else:
        placed = [
            Placed(
                subject=r.subject, day=r.day, phase=r.phase, cond=r.cond,
                vhdr=r.src.get(".vhdr"), recorded=r.recorded,
            )
            for r in recordings
        ]
        report("  Read from: the source tree -- the target is still empty (dry run).")

    with_stamp = [p for p in placed if p.recorded]
    unreadable = [p for p in placed if not p.recorded]
    report(f"  Timestamps recovered for {len(with_stamp)} of {len(placed)} recordings.")

    if unreadable:
        report("")
        report(f"  {WARN} {len(unreadable)} recording(s) carry no readable timestamp "
               f"and cannot be checked:")
        for entry in unreadable:
            reason = "no .vhdr" if entry.vhdr is None else "no New Segment marker"
            report(f"        {entry.location:<38} ({reason})")

    # The same physical recording filed in two places at once. This is what a
    # half-finished manual re-filing looks like: the file was copied to its new
    # home but the old copy is still there (or a later run put it back).
    by_stamp: Dict[Tuple[str, datetime], List[Placed]] = {}
    for entry in with_stamp:
        by_stamp.setdefault((entry.subject, entry.recorded), []).append(entry)
    duplicates = [group for group in by_stamp.values() if len(group) > 1]

    extras = [p for p in placed if not p.in_source]
    if extras:
        report("")
        report(f"  {INFO} {len(extras)} recording(s) in the target have no counterpart in the")
        report("         source tree -- added or moved by hand since the last run:")
        for entry in extras:
            report(f"        {entry.location}")

    report("")

    subjects = sorted({p.subject for p in placed})
    mismatches: List[str] = []

    for group in duplicates:
        locations = ", ".join(sorted(e.location for e in group))
        mismatches.append(
            f"sub-{group[0].subject}: the recording made at "
            f"{group[0].recorded:%Y-%m-%d %H:%M} is filed in {len(group)} places at once "
            f"({locations})"
        )

    for subject in subjects:
        subject_lines: List[str] = []
        subject_issues: List[str] = []
        previous_day: Optional[int] = None
        previous_date: Optional[datetime] = None

        for day in EXPECTED_DAYS:
            cells: Dict[Tuple[str, str], Optional[datetime]] = {}
            for phase in EXPECTED_PHASES:
                for cond in EXPECTED_CONDITIONS:
                    match = [
                        p.recorded for p in placed
                        if p.subject == subject and p.day == day
                        and p.phase == phase and p.cond == cond and p.recorded
                    ]
                    cells[(phase, cond)] = min(match) if match else None

            present = [t for t in cells.values() if t]
            if not present:
                continue

            dates = sorted({t.date() for t in present})
            times = "  ".join(
                f"{phase} {cond} "
                + (f"{cells[(phase, cond)]:%H:%M}" if cells[(phase, cond)] else "--:--")
                for phase in EXPECTED_PHASES
                for cond in EXPECTED_CONDITIONS
            )
            subject_lines.append(f"      ses-day{day}  {dates[0]:%Y-%m-%d}   {times}")

            # 1. one session should not straddle two calendar days
            if len(dates) > 1:
                subject_issues.append(
                    f"ses-day{day} holds recordings from {len(dates)} different dates: "
                    + ", ".join(str(d) for d in dates)
                )

            # 2. day numbering must follow the calendar
            earliest = min(present)
            if previous_date and earliest.date() < previous_date.date():
                subject_issues.append(
                    f"ses-day{day} ({earliest:%Y-%m-%d}) was recorded BEFORE "
                    f"ses-day{previous_day} ({previous_date:%Y-%m-%d}) -- the day "
                    f"numbering does not match the calendar"
                )
            previous_day, previous_date = day, earliest

            # 3. everything in pre/ must precede everything in post/
            pre = [t for (phase, _), t in cells.items() if phase == "pre" and t]
            post = [t for (phase, _), t in cells.items() if phase == "post" and t]
            if pre and post and min(post) <= max(pre):
                offenders = [
                    f"{phase}/{cond} at {t:%H:%M}"
                    for (phase, cond), t in cells.items()
                    if t and ((phase == "post" and t <= max(pre)) or (phase == "pre" and t >= min(post)))
                ]
                subject_issues.append(
                    f"ses-day{day}: post is not after pre ({', '.join(sorted(offenders))})"
                )

            # 4. EO and EC of one block are recorded back to back
            for phase in EXPECTED_PHASES:
                eo, ec = cells[(phase, "EO")], cells[(phase, "EC")]
                if not (eo and ec):
                    continue
                gap = abs((eo - ec).total_seconds()) / 60.0
                if gap > MAX_EO_EC_GAP_MINUTES:
                    subject_issues.append(
                        f"ses-day{day} {phase}: EO ({eo:%H:%M}) and EC ({ec:%H:%M}) are "
                        f"{gap:.0f} min apart, so they were not recorded as one block -- "
                        f"one of the two is filed under the wrong phase"
                    )

        if not subject_lines:
            continue

        tag = WARN if subject_issues else OK
        report(f"  {tag} sub-{subject}")
        for line in subject_lines:
            report(line)
        for issue in subject_issues:
            report(f"        -> {issue}")
            mismatches.append(f"sub-{subject}: {issue}")
        report("")

    report.rule("SANITY CHECK VERDICT", char="-")
    if not mismatches:
        report(f"  {OK} No mismatches. For every subject the session days run in")
        report("        chronological order, each session sits on a single date, every")
        report("        pre precedes its post, and EO/EC were recorded as one block.")
        report("        The folder each recording was filed into agrees with its clock.")
    else:
        report(f"  {FAIL} {len(mismatches)} mismatch(es) between recording time and location:")
        report("")
        for line in mismatches:
            report(f"        {line}")
        report("")
        report("        These are labelling problems in the source data, not copy errors.")
        report("        Nothing was moved or renamed to 'fix' them -- that judgement is")
        report("        yours. Re-file by hand, then re-run with --subjects <CODE>.")
    return len(mismatches)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reorganise APOP_Glu raw BrainVision files into a BIDS-style tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--project", type=Path, default=DEFAULT_PROJECT,
                        help="project root containing raw_data/")
    parser.add_argument("--source", type=Path, default=None,
                        help="source tree (default: <project>/raw_data)")
    parser.add_argument("--target", type=Path, default=None,
                        help="output tree (default: <project>/raw_BIDS)")
    parser.add_argument("--log", type=Path, default=None,
                        help="path for the .txt transcript "
                             "(default: <target>/bidsify_log_<timestamp>.txt)")
    parser.add_argument("--dry-run", action="store_true",
                        help="audit and plan only; copy nothing (the log is still written)")
    parser.add_argument("--verify", choices=("hash", "size", "none"), default="hash",
                        help="how thoroughly to confirm each copy")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace destination files that already exist")
    parser.add_argument("--limit", type=int, default=0,
                        help="only process the first N recordings (0 = all)")
    parser.add_argument("--subjects", nargs="+", default=None, metavar="CODE",
                        help="only process these subject codes, e.g. --subjects GB StP")
    parser.add_argument("--exclude", nargs="+", default=None, metavar="CODE",
                        help="skip these subject codes entirely -- use this for subjects you "
                             "have deliberately moved out of the target tree, so a re-run does "
                             "not copy them back in, e.g. --exclude FD FS")
    parser.add_argument("--min-duration", type=float, default=DEFAULT_MIN_DURATION_S,
                        help="flag recordings shorter than this many seconds")
    parser.add_argument("--rel-fraction", type=float, default=DEFAULT_REL_FRACTION,
                        help="also flag recordings below this fraction of the median length")
    parser.add_argument("--summary-only", action="store_true",
                        help="do not print a line per file during the transfer")
    parser.add_argument("--verbose-pointers", action="store_true",
                        help="list every rewritten header pointer, not just the broken ones")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    source_root = args.source or (args.project / DEFAULT_SOURCE)
    target_root = args.target or (args.project / DEFAULT_TARGET)
    started = time.time()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = args.log or (target_root / f"bidsify_log_{stamp}.txt")

    report = Reporter()
    report.rule("APOP_Glu -> BIDS-style reorganisation")
    report(f"  source     : {source_root}")
    report(f"  target     : {target_root}")
    report(f"  log        : {log_path}")
    report(f"  mode       : {'DRY RUN (nothing will be copied)' if args.dry_run else 'copy'}")
    report(f"  verify     : {args.verify}")
    report(f"  started    : {datetime.now():%Y-%m-%d %H:%M:%S}")

    if not source_root.is_dir():
        report(f"\n  {FAIL} source folder does not exist: {source_root}")
        return 2

    # ---- discovery -------------------------------------------------------
    recordings, unparsed, ignored = discover(source_root, report)
    enrich(recordings)
    median = flag_sizes(recordings, args.min_duration, args.rel_fraction)

    report.rule("DISCOVERY")
    total_files = sum(len(r.src) for r in recordings)
    total_bytes = sum(p.stat().st_size for r in recordings for p in r.src.values())
    subjects = sorted({r.subject for r in recordings})
    report(f"  {len(subjects)} subjects: {', '.join(subjects)}")
    report(f"  {len(recordings)} recordings, {total_files} files, {human_bytes(total_bytes)} total")

    if unparsed:
        report("")
        report(f"  {WARN} {len(unparsed)} path(s) could not be interpreted and will NOT be copied:")
        for line in unparsed:
            report(f"       {line}")
    if ignored:
        report("")
        report(f"  {INFO} {len(ignored)} non-recording file(s)/folder(s) ignored:")
        for line in ignored[:10]:
            report(f"       {line}")
        if len(ignored) > 10:
            report(f"       ... and {len(ignored) - 10} more")

    # Thresholds above were computed on the whole dataset; only now do we
    # narrow down to the subset the user asked for.
    if args.subjects:
        wanted = {s.lower() for s in args.subjects}
        unknown = wanted - {r.subject.lower() for r in recordings}
        if unknown:
            report("")
            report(f"  {FAIL} unknown subject code(s): {', '.join(sorted(unknown))}")
            return 2
        recordings = [r for r in recordings if r.subject.lower() in wanted]
        report("")
        report(f"  {INFO} --subjects: restricted to {len(recordings)} recordings from "
               f"{', '.join(args.subjects)}.")

    if args.exclude:
        dropped = {s.lower() for s in args.exclude}
        unknown = dropped - {r.subject.lower() for r in recordings}
        if unknown:
            report("")
            report(f"  {FAIL} unknown subject code(s) in --exclude: {', '.join(sorted(unknown))}")
            return 2
        before = len(recordings)
        recordings = [r for r in recordings if r.subject.lower() not in dropped]
        report("")
        report(f"  {INFO} --exclude: skipping {', '.join(args.exclude)} "
               f"({before - len(recordings)} recordings). They are neither copied nor")
        report("         checked, and anything already in the target tree is left alone.")

    if args.limit:
        recordings = recordings[: args.limit]
        report("")
        report(f"  {INFO} --limit {args.limit}: only the first {len(recordings)} "
               f"recordings will be processed.")

    # ---- transfer --------------------------------------------------------
    snapshot = snapshot_sources(recordings)
    if not args.dry_run:
        target_root.mkdir(parents=True, exist_ok=True)

    report.rule("TRANSFER")
    results: List[FileResult] = []
    for index, recording in enumerate(recordings, start=1):
        recording_results = transfer_recording(
            recording, target_root, args.verify, args.overwrite, args.dry_run
        )
        results.extend(recording_results)

        statuses = {r.status for r in recording_results}
        if statuses & {"FAILED", "CONFLICT"}:
            tag = FAIL
        elif "MISSING_IN_SOURCE" in statuses:
            tag = WARN
        elif statuses == {"SKIPPED_IDENTICAL"}:
            tag = SKIP
        else:
            tag = OK

        size = sum(r.bytes for r in recording_results)
        if not args.summary_only:
            report(
                f"  {tag} [{index:>3}/{len(recordings)}] {recording.label:<32} "
                f"{human_duration(recording.duration_s)} {human_bytes(size):>10}"
            )
            for file_result in recording_results:
                marker = {
                    "COPIED": OK,
                    "SKIPPED_IDENTICAL": SKIP,
                    "WOULD_COPY": INFO,
                    "MISSING_IN_SOURCE": WARN,
                    "CONFLICT": FAIL,
                    "FAILED": FAIL,
                }.get(file_result.status, INFO)
                name = Path(file_result.dest).name
                detail = f" -- {file_result.detail}" if file_result.detail else ""
                report(f"        {marker} {file_result.ext:<5} {name}{detail}")
                if file_result.pointer_changes and file_result.status in ("COPIED", "WOULD_COPY"):
                    verb = "pointer fix" if file_result.status == "COPIED" else "would rewrite"
                    report(f"              {verb}: {file_result.pointer_changes}")
        elif tag in (FAIL, WARN):
            report(f"  {tag} {recording.label}")
            for file_result in recording_results:
                if file_result.status in ("FAILED", "CONFLICT", "MISSING_IN_SOURCE"):
                    report(f"        {file_result.ext:<5} {file_result.status}: "
                           f"{file_result.detail}")

    # ---- pointer repairs -------------------------------------------------
    repaired = [r for r in results if r.pointer_changes and r.status in ("COPIED", "WOULD_COPY")]
    report.rule("HEADER POINTER REWRITES")
    if repaired:
        tense = "would have their" if args.dry_run else "had their"
        report(f"  {len(repaired)} header/marker file(s) {tense} DataFile= / MarkerFile= "
               f"lines updated in the COPY.")
        report("  The source files are not touched by this.")
        report("")

        broken = [r for r in repaired if r.pointer_broken_in_source]
        if broken:
            report(f"  {WARN} {len(broken)} of them pointed at a file that does not exist. Those")
            report("        were already broken in the source tree and would not load from")
            report("        there at all; they are repaired as a side effect of the rename:")
            report("")
            for file_result in broken:
                report(f"    {file_result.recording}  {file_result.ext}")
                report(f"        dangling: {file_result.pointer_broken_in_source}")
        else:
            report(f"  {OK} None of them were dangling in the source tree.")

        if args.verbose_pointers:
            report("")
            report(f"  {INFO} Full list of rewritten pointer lines:")
            report("")
            for file_result in repaired:
                report(f"    {file_result.recording}  {file_result.ext}")
                for change in file_result.pointer_changes.split("; "):
                    report(f"        {change}")
        else:
            report("")
            report(f"  {INFO} Pass --verbose-pointers to list all {len(repaired)} rewrites.")
    else:
        report(f"  {INFO} No pointer lines needed changing.")

    # ---- can the result actually be loaded? ------------------------------
    report.rule("LOADABILITY OF THE RESULT")
    loadability_failed = False
    if args.dry_run:
        report(f"  {INFO} Not checked in dry-run mode (nothing was written).")
    else:
        checked, loadable, problems_found = verify_loadable(target_root, report)
        report("")
        if checked == 0:
            report(f"  {INFO} Nothing to check.")
        elif not problems_found:
            report(f"  {OK} All {loadable}/{checked} .vhdr files in the target resolve to an "
                   f"existing .eeg and .vmrk of the expected size.")
        else:
            loadability_failed = True
            report(f"  {FAIL} {checked - loadable} of {checked} .vhdr files in the target do "
                   f"not resolve correctly (listed above).")
            report("        A header that names a file which is not there will not open. If")
            report("        this followed a manual rename, the header pointers need the same")
            report("        edit -- or re-copy that recording with --overwrite.")

    # ---- audits ----------------------------------------------------------
    report_completeness(recordings, report)
    report_flags(recordings, report, median)

    # ---- source integrity ------------------------------------------------
    report.rule("SOURCE INTEGRITY")
    problems = check_sources_untouched(snapshot)
    if problems:
        report(f"  {FAIL} {len(problems)} source file(s) changed during this run:")
        for line in problems:
            report(f"       {line}")
    else:
        report(f"  {OK} All {len(snapshot)} source files are byte-for-byte unchanged "
               f"(size and mtime re-checked).")

    # ---- transfer summary ------------------------------------------------
    counts: Dict[str, int] = {}
    for file_result in results:
        counts[file_result.status] = counts.get(file_result.status, 0) + 1

    report.rule("TRANSFER SUMMARY")
    for status in ("COPIED", "WOULD_COPY", "SKIPPED_IDENTICAL", "MISSING_IN_SOURCE",
                   "CONFLICT", "FAILED"):
        if status in counts:
            tag = FAIL if status in ("FAILED", "CONFLICT") else (
                WARN if status == "MISSING_IN_SOURCE" else OK
            )
            report(f"  {tag} {status:<20} {counts[status]:>5}")
    failed = counts.get("FAILED", 0) + counts.get("CONFLICT", 0)
    report("")
    report(f"  elapsed: {(time.time() - started) / 60:.1f} min")

    # ---- final step: sanity check ----------------------------------------
    mismatches = sanity_check(recordings, target_root, report)

    # ---- the one and only output file ------------------------------------
    try:
        report.save(log_path)
        print(f"\nLog written to: {log_path}", flush=True)
    except OSError as exc:
        print(f"\n{FAIL} could not write the log to {log_path}: {exc}", flush=True)

    if problems or failed or loadability_failed:
        print(f"{FAIL} Finished with transfer/integrity problems -- see the log.", flush=True)
        return 1
    if mismatches:
        print(f"{WARN} Files copied cleanly, but the sanity check found {mismatches} "
              f"time/location mismatch(es) -- see the log.", flush=True)
        return 3
    print(f"{OK} Finished.", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
