"""Public reference data and CSV format checks for Changjiang Cup 2026."""

import argparse
import csv
from dataclasses import dataclass
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
import sys
from tempfile import NamedTemporaryFile
from textwrap import fill
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
TARGETS = ROOT / "data" / "targets"
SUBMISSIONS = ROOT / "submissions"
AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
LENGTHS = {"VHH": (100, 150), "minibinder": (40, 100)}
MAX_ROWS = 20
COLUMNS = (
    "sequence_id", "target", "track", "molecule_type", "sequence", "length",
    "score", "score_method", "rank_in_track", "notes",
)
COMMON_TAGS = (
    "HHHHHH", "DYKDDDDK", "YPYDVPDYA", "EQKLISEEDL", "WSHPQFEK",
)


@dataclass(frozen=True)
class Track:
    track: int
    target: str
    molecule_type: str
    accession: str
    start: int
    end: int
    first: str
    last: str
    pdb: str
    chains: list[str]
    condition: str
    sha256: str

    @property
    def csv_name(self):
        return f"track{self.track}_{self.target}_{self.molecule_type}.csv"

    @property
    def fasta_name(self):
        return f"track{self.track}_{self.accession}_{self.start}-{self.end}.fasta"

    @property
    def fasta_id(self):
        return f"{self.accession}:{self.start}-{self.end}"


def read_tracks():
    return [Track(**item) for item in json.loads((TARGETS / "tracks.json").read_text())]


def read_fasta(text):
    header, *lines = text.splitlines()
    sequence = "".join(lines)
    if not header.startswith(">") or not sequence or set(sequence) - AMINO_ACIDS:
        raise ValueError("expected a single FASTA record with standard amino acids")
    return header[1:], sequence


def check_target_sequence(track, sequence):
    if len(sequence) != track.end - track.start + 1:
        raise ValueError(f"{track.fasta_id}: incorrect region length")
    if (sequence[0], sequence[-1]) != (track.first, track.last):
        raise ValueError(f"{track.fasta_id}: incorrect endpoint residues")
    if sha256(sequence.encode("ascii")).hexdigest() != track.sha256:
        raise ValueError(f"{track.fasta_id}: sequence differs from the pinned snapshot")


def check_cached_targets(tracks, directory):
    for track in tracks:
        path = directory / track.fasta_name
        header, sequence = read_fasta(path.read_text(encoding="ascii"))
        if header.split(maxsplit=1)[0] != track.fasta_id:
            raise ValueError(f"{path}: FASTA identifier does not match the track")
        check_target_sequence(track, sequence)


def fetch_targets(tracks, directory):
    records = {}
    for track in tracks:
        url = f"https://rest.uniprot.org/uniprotkb/{track.accession}.fasta"
        with urlopen(url, timeout=30) as response:
            header, full = read_fasta(response.read().decode("ascii"))
        fields = header.split("|", maxsplit=2)
        if len(fields) != 3 or fields[1] != track.accession or "OX=9606" not in fields[2].split():
            raise ValueError(f"{url}: unexpected UniProt accession or taxonomy")
        sequence = full[track.start - 1:track.end]
        check_target_sequence(track, sequence)
        records[track.fasta_name] = f">{track.fasta_id} {url}\n{fill(sequence, 60)}\n"

    directory.mkdir(parents=True, exist_ok=True)
    for name, record in records.items():
        with NamedTemporaryFile(
            mode="w", encoding="ascii", dir=directory, delete=False
        ) as handle:
            handle.write(record)
        Path(handle.name).replace(directory / name)


def write_templates(tracks, directory):
    directory.mkdir(parents=True, exist_ok=True)
    for track in tracks:
        path = directory / track.csv_name
        if path.exists():
            raise FileExistsError(f"{path} already exists")
    for track in tracks:
        with (directory / track.csv_name).open("x", encoding="utf-8", newline="") as handle:
            csv.writer(handle, lineterminator="\n").writerow(COLUMNS)


def validate_sequence(sequence, molecule_type):
    minimum, maximum = LENGTHS[molecule_type]
    if not minimum <= len(sequence) <= maximum:
        raise ValueError(f"{molecule_type} length must be {minimum}–{maximum} residues")
    invalid = sorted(set(sequence) - AMINO_ACIDS)
    if invalid:
        raise ValueError(f"invalid amino acid characters {invalid!r}")
    for tag in COMMON_TAGS:
        if sequence.startswith(tag) or sequence.endswith(tag):
            raise ValueError(f"common terminal tag detected ({tag})")


def validate_csvs(tracks, directory, allow_empty=False):
    sequence_ids = set()
    sequences = set()
    total = 0
    for track in tracks:
        path = directory / track.csv_name
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, strict=True)
            if tuple(next(reader, ())) != COLUMNS:
                raise ValueError(f"{path}: unexpected CSV header")

            count = 0
            for count, values in enumerate(reader, start=1):
                try:
                    if count > MAX_ROWS:
                        raise ValueError(f"track exceeds {MAX_ROWS} rows")
                    if len(values) != len(COLUMNS):
                        raise ValueError(f"expected {len(COLUMNS)} columns")
                    row = dict(zip(COLUMNS, values))
                    if (row["target"], row["track"], row["molecule_type"]) != (
                        track.target, str(track.track), track.molecule_type
                    ):
                        raise ValueError("target, track or molecule_type does not match filename")
                    sequence_id = row["sequence_id"]
                    if not sequence_id or sequence_id.strip() != sequence_id:
                        raise ValueError("missing or whitespace-padded sequence_id")
                    if sequence_id in sequence_ids:
                        raise ValueError(f"duplicate sequence_id {sequence_id!r}")

                    sequence = row["sequence"]
                    validate_sequence(sequence, track.molecule_type)
                    if sequence in sequences:
                        raise ValueError("duplicate sequence within or across tracks")
                    if row["length"] != str(len(sequence)):
                        raise ValueError("length differs from sequence length")
                    if row["rank_in_track"] != str(count):
                        raise ValueError("rank_in_track must match the 1-based row number")
                    label, _, method = row["score_method"].partition(":")
                    if label not in {"proxy", "official-method"} or not method.strip():
                        raise ValueError("score_method requires proxy: or official-method: and a name")
                    if not isfinite(float(row["score"])):
                        raise ValueError("score must be finite")
                except ValueError as error:
                    raise ValueError(f"{path}:{reader.line_num}: {error}") from error

                sequence_ids.add(sequence_id)
                sequences.add(sequence)
            if count == 0 and not allow_empty:
                raise ValueError(f"{path}: empty track CSV (use --allow-empty for templates)")
            total += count
    return total


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command, default, help_text in (
        ("check-targets", TARGETS, "verify the pinned FASTA cache offline"),
        ("fetch-targets", TARGETS, "restore the pinned FASTA cache from UniProt"),
        ("write-templates", SUBMISSIONS, "create seven header-only CSVs in a fresh directory"),
        ("validate", SUBMISSIONS, "check CSV headers and any populated rows"),
    ):
        subparser = commands.add_parser(command, help=help_text)
        subparser.add_argument("--directory", type=Path, default=default)
        if command == "validate":
            subparser.add_argument("--allow-empty", action="store_true")

    args = parser.parse_args(argv)
    try:
        tracks = read_tracks()
        if args.command == "check-targets":
            check_cached_targets(tracks, args.directory)
            print(f"Verified {len(tracks)} pinned UniProt regions.")
        elif args.command == "fetch-targets":
            fetch_targets(tracks, args.directory)
            print(f"Restored {len(tracks)} pinned UniProt regions.")
        elif args.command == "write-templates":
            write_templates(tracks, args.directory)
            print(f"Wrote {len(tracks)} header-only templates.")
        elif args.command == "validate":
            total = validate_csvs(tracks, args.directory, args.allow_empty)
            print(f"Checked {len(tracks)} CSVs and {total} rows.")
    except (ValueError, OSError, csv.Error) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
