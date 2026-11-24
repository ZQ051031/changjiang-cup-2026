import csv
import json
from urllib.error import URLError

import pytest

import cup


@pytest.fixture
def tracks():
    return cup.read_tracks()


@pytest.fixture
def templates(tracks, tmp_path):
    cup.write_templates(tracks, tmp_path)
    return tmp_path


def fixture_row(track, *, residue="A", rank=1, sequence_id=None):
    """Use homopolymer-based data solely to exercise the CSV format."""
    sequence = residue * cup.LENGTHS[track.molecule_type][0]
    return {
        "sequence_id": sequence_id or f"fixture-{track.track}-{rank}",
        "target": track.target,
        "track": track.track,
        "molecule_type": track.molecule_type,
        "sequence": sequence,
        "length": len(sequence),
        "score": 0,
        "score_method": "proxy:format test fixture",
        "rank_in_track": rank,
        "notes": 'Test data, with "quoted" punctuation.',
    }


def append_rows(path, rows):
    with path.open("a", newline="", encoding="utf-8") as handle:
        csv.DictWriter(handle, fieldnames=cup.COLUMNS).writerows(rows)


def test_track_metadata(tracks):
    assert [
        (item.track, item.target, item.molecule_type, item.accession,
         item.start, item.end, item.first, item.last, item.pdb, item.chains)
        for item in tracks
    ] == [
        (1, "EGFR", "VHH", "P00533-1", 1, 645, "M", "S", "6ARU", ["A"]),
        (2, "TNF-alpha", "VHH", "P01375", 77, 233, "V", "L", "1TNF", ["A", "B", "C"]),
        (3, "GLP-1R", "VHH", "P43220", 24, 432, "R", "S", "6X18", ["R"]),
        (4, "TfR1", "VHH", "P02786", 89, 760, "C", "F", "11QC", ["A", "B"]),
        (5, "TL1A", "minibinder", "O95150", 72, 251, "L", "L", "2RE9", ["A", "B", "C"]),
        (6, "TNFR1", "VHH", "P19438-1", 30, 211, "L", "T", "9S3I", ["A"]),
        (7, "IL-7Ralpha", "minibinder", "P16871-1", 21, 239, "E", "D", "3DI3", ["B"]),
    ]
    assert "pH 6.0" in tracks[0].condition
    assert tracks[2].condition == "Agonist required."


def test_pinned_cache_is_valid(tracks):
    cup.check_cached_targets(tracks, cup.TARGETS)


def test_changed_region_fails_digest_check(tracks):
    track = tracks[0]
    path = cup.TARGETS / track.fasta_name
    _, sequence = cup.read_fasta(path.read_text(encoding="ascii"))
    replacement = "A" if sequence[1] != "A" else "C"
    with pytest.raises(ValueError, match="pinned snapshot"):
        cup.check_target_sequence(track, sequence[0] + replacement + sequence[2:])


def test_network_failure_does_not_replace_cache(tracks, tmp_path, monkeypatch):
    path = tmp_path / tracks[0].fasta_name
    path.write_bytes(b"existing cache")

    def unavailable(*_args, **_kwargs):
        raise URLError("offline")

    monkeypatch.setattr(cup, "urlopen", unavailable)
    with pytest.raises(URLError, match="offline"):
        cup.fetch_targets(tracks, tmp_path)
    assert path.read_bytes() == b"existing cache"
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize(
    ("molecule_type", "length"),
    [("VHH", 100), ("VHH", 150), ("minibinder", 40), ("minibinder", 100)],
)
def test_accepts_length_boundaries(molecule_type, length):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    cup.validate_sequence((alphabet * 8)[:length], molecule_type)


@pytest.mark.parametrize(
    ("molecule_type", "length"),
    [("VHH", 99), ("VHH", 151), ("minibinder", 39), ("minibinder", 101)],
)
def test_rejects_out_of_bounds_length(molecule_type, length):
    with pytest.raises(ValueError, match="length must be"):
        cup.validate_sequence("A" * length, molecule_type)


@pytest.mark.parametrize("character", ["X", "U", "a", "-", " ", "*"])
def test_rejects_nonstandard_residue(character):
    with pytest.raises(ValueError, match="invalid amino acid"):
        cup.validate_sequence(character + "A" * 99, "VHH")


@pytest.mark.parametrize("sequence", ["HHHHHH" + "A" * 94, "A" * 92 + "DYKDDDDK"])
def test_rejects_common_terminal_tags(sequence):
    with pytest.raises(ValueError, match="terminal tag"):
        cup.validate_sequence(sequence, "VHH")


def test_templates_match_contract(tracks, templates):
    expected = {
        "track1_EGFR_VHH.csv", "track2_TNF-alpha_VHH.csv", "track3_GLP-1R_VHH.csv",
        "track4_TfR1_VHH.csv", "track5_TL1A_minibinder.csv", "track6_TNFR1_VHH.csv",
        "track7_IL-7Ralpha_minibinder.csv",
    }
    assert {path.name for path in templates.glob("*.csv")} == expected
    header = (
        "sequence_id,target,track,molecule_type,sequence,length,score,"
        "score_method,rank_in_track,notes\n"
    )
    for name in expected:
        assert (templates / name).read_text() == header
        assert (cup.SUBMISSIONS / name).read_text() == header
    assert cup.validate_csvs(tracks, templates, allow_empty=True) == 0
    with pytest.raises(ValueError, match="empty track CSV"):
        cup.validate_csvs(tracks, templates)


def test_templates_preserve_existing_files(tracks, templates):
    path = templates / tracks[0].csv_name
    content = path.read_bytes()
    with pytest.raises(FileExistsError):
        cup.write_templates(tracks, templates)
    assert path.read_bytes() == content


def test_population_across_all_tracks(tracks, templates):
    for track, residue in zip(tracks, "ACDEFGI", strict=True):
        append_rows(templates / track.csv_name, [fixture_row(track, residue=residue)])
    assert cup.validate_csvs(tracks, templates) == 7


def test_accepts_spreadsheet_csv_with_utf8_bom(tracks, templates):
    track = tracks[0]
    path = templates / track.csv_name
    append_rows(path, [fixture_row(track)])
    original = path.read_text(encoding="utf-8")
    path.write_text(original, encoding="utf-8-sig")
    assert cup.validate_csvs(tracks, templates, allow_empty=True) == 1


@pytest.mark.parametrize("count", [20, 21])
def test_enforces_candidate_cap(tracks, templates, count):
    track = tracks[0]
    rows = [fixture_row(track, rank=rank) for rank in range(1, count + 1)]
    for rank, row in enumerate(rows, start=1):
        row["sequence"] = "C" * rank + "A" * (100 - rank)
    append_rows(templates / track.csv_name, rows)
    if count == 20:
        assert cup.validate_csvs(tracks, templates, allow_empty=True) == 20
    else:
        with pytest.raises(ValueError, match="exceeds 20 rows"):
            cup.validate_csvs(tracks, templates, allow_empty=True)


@pytest.mark.parametrize("cross_track", [False, True])
def test_rejects_duplicate_sequences(tracks, templates, cross_track):
    first, second = tracks[:2]
    append_rows(templates / first.csv_name, [fixture_row(first)])
    other = second if cross_track else first
    rank = 1 if cross_track else 2
    append_rows(templates / other.csv_name, [fixture_row(other, rank=rank)])
    with pytest.raises(ValueError, match="duplicate sequence within or across tracks"):
        cup.validate_csvs(tracks, templates, allow_empty=True)


@pytest.mark.parametrize(
    ("field", "value", "error"),
    [
        ("target", "Wrong target", "does not match filename"),
        ("track", 2, "does not match filename"),
        ("molecule_type", "minibinder", "does not match filename"),
        ("sequence_id", " ", "whitespace-padded sequence_id"),
        ("length", 101, "length differs"),
        ("rank_in_track", 2, "1-based row number"),
        ("score", "nan", "score must be finite"),
        ("score_method", "ipSAE", "requires proxy"),
    ],
)
def test_rejects_bad_fields(tracks, templates, field, value, error):
    track = tracks[0]
    row = fixture_row(track)
    row[field] = value
    append_rows(templates / track.csv_name, [row])
    with pytest.raises(ValueError, match=error):
        cup.validate_csvs(tracks, templates, allow_empty=True)


def test_rejects_duplicate_ids_across_tracks(tracks, templates):
    for track, residue in zip(tracks[:2], "AC", strict=True):
        row = fixture_row(track, residue=residue, sequence_id="duplicate")
        append_rows(templates / track.csv_name, [row])
    with pytest.raises(ValueError, match="duplicate sequence_id"):
        cup.validate_csvs(tracks, templates, allow_empty=True)


@pytest.mark.parametrize(
    ("bad_csv", "error"),
    [("bad-header\n", "unexpected CSV header"), ("one,two\n", "expected 10 columns")],
)
def test_rejects_bad_csv_shape(tracks, templates, bad_csv, error):
    path = templates / tracks[0].csv_name
    mode = "w" if bad_csv.startswith("bad-header") else "a"
    with path.open(mode, encoding="utf-8") as handle:
        handle.write(bad_csv)
    with pytest.raises(ValueError, match=error):
        cup.validate_csvs(tracks, templates, allow_empty=True)


def test_cli_reports_bad_csv(templates, capsys):
    assert cup.main(["validate", "--directory", str(templates)]) == 1
    assert "empty track CSV" in capsys.readouterr().err


def test_describe_fasta_reports_counts_and_mean_hydropathy(tmp_path):
    path = tmp_path / "example.fasta"
    path.write_text(">example composition fixture\nACCDE\nHIKR\n")
    assert cup.describe_fasta(path) == {
        "id": "example",
        "length": 9,
        "cysteine_count": 2,
        "acidic_residue_count": 2,
        "basic_residue_count": 2,
        "basic_minus_acidic_count": 0,
        "mean_kyte_doolittle_hydropathy": pytest.approx(-7.3 / 9),
    }


def test_describe_cli_preserves_input_order_and_is_deterministic(tmp_path, capsys, monkeypatch):
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    first.write_text(">first\nKRR\n")
    second.write_text(">second\nDEE\n", encoding="utf-8-sig")
    # Description of standalone FASTAs should work without any track data.
    monkeypatch.setattr(cup, "TARGETS", tmp_path / "missing")
    args = ["describe", str(second), str(first)]
    assert cup.main(args) == 0
    original = capsys.readouterr().out
    reports = json.loads(original)
    assert [report["id"] for report in reports] == ["second", "first"]
    assert [report["basic_minus_acidic_count"] for report in reports] == [-3, 3]
    assert cup.main(args) == 0
    assert capsys.readouterr().out == original


@pytest.mark.parametrize("contents", ["", ">\nAA\n", ">id\n", ">id\nAX\n", ">a\nA\n>b\nA\n"])
def test_describe_rejects_invalid_single_record_fasta(tmp_path, capsys, contents):
    good = tmp_path / "good.fasta"
    good.write_text(">good\nAA\n")
    bad = tmp_path / "bad.fasta"
    bad.write_text(contents)
    assert cup.main(["describe", str(good), str(bad)]) == 1
    output = capsys.readouterr()
    assert output.out == ""
    assert f"{bad}: expected a single FASTA record" in output.err
