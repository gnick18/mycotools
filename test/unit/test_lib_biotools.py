#! /usr/bin/env python3
"""Unit tests for the pure parsers in ``mycotools.lib.biotools``.

These are login-free, filesystem-light functions - ideal unit-test targets and a
worked example of how to grow the per-script scaffolds beyond CLI smoke tests.
"""
import pytest

from mycotools.lib.biotools import (
    fa2dict, dict2fa, gff2list, gff2features, list2gff, Feature, FeatureCol,
)


@pytest.fixture
def sample_fasta(tmp_path):
    fa = tmp_path / "sample.fa"
    fa.write_text(
        ">seq1 a description here\n"
        "ACGTACGT\n"
        "AAAA\n"  # multi-line sequence should be concatenated
        ">seq2\n"
        "TTTT\n"
    )
    return fa


def test_fa2dict_parses_headers_and_sequences(sample_fasta):
    parsed = fa2dict(str(sample_fasta))
    assert set(parsed) == {"seq1", "seq2"}
    assert parsed["seq1"]["sequence"] == "ACGTACGTAAAA"  # wrapped lines joined
    assert parsed["seq1"]["description"] == "a description here"
    assert parsed["seq2"]["sequence"] == "TTTT"
    assert parsed["seq2"]["description"] == ""  # no description


def test_dict2fa_emits_expected_format(sample_fasta):
    parsed = fa2dict(str(sample_fasta))
    out = dict2fa(parsed)
    assert ">seq1 a description here" in out
    assert ">seq2" in out
    assert "ACGTACGTAAAA" in out


def test_fa2dict_dict2fa_round_trip(sample_fasta, tmp_path):
    parsed = fa2dict(str(sample_fasta))
    roundtrip_path = tmp_path / "roundtrip.fa"
    roundtrip_path.write_text(dict2fa(parsed))
    assert fa2dict(str(roundtrip_path)) == parsed


# --------------------------------------------------------------------------- #
# TODO(scaffold): extend coverage of the other pure biotools parsers.
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# gff2list / list2gff
#
# The Feature migration made `gff2list` return a FeatureCol built without
# `ingest`, so it raised "No ID obtained for feature" on any GFF3, and the
# consumers (predb curation, `mtdb accession gff|locus|gbk`, 32 call sites) all
# index rows as dicts. Regenerating the 6,877-protein proteome of antflo1 from its
# GFF3 and FNA gave back the database's proteome exactly once this was repaired.
# --------------------------------------------------------------------------- #
GFF = (
    "##gff-version 3\n"
    "ctg1\tsrc\tregion\t1\t500\t.\t+\t.\tID=r1\n"
    "ctg1\tsrc\tgene\t11\t200\t.\t-\t.\tID=g1;Name=x\n"
    "ctg1\tsrc\tmRNA\t11\t200\t.\t-\t.\tID=m1;Parent=g1\n"
    "ctg1\tsrc\tCDS\t31\t180\t.\t-\t0\tID=c1;Parent=m1\n"
)
GFF_LINES = [l for l in GFF.split("\n") if l and not l.startswith("#")]


def test_gff2list_returns_column_dicts():
    rows = gff2list(GFF, path=False)
    assert len(rows) == 4 and all(isinstance(r, dict) for r in rows)
    assert rows[1] == {
        "seqid": "ctg1", "source": "src", "type": "gene", "start": 11, "end": 200,
        "score": ".", "strand": "-", "phase": ".", "attributes": "ID=g1;Name=x",
    }


def test_gff2list_rows_can_be_indexed_as_every_consumer_does():
    row = gff2list(GFF, path=False)[2]
    assert row["type"] == "mRNA" and "Parent=g1" in row["attributes"]


def test_gff2list_reads_a_file_and_stops_at_the_fasta_directive(tmp_path):
    path = tmp_path / "a.gff3"
    path.write_text(GFF + "##FASTA\n>ctg1\nACGT\n")
    assert len(gff2list(str(path))) == 4


def test_gff2list_reports_a_short_row_and_a_bad_coordinate():
    with pytest.raises(IndexError, match="expected tab-"):
        gff2list("ctg1\tsrc\tgene\t1\t2\n", path=False)
    with pytest.raises(ValueError, match="invalid integer"):
        gff2list("ctg1\tsrc\tgene\tx\t2\t.\t+\t.\tID=a\n", path=False)


def test_gff2list_list2gff_round_trip():
    text = list2gff(gff2list(GFF, path=False), ver=3)
    assert text == GFF.rstrip("\n")


def test_list2gff_accepts_features_as_well_as_rows():
    features = gff2features(GFF, path=False)
    assert list2gff(features, ver=0).split("\n") == GFF_LINES


def test_gff2features_derives_ids_parents_and_zero_based_coordinates():
    features = gff2features(GFF, path=False)
    assert isinstance(features, FeatureCol) and len(features) == 4
    gene = features[1]
    assert isinstance(gene, Feature)
    assert (gene.fid, gene.pid) == ("g1", None)
    assert (gene.start, gene.end) == (10, 200)  # GFF3 11..200, 0-based half-open
    assert features[3].pid == "m1" and features[3].phase == 0


@pytest.mark.skip(reason="TODO: scaffold - add gff2list / list2gff round-trip tests")
def test_gff2list_list2gff_round_trip(): ...


# ---- fa2dict gathers a sequence as a list of lines (it appended to a dictionary's string, quadratic in a contig) ----

def test_fa2dict_reads_records_descriptions_and_skips_comments(tmp_path):
    f = tmp_path / "g.fna"
    f.write_text(">c1 first contig\nACGT\nTT\n# a comment\nGG\n\n>c2\nCC\n>c3 with  two spaces\nA\n")
    from mycotools.lib.biotools import fa2dict
    d = fa2dict(str(f))
    assert list(d) == ["c1", "c2", "c3"]
    assert d["c1"] == {"sequence": "ACGTTTGG", "description": "first contig"}
    assert d["c2"] == {"sequence": "CC", "description": ""}
    assert d["c3"] == {"sequence": "A", "description": "with  two spaces"}


def test_fa2dict_keeps_an_empty_record_and_the_last_of_a_repeated_name(tmp_path):
    f = tmp_path / "g.fna"
    f.write_text(">e\n>r\nAAA\n>r\nCC\n")
    from mycotools.lib.biotools import fa2dict
    d = fa2dict(str(f))
    assert d["e"]["sequence"] == "" and d["r"]["sequence"] == "CC"


def test_fa2dict_refuses_sequence_before_a_header(tmp_path):
    import pytest
    from mycotools.lib.biotools import fa2dict
    f = tmp_path / "g.fna"
    f.write_text("ACGT\n>c1\nAA\n")
    with pytest.raises(UnboundLocalError):
        fa2dict(str(f))


def test_fa2dict_str_matches_fa2dict():
    from mycotools.lib.biotools import fa2dict_str
    d = fa2dict_str(">c1 x\nAC\nGT\n>c2\nTT")
    assert d == {"c1": {"sequence": "ACGT", "description": "x"}, "c2": {"sequence": "TT", "description": ""}}


def test_fa2dict_reads_a_long_contig_in_linear_time(tmp_path):
    # 8 Mb in 100,000 lines of 80: the quadratic version needs about ten seconds for this (hours for a 346 Mb genome)
    import time
    from mycotools.lib.biotools import fa2dict
    f = tmp_path / "long.fna"
    f.write_text(">chr1 long\n" + ("ACGT" * 20 + "\n") * 100000)
    t = time.time()
    d = fa2dict(str(f))
    assert len(d["chr1"]["sequence"]) == 8_000_000
    assert time.time() - t < 5
