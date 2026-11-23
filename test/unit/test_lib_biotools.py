#! /usr/bin/env python3
"""Unit tests for the pure parsers in ``mycotools.lib.biotools``.

These are login-free, filesystem-light functions - ideal unit-test targets and a
worked example of how to grow the per-script scaffolds beyond CLI smoke tests.
"""
import pytest

from mycotools.lib.biotools import fa2dict, dict2fa


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
