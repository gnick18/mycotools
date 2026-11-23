#! /usr/bin/env python3
"""The NCBI half of an update, from the redundancy check to the omes the genomes
are given, with the download itself stubbed.

* An updated genome is carried forward under its OLD ome code: the old row is
  taken out of the reference database before the omes are generated, so
  `gen_omes` has no reference row to bump (aspnig1 stays aspnig1, aspnig7.1
  stays aspnig7.1). It replaces its predecessor in place.
* When the download of an updated genome fails, its old row must come back.
  `ref_db.append(...)` raised AttributeError under pandas 2, so one failed
  download stopped the whole update.
"""
import pandas as pd
import pytest

from mycotools.lib.dbtools import mtdb
from mycotools.mtdb import predb as predb_mod
from mycotools.utils import ncbi2db

COLUMNS = ["ome", "taxonomy", "genus", "species", "strain", "version", "source",
           "biosample", "assembly_acc", "acquisition_date", "published", "fna",
           "faa", "gff3"]


def database(rows):
    return pd.DataFrame(
        [{**{c: "" for c in COLUMNS}, "genus": "Aspergillus", "species": "niger",
          "strain": "X", "source": "ncbi", "published": "2020/01/29", "ome": ome,
          "assembly_acc": acc, "version": version} for ome, acc, version in rows]
    )


def listing(accessions):
    return pd.DataFrame(
        [{"assembly_acc": acc, "genus": "Aspergillus", "species": "niger", "strain": "X",
          "version": pd.Timestamp("2026-01-01"), "Release Date": "2026/01/01",
          "BioSample Accession": "SAMN1"} for acc in accessions]
    )


REF = [("aspnig1", "GCA_000001.1", "20200101000000"),
       ("aspnig2", "GCA_000002.1", "20200101000000"),
       ("aspfis2", "GCA_000003", "4"),
       ("aspnig7.1", "GCA_000007.2", "20200101000000")]


def fetched_all(**kw):
    df = kw["ncbi_df"].copy()
    df["assemblyPath"], df["gffPath"] = "x.fna", "x.gff"
    return df, []


@pytest.fixture
def run(tmp_path, monkeypatch):
    def _run(accessions, downloader=fetched_all):
        monkeypatch.setattr(ncbi2db, "ncbi_dwnld", downloader)
        monkeypatch.chdir(tmp_path)
        return ncbi2db.main(str(tmp_path) + "/", ncbi_df=listing(accessions),
                            ref_db=database(REF), date="20261123")
    return _run


def test_an_updated_genome_keeps_its_ome_and_the_unchanged_ones_are_not_sent(run):
    ncbi_df, ref_db, failed = run(
        ["GCA_000001.2", "GCA_000007.3", "GCA_000002.1", "GCA_000003.4", "GCA_000009.1"])
    assert list(ncbi_df["assembly_acc"]) == ["GCA_000001.2", "GCA_000007.3", "GCA_000009.1"]
    assert list(ncbi_df["ome"].fillna("")) == ["aspnig1", "aspnig7.1", ""]
    # the replaced rows are out of the reference database, the rest stay
    assert sorted(ref_db["ome"]) == ["aspfis2", "aspnig2"]
    assert failed == []


def test_the_omes_generated_for_an_update_are_the_old_ones(run):
    ncbi_df, ref_db, failed = run(["GCA_000001.2", "GCA_000007.3", "GCA_000009.1"])
    predb = ncbi_df.fillna("").copy()
    predb["version"] = predb["version"].astype(str)
    omedb, bad = predb_mod.gen_omes(
        predb_mod.predb2mtdb(predb.to_dict(orient="list")),
        mtdb.pd2mtdb(ref_db), spacer="")
    got = dict(zip(omedb["assembly_acc"], omedb["ome"]))
    assert got["GCA_000001.2"] == "aspnig1"
    assert got["GCA_000007.3"] == "aspnig7.1"
    assert got["GCA_000009.1"] not in {"aspnig1", "aspnig2", "aspnig7.1", "aspfis2"}


def partial(**kw):
    df = kw["ncbi_df"].copy()
    df = df[df["assembly_acc"] != "GCA_000001.2"].copy()
    df["assemblyPath"], df["gffPath"] = "x.fna", "x.gff"
    return df, [["GCA_000001.2", "20260101"]]


def test_a_failed_download_of_an_updated_genome_puts_the_old_row_back(run):
    ncbi_df, ref_db, failed = run(["GCA_000001.2", "GCA_000009.1"], downloader=partial)
    assert failed == [["GCA_000001.2", "20260101"]]
    assert list(ncbi_df["assembly_acc"]) == ["GCA_000009.1"]
    assert "ome" not in ncbi_df or not ncbi_df["ome"].dropna().tolist()
    assert sorted(ref_db["ome"]) == ["aspfis2", "aspnig1", "aspnig2", "aspnig7.1"]
    restored = ref_db.set_index("ome").loc["aspnig1"]
    assert restored["assembly_acc"] == "GCA_000001.1"
    assert restored["genus"] == "Aspergillus" and restored["version"] == "20200101000000"
    assert restored["versionless"] == "GCA_000001" and "index" not in ref_db.columns


def test_a_failed_new_genome_leaves_the_reference_alone(run):
    def drop_new(**kw):
        df = kw["ncbi_df"].copy()
        df = df[df["assembly_acc"] != "GCA_000009.1"].copy()
        df["assemblyPath"], df["gffPath"] = "x.fna", "x.gff"
        return df, [["GCA_000009.1", "20260101"]]
    ncbi_df, ref_db, failed = run(["GCA_000001.2", "GCA_000009.1"], downloader=drop_new)
    assert list(ncbi_df["ome"].dropna()) == ["aspnig1"]
    assert sorted(ref_db["ome"]) == ["aspfis2", "aspnig2", "aspnig7.1"]


def test_an_update_the_downloader_neither_returned_nor_reported_adds_no_row(run):
    def vanish(**kw):
        df = kw["ncbi_df"].copy()
        df = df[df["assembly_acc"] != "GCA_000001.2"].copy()
        df["assemblyPath"], df["gffPath"] = "x.fna", "x.gff"
        return df, []
    ncbi_df, ref_db, failed = run(["GCA_000001.2", "GCA_000009.1"], downloader=vanish)
    assert list(ncbi_df["assembly_acc"]) == ["GCA_000009.1"]


# --------------------------------------------------------------------------- #
# the database is de-duplicated before anything else in an update

def test_distinct_accessions_without_a_version_are_not_taken_for_versions_of_one_genome():
    """no_biosample_yet_8954-AG-0001 and -0007 are two isolates. Cut at the first
    "." both became no_biosample_yet_8954-AG-000 and one was deleted; seven genomes
    of the 2025-04-30 database went that way."""
    from mycotools.mtdb import update as update_mod
    rows = [(f"asplat{n}", f"no_biosample_yet_8954-AG-00{n:02d}", "1") for n in (1, 7, 8, 13, 14, 15, 17, 18, 19)]
    kept = update_mod.internal_redundancy_check(database(rows))
    assert sorted(kept["ome"]) == sorted(o for o, _, _ in rows)


def test_two_versions_of_one_assembly_in_the_database_collapse_to_the_newest():
    from mycotools.mtdb import update as update_mod
    rows = [("aspnig1", "GCA_000001.1", "20200101000000"), ("aspnig2", "GCA_000001.2", "20200101000000"),
            ("aspnig3", "GCA_000002", "4"), ("aspnig4", "GCA_000003.1", "20200101000000")]
    kept = update_mod.internal_redundancy_check(database(rows))
    assert sorted(kept["ome"]) == ["aspnig2", "aspnig3", "aspnig4"]


# --------------------------------------------------------------------------- #
# a separately curated unversioned copy of an assembly is never deleted by an update
# (paevar10 = GCA_030061375.1 and paevar46 = the same assembly and biosample under
# GCA_030061375, with its own annotation: 26 such pairs sit in the 2025-04-30 database)

TWINS = [("paevar10", "GCA_030061375.1", "20230522000000"), ("paevar46", "GCA_030061375", "1")]


def test_the_ncbi_check_leaves_an_unversioned_twin_alone_when_nothing_is_newer():
    new, kept, updates, old = ncbi2db.redundancy_check(database(TWINS), listing(["GCA_030061375.1"]), "assembly_acc")
    assert updates == [] and len(new) == 0
    assert sorted(kept["ome"]) == ["paevar10", "paevar46"]


def test_the_ncbi_check_replaces_only_the_row_that_is_updated():
    new, kept, updates, old = ncbi2db.redundancy_check(database(TWINS), listing(["GCA_030061375.2"]), "assembly_acc")
    assert [u[3] for u in updates] == ["GCA_030061375.2"] and len(new) == 1
    assert len(updates) == 1 and updates[0][4] == "paevar10"
    assert sorted(kept["ome"]) == ["paevar46"]


def test_the_ncbi_check_still_drops_a_superseded_versioned_row():
    rows = [("aspnig1", "GCA_000001.1", "20200101000000"), ("aspnig2", "GCA_000001.2", "20200101000000")]
    new, kept, updates, old = ncbi2db.redundancy_check(database(rows), listing(["GCA_000001.2"]), "assembly_acc")
    assert updates == [] and list(kept["ome"]) == ["aspnig2"]


def test_the_database_check_keeps_unversioned_twins_and_drops_only_a_superseded_version():
    from mycotools.mtdb import update as update_mod
    kept = update_mod.internal_redundancy_check(database(TWINS))
    assert sorted(kept["ome"]) == ["paevar10", "paevar46"]
    rows = [("aspnig1", "GCA_000001.1", "20200101000000"), ("aspnig2", "GCA_000001.2", "20200101000000"),
            ("aspnig3", "GCA_000001", "1")]
    kept = update_mod.internal_redundancy_check(database(rows))
    assert sorted(kept["ome"]) == ["aspnig2", "aspnig3"]  # the older versioned copy goes, the twin stays


def test_the_database_check_collapses_an_exact_repeat_of_one_accession():
    from mycotools.mtdb import update as update_mod
    rows = [("aspnig1", "GCA_000001.1", "20200101000000"), ("aspnig2", "GCA_000001.1", "20200101000000")]
    kept = update_mod.internal_redundancy_check(database(rows))
    assert list(kept["ome"]) == ["aspnig2"]
