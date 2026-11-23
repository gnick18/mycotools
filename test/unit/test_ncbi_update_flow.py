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
