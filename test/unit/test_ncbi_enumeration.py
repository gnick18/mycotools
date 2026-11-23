#! /usr/bin/env python3
"""Offline tests for listing NCBI genomes from the Datasets API.

An update used to enumerate `GENOME_REPORTS/eukaryotes.txt`, which NCBI stopped
regenerating on 2024-10-10 (its newest genome was released 2024-08-30), so no
genome deposited afterwards could be found and nothing said so. On 2026-11-22
NCBI listed 6,167 annotated current GenBank fungal assemblies, 815 of them
released after the 2025-04-30 database, and 803 of those were invisible to that
file. These tests pin the replacement:

* a report becomes the row the rest of the update expects;
* the listing is paged, retried when the API is busy, and never returned short;
* it is cached only once complete;
* an assembly already in the database is an update only when NCBI holds a newer
  VERSION of it. Compared by release date, 711 of the 4,254 genomes of the
  2025-04-30 database that NCBI still lists came out as updated.
"""
import pandas as pd
import pytest
import requests

from mycotools.mtdb import update as mod
from mycotools.utils.ncbi2db import accession_version, redundancy_check


def report(accession, release="2025-07-29", annotated="2025-07-29", genes=15310,
           proteins=15024, biosample="SAMN52443102"):
    return {
        "accession": accession,
        "organism": {"organism_name": "Neocosmospora sp.", "tax_id": 3469074},
        "assembly_info": {
            "assembly_level": "Scaffold",
            "bioproject_accession": "PRJNA1255138",
            "biosample": {"accession": biosample},
            "release_date": release,
        },
        "annotation_info": {
            "release_date": annotated,
            "stats": {"gene_counts": {"total": genes, "protein_coding": proteins}},
        },
    }


class Response:
    def __init__(self, payload=None, status=200):
        self.payload, self.status_code = payload, status

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def json(self):
        if self.payload is None:
            raise ValueError("not json")
        return self.payload


@pytest.fixture
def api(monkeypatch):
    """Replace `requests.get` with a scripted queue, recording each call."""
    calls, queue = [], []

    def fake_get(url, params=None, headers=None, timeout=None):
        calls.append({"url": url, "params": dict(params), "headers": headers})
        return queue.pop(0)

    monkeypatch.setattr(mod.requests, "get", fake_get)
    monkeypatch.setattr(mod.time, "sleep", lambda s: None)
    return calls, queue


# --------------------------------------------------------------------------- #
# a report as a row of the old table


def test_report_row_has_every_column_the_update_reads():
    row = mod.ncbi_report_row(report("GCA_057674295.1"), "Fungi")
    assert list(row) == mod.NCBI_REPORT_COLUMNS
    assert row["Assembly Accession"] == "GCA_057674295.1"
    assert row["Group"] == "Fungi"
    assert row["BioSample Accession"] == "SAMN52443102"
    assert (row["Genes"], row["Proteins"]) == (15310, 15024)


def test_report_row_dates_use_slashes_and_modify_is_the_later_release():
    row = mod.ncbi_report_row(
        report("GCA_1.1", release="2019-02-20", annotated="2023-11-22"), "Fungi"
    )
    assert row["Release Date"] == "2019/02/20"
    assert row["Modify Date"] == "2023/11/22"


def test_report_row_without_annotation_counts_is_dropped_by_clean_ncbi_df():
    bare = report("GCA_2.1")
    del bare["annotation_info"]
    row = mod.ncbi_report_row(bare, "Fungi")
    assert row["Genes"] == "-" and row["Proteins"] == "-"


# --------------------------------------------------------------------------- #
# paging, retries, never short


def test_listing_follows_page_tokens(api):
    calls, queue = api
    queue.append(Response({"total_count": 3, "reports": [report("GCA_1.1"), report("GCA_2.1")],
                           "next_page_token": "tok"}))
    queue.append(Response({"total_count": 3, "reports": [report("GCA_3.1")]}))
    reports = mod.fetch_ncbi_reports(4751)
    assert [r["accession"] for r in reports] == ["GCA_1.1", "GCA_2.1", "GCA_3.1"]
    assert "page_token" not in calls[0]["params"]
    assert calls[1]["params"]["page_token"] == "tok"
    assert calls[0]["url"].endswith("/genome/taxon/4751/dataset_report")
    assert calls[0]["params"]["filters.assembly_source"] == "genbank"
    assert calls[0]["params"]["filters.has_annotation"] == "true"


def test_api_key_travels_as_a_header(api):
    calls, queue = api
    queue.append(Response({"total_count": 0, "reports": []}))
    mod.fetch_ncbi_reports(4751, api="abc123")
    assert calls[0]["headers"] == {"api-key": "abc123"}


def test_a_busy_api_is_retried(api):
    calls, queue = api
    queue.extend([Response(status=503), Response(status=429),
                  Response({"total_count": 1, "reports": [report("GCA_1.1")]})])
    assert len(mod.fetch_ncbi_reports(4751)) == 1
    assert len(calls) == 3


def test_a_bad_request_is_not_retried(api):
    calls, queue = api
    queue.append(Response(status=400))
    with pytest.raises(RuntimeError):
        mod.fetch_ncbi_reports(4751)
    assert len(calls) == 1


def test_an_api_that_never_answers_raises(api):
    calls, queue = api
    queue.extend([Response(status=503)] * 4)
    with pytest.raises(RuntimeError, match="did not answer"):
        mod.fetch_ncbi_reports(4751, max_attempts=4)
    assert len(calls) == 4


def test_a_short_listing_raises_instead_of_reading_as_nothing_new(api):
    calls, queue = api
    queue.append(Response({"total_count": 5, "reports": [report("GCA_1.1")]}))
    with pytest.raises(RuntimeError, match="listed 5"):
        mod.fetch_ncbi_reports(4751)


# --------------------------------------------------------------------------- #
# the cached table


def test_listing_is_written_once_complete_and_reused(api, tmp_path):
    calls, queue = api
    queue.append(Response({"total_count": 2, "reports": [report("GCA_1.1"), report("GCA_2.1")]}))
    path = str(tmp_path / "20261122.ncbi.tsv")
    df = mod.dwnld_ncbi_metadata(path, kingdom="fungi")
    assert list(df["Assembly Accession"]) == ["GCA_1.1", "GCA_2.1"]
    assert set(df["Group"]) == {"Fungi"}
    again = mod.dwnld_ncbi_metadata(path, kingdom="fungi")  # no further request
    assert len(calls) == 1 and len(again) == 2


def test_the_listing_survives_the_rest_of_the_update_unchanged(api, tmp_path):
    """The table must go through `clean_ncbi_df`'s own selection."""
    calls, queue = api
    queue.append(Response({"total_count": 2, "reports": [report("GCA_1.1"), report("GCA_2.1")]}))
    df = mod.dwnld_ncbi_metadata(str(tmp_path / "x.tsv"), kingdom="fungi").astype(str)
    df = df[df["Group"] == "Fungi"]
    df = df[~df["Genes"].isin({"-", "", "0"})]
    assert len(df) == 2
    assert pd.to_datetime(df["Modify Date"]).notna().all()


def test_a_failed_listing_leaves_no_cache_and_exits(api, tmp_path):
    calls, queue = api
    queue.extend([Response(status=503)] * 4)
    path = tmp_path / "x.tsv"
    with pytest.raises(SystemExit) as stop:
        mod.dwnld_ncbi_metadata(str(path), kingdom="fungi")
    assert stop.value.code == mod.NCBI_REPORT_EXIT
    assert not path.exists() and not (tmp_path / "x.tsv.tmp").exists()


def test_an_empty_listing_exits_rather_than_report_no_updates(api, tmp_path):
    calls, queue = api
    queue.append(Response({"total_count": 0, "reports": []}))
    with pytest.raises(SystemExit) as stop:
        mod.dwnld_ncbi_metadata(str(tmp_path / "x.tsv"), kingdom="fungi")
    assert stop.value.code == mod.NCBI_REPORT_EXIT


def test_each_kingdom_has_a_taxon_and_the_group_clean_ncbi_df_selects_on():
    assert mod.NCBI_KINGDOM_TAXA["fungi"] == (4751, "Fungi")
    for kingdom, group in {"plants": "Plants", "animals": "Animals", "archaea": "Archaea"}.items():
        assert mod.NCBI_KINGDOM_TAXA[kingdom][1] == group
    assert mod.NCBI_KINGDOM_TAXA["metazoa"] == mod.NCBI_KINGDOM_TAXA["animals"]


def test_an_unknown_kingdom_is_refused(tmp_path):
    with pytest.raises(KeyError):
        mod.dwnld_ncbi_metadata(str(tmp_path / "x.tsv"), kingdom="protists")


# --------------------------------------------------------------------------- #
# an assembly in the database is an update only when its VERSION is newer


def database(rows):
    return pd.DataFrame(
        [{"ome": ome, "genus": "Aspergillus", "species": "niger", "strain": "X",
          "assembly_acc": acc, "version": version, "published": "2020/01/29"}
         for ome, acc, version in rows]
    )


def listing(accessions, date="2026-01-01"):
    return pd.DataFrame(
        [{"assembly_acc": acc, "genus": "Aspergillus", "species": "niger",
          "strain": "X", "version": pd.Timestamp(date)} for acc in accessions]
    )


def test_accession_version():
    assert accession_version("GCA_000001405.29") == 29
    assert accession_version("GCA_000001405.1") == 1
    assert accession_version("SRR123") == 0
    assert accession_version("Not available") == 0


def test_same_version_with_a_later_release_date_is_not_an_update():
    """The case that re-versioned 711 genomes when compared by date: the record
    was released after the date the database stored for it."""
    db = database([("aspnig1", "GCA_910592055.1", "20211019000000")])
    new, kept, updates, old = redundancy_check(db, listing(["GCA_910592055.1"], "2021-11-12"), "assembly_acc")
    assert updates == [] and len(new) == 0 and list(kept["ome"]) == ["aspnig1"]


def test_a_newer_version_of_a_database_assembly_is_an_update():
    db = database([("aspnig1", "GCA_000001.1", "20200101000000")])
    new, kept, updates, old = redundancy_check(db, listing(["GCA_000001.2"]), "assembly_acc")
    assert [u[0] for u in updates] == ["GCA_000001.1"] and [u[3] for u in updates] == ["GCA_000001.2"]
    assert len(new) == 1 and "aspnig1" not in set(kept["ome"])


def test_an_older_version_is_not_an_update():
    db = database([("aspnig1", "GCA_000001.2", "20200101000000")])
    new, kept, updates, old = redundancy_check(db, listing(["GCA_000001.1"], "2030-01-01"), "assembly_acc")
    assert updates == [] and len(new) == 0


def test_an_unknown_assembly_is_new():
    db = database([("aspnig1", "GCA_000001.1", "20200101000000")])
    new, kept, updates, old = redundancy_check(db, listing(["GCA_999999.1"]), "assembly_acc")
    assert [u[:2] for u in updates] == [[None, None]] and len(new) == 1


def test_a_database_version_that_is_not_a_date_no_longer_blocks_a_version_update():
    db = database([("aspnig1", "GCA_000001.1", "1")])
    new, kept, updates, old = redundancy_check(db, listing(["GCA_000001.2"]), "assembly_acc")
    assert len(updates) == 1
