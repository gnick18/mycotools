#! /usr/bin/env python3
"""The MycoCosm half of an update: a portal already in the database whose catalog
version moved on, from the redundancy check to the rows that survive.

* The database stores a portal's version with the dot removed (1.0 -> "10",
  2.0 -> "20", 3.16 -> "316", none -> "0"; 1,536 of 1,545 MycoCosm genomes of the
  2025-04-30 database), the catalog reports "1.0". Compared as numbers, "2.0"
  was never above "10", so a portal moving from 1.0 to 2.0 was never an update
  (958 of the 1,536 genomes could not have been flagged by any catalog change).
* An update that cannot be fetched (failed, or still waiting on JGI's tape
  restore, which is where nearly every genome starts) must leave the genome it
  was to replace in the database. The old row was filed under the database's
  ome code but looked up by portal, so it was never found, and `DataFrame.append`
  no longer exists under pandas 2.
"""
import pandas as pd
import pytest

from mycotools.utils import jgi2db
from mycotools.utils.jgi2db import jgi_redundancy_check, parse_version

COLUMNS = ["ome", "taxonomy", "genus", "species", "strain", "version", "source",
           "biosample", "assembly_acc", "acquisition_date", "published", "fna",
           "faa", "gff3"]


def database(rows):
    """rows: (ome, portal, stored version, source)"""
    return pd.DataFrame(
        [{**{c: "" for c in COLUMNS}, "genus": "Aspergillus", "species": "niger",
          "strain": "X", "published": "1", "ome": ome, "assembly_acc": portal,
          "version": version, "source": source} for ome, portal, version, source in rows]
    )


def catalog(rows):
    """rows: (portal, catalog version)"""
    return pd.DataFrame(
        [{"portal": portal, "name": f"Aspergillus niger {portal} v{version}",
          "version": version, "biosample": "", "is published": "Y", "is public": "Y",
          "publication(s)": "Doe 2024", "genus": "Aspergillus", "species": "niger",
          "strain": "X"} for portal, version in rows]
    )


def check(db_rows, cat_rows):
    db = database(db_rows)
    db["index"] = db["assembly_acc"].copy()
    db = db.set_index("index")
    return jgi_redundancy_check(db, catalog(cat_rows), ome_col="portal")


# --------------------------------------------------------------------------- #
# reading a version


@pytest.mark.parametrize("raw,expected", [
    ("1.0", (1, 0)), ("v2.0", (2, 0)), ("3.16", (3, 16)), ("3.2", (3, 2)),
    ("10", (1, 0)), ("20", (2, 0)), ("11", (1, 1)), ("316", (3, 16)), ("110", (1, 10)),
    ("1", (1, 0)), ("0", (0, 0)), ("", (0, 0)), (None, (0, 0)), (float("nan"), (0, 0)),
    ("0.0", (0, 0)), (1.0, (1, 0)), ("unknown", (0, 0)),
])
def test_parse_version(raw, expected):
    assert parse_version(raw) == expected


def test_a_stored_version_and_the_catalog_version_it_came_from_read_the_same():
    for catalog_version in ["0.0", "1.0", "1.1", "1.2", "2.0", "2.2", "3.16", "4.0", "6.0"]:
        stored = catalog_version.replace(".", "")
        if catalog_version == "0.0":
            stored = "0"
        assert parse_version(stored) == parse_version(catalog_version)


# --------------------------------------------------------------------------- #
# which portals are updates


def test_a_portal_that_moved_from_1_0_to_2_0_is_an_update():
    jgi_df, db, updates, old = check([("aspnig3", "Aspni1", "10", "jgi")], [("Aspni1", "2.0")])
    assert list(updates) == ["Aspni1"] and len(jgi_df) == 1
    assert "aspnig3" not in set(db["ome"])


def test_a_portal_at_the_catalog_version_is_not_an_update():
    jgi_df, db, updates, old = check([("aspnig3", "Aspni1", "20", "jgi")], [("Aspni1", "2.0")])
    assert updates == {} and len(jgi_df) == 0
    assert list(db["ome"]) == ["aspnig3"]


def test_a_portal_recorded_without_a_version_is_updated_by_any_catalog_version():
    jgi_df, db, updates, old = check([("aspnig3", "Aspni1", "0", "jgi")], [("Aspni1", "1.0")])
    assert list(updates) == ["Aspni1"]


def test_minor_versions_compare_as_numbers():
    # 3.16 is above 3.2, 1.10 above 1.9
    _, _, updates, _ = check([("a", "P1", "32", "jgi")], [("P1", "3.16")])
    assert list(updates) == ["P1"]
    _, _, updates, _ = check([("a", "P1", "316", "jgi")], [("P1", "3.2")])
    assert updates == {}
    _, _, updates, _ = check([("a", "P1", "19", "jgi")], [("P1", "1.10")])
    assert list(updates) == ["P1"]


def test_an_older_catalog_version_is_not_an_update():
    _, _, updates, _ = check([("aspnig3", "Aspni1", "30", "jgi")], [("Aspni1", "2.0")])
    assert updates == {}


def test_an_ncbi_row_sharing_a_portal_code_is_not_updated_as_a_jgi_genome():
    jgi_df, db, updates, old = check([("aspnig3", "Aspni1", "10", "ncbi")], [("Aspni1", "2.0")])
    assert updates == {} and len(jgi_df) == 0 and list(db["ome"]) == ["aspnig3"]


def test_an_update_is_recorded_as_text_with_the_old_ome_and_the_new_name():
    jgi_df, db, updates, old = check(
        [("aspnig3", "Aspni1", "10", "jgi")], [("Aspni1", "2.0"), ("Other1", "1.0"), ("Other2", "1.0")])
    old_ome, db_organism, organism = updates["Aspni1"]
    assert old_ome == "aspnig3" and db_organism == "Aspergillus_niger_X"
    assert organism == "Aspergillus_niger_Aspni1_v2.0"  # this portal's name, not the whole column
    assert "\n" not in str(organism)


def test_a_database_row_with_no_strain_does_not_stop_the_check():
    db = database([("aspnig3", "Aspni1", "10", "jgi")])
    db["strain"] = float("nan")
    db["index"] = db["assembly_acc"].copy()
    jgi_df, new_db, updates, old = jgi_redundancy_check(db.set_index("index"), catalog([("Aspni1", "2.0")]))
    assert list(updates) == ["Aspni1"]


def test_the_old_row_of_an_update_is_kept_under_the_portal_the_caller_looks_up(): 
    jgi_df, db, updates, old = check([("aspnig3", "Aspni1", "10", "jgi")], [("Aspni1", "2.0")])
    assert "Aspni1" in old and old["Aspni1"]["ome"] == "aspnig3"


# --------------------------------------------------------------------------- #
# an update that is not fetched leaves the old genome in place


def run_main(tmp_path, monkeypatch, db_rows, cat_rows, outcome):
    """`outcome`: 'ok', 'failed' or 'deferred' for every portal sent to download."""
    def stub(jgi_df, token, ome_col, output, log, log_path, dwnlds, failed, rerun,
             masked, spacer="\t\t", restore_wait=None):
        jgi_df = jgi_df.copy()
        failed, deferred = list(failed), []
        for i, row in list(jgi_df.iterrows()):
            if outcome == "ok":
                jgi_df.at[i, "fna_path"] = f"{row[ome_col]}.fna"
                jgi_df.at[i, "gff3_path"] = f"{row[ome_col]}.gff3"
            else:
                (failed if outcome == "failed" else deferred).append([row[ome_col], row["version"]])
                jgi_df = jgi_df.drop(i)
        return jgi_df, log, failed, deferred
    monkeypatch.setattr(jgi2db, "runjgi_dwnld", stub)
    monkeypatch.chdir(tmp_path)
    return jgi2db.main(catalog(cat_rows), database(db_rows), str(tmp_path), "token", date="20261123")


@pytest.mark.parametrize("outcome", ["failed", "deferred"])
def test_an_update_that_is_not_fetched_leaves_the_old_genome_in_the_database(tmp_path, monkeypatch, outcome):
    sent, new_ref, failed, deferred = run_main(
        tmp_path, monkeypatch, [("aspnig3", "Aspni1", "10", "jgi"), ("aspnig4", "Aspni2", "10", "jgi")],
        [("Aspni1", "2.0"), ("Aspni2", "1.0")], outcome)
    assert sorted(new_ref["ome"].dropna()) == ["aspnig3", "aspnig4"]
    restored = new_ref.set_index("assembly_acc").loc["Aspni1"]
    assert restored["ome"] == "aspnig3" and restored["version"] == "10"
    assert len(sent) == 0  # nothing was fetched, and no row was made up for it
    assert len(new_ref) == 2


def test_a_fetched_update_replaces_the_old_row_and_arrives_for_curation(tmp_path, monkeypatch):
    sent, new_ref, failed, deferred = run_main(
        tmp_path, monkeypatch, [("aspnig3", "Aspni1", "10", "jgi"), ("aspnig4", "Aspni2", "10", "jgi")],
        [("Aspni1", "2.0"), ("Aspni2", "1.0")], "ok")
    assert sorted(new_ref["ome"]) == ["aspnig4"]
    assert list(sent["assembly_acc"]) == ["Aspni1"] and failed == [] and deferred == []
    assert sent.iloc[0]["version"] == "2.0"
