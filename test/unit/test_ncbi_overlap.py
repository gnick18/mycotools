#! /usr/bin/env python3
"""Which NCBI genomes the NCBI half of an update goes on to download.

`rm_ncbi_overlap` sorts the listing into NCBI genomes that are also MycoCosm
genomes (left to the MycoCosm half) and the rest. `exec_rm_overlap` returned
(overlap, remaining) and both callers unpacked (ncbi_df, overlap), so the NCBI
half downloaded the MycoCosm copies and ignored the NCBI-only genomes. The
2024-12-23 update added 173 genomes through its NCBI half and every one was on
its own list of MycoCosm copies (redundant_ncbi.tsv); only ~10 MycoCosm
downloads had failed, so the fallback that re-adds a failed copy cannot explain
it. Version 0.31.17 has the same code.
"""
import pandas as pd

from mycotools.mtdb import update as mod


def listing(rows):
    df = pd.DataFrame(
        [{"assembly_acc": acc, "genus": genus, "species": species, "strain": strain,
          "BioSample Accession": f"SAMN{n}"} for n, (acc, genus, species, strain) in enumerate(rows)]
    )
    return df.set_index("assembly_acc", drop=False)


def catalog(rows):
    return pd.DataFrame([{"portal": portal, "genus": genus, "species": species, "strain": strain}
                         for portal, genus, species, strain in rows])


def meta(accessions, name="Asm", submitter="University"):
    return {acc: {"accession": f"{name}{i}", "submitter": submitter} for i, acc in enumerate(accessions)}


def test_exec_rm_overlap_returns_what_remains_first_then_the_overlap():
    df = pd.DataFrame({"assembly_acc": ["A", "B", "C", "D"]}).set_index("assembly_acc", drop=False)
    remaining, overlap = mod.exec_rm_overlap(df, ["B", "D"])
    assert list(remaining.index) == ["A", "C"]
    assert list(overlap.index) == ["B", "D"]


def test_the_ncbi_genomes_that_remain_are_the_ones_that_are_not_mycocosm_copies():
    ncbi = listing([("GCA_1.1", "Aspergillus", "niger", "ATCC1015"),    # a MycoCosm genome (same genus/species/strain)
                    ("GCA_2.1", "Fusarium", "oxysporum", "Fo47"),     # submitted by JGI
                    ("GCA_3.1", "Candida", "albicans", "SC5314")])    # NCBI only
    jgi = catalog([("Aspni5", "Aspergillus", "niger", "ATCC1015")])
    acc2meta = {"GCA_1.1": {"accession": "ASM1", "submitter": "Univ"},
                "GCA_2.1": {"accession": "ASM2", "submitter": "DOE Joint Genome Institute"},
                "GCA_3.1": {"accession": "ASM3", "submitter": "Univ"}}
    remaining, jgi2ncbi, jgi2bio, fails, overlap, todel = mod.rm_ncbi_overlap(ncbi, jgi, {}, set(), acc2meta)
    assert list(remaining["assembly_acc"]) == ["GCA_3.1"]
    assert sorted(overlap["assembly_acc"]) == ["GCA_1.1", "GCA_2.1"]


def test_a_genome_the_table_holds_as_ncbi_is_not_taken_for_a_mycocosm_copy_in_its_next_version():
    """The ledger of NCBI-only accessions has GCA_1.1; the listing now has GCA_1.2 of the same
    genome, which also matches a MycoCosm portal by name. It is the table's own genome."""
    db = pd.DataFrame([{"assembly_acc": "GCA_1.1", "source": "ncbi"}, {"assembly_acc": "Aspni5", "source": "jgi"}])
    ncbi = listing([("GCA_1.2", "Aspergillus", "niger", "ATCC1015"), ("GCA_9.1", "Aspergillus", "niger", "Other")])
    known = mod.known_ncbi_accessions(db, ncbi, {"GCA_5.1"})
    assert known == {"GCA_5.1", "GCA_1.2"}
    jgi = catalog([("Aspni5", "Aspergillus", "niger", "ATCC1015"), ("Aspni6", "Aspergillus", "niger", "Other")])
    remaining, _, _, _, overlap, _ = mod.rm_ncbi_overlap(ncbi, jgi, {}, known, meta(ncbi["assembly_acc"]))
    assert list(remaining["assembly_acc"]) == ["GCA_1.2"]   # kept, so redundancy_check can see it is newer
    assert list(overlap["assembly_acc"]) == ["GCA_9.1"]     # the other is still a MycoCosm copy


def test_a_jgi_genome_held_in_the_table_does_not_exempt_an_ncbi_record():
    db = pd.DataFrame([{"assembly_acc": "Aspni5", "source": "jgi"}])
    ncbi = listing([("GCA_1.1", "Aspergillus", "niger", "ATCC1015")])
    assert mod.known_ncbi_accessions(db, ncbi, set()) == set()
