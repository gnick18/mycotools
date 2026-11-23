#! /usr/bin/env python3
"""GFF3s that hang features off a gene instead of its RNA must curate, not fail the genome.

`cur_gff3.rename_and_organize` expects gene > RNA > exon/CDS and raises on anything else, and
`predb.main` answers a raise with "failed GFF3 curation" and drops the whole genome. In a 75-genome NCBI
test 13 genomes failed this way: 11 of the 12 recent ENA/EMBL-sourced ones (accessions GCA_9...) and 2 older
GenBank ones. The shapes below are cut down from those files (GCA_965637885.1, GCA_900074735.1,
GCA_965358675.1, GCA_905071755.1, GCA_000976875.3 and GCA_000182895.1).
"""
import re

from mycotools.utils import cur_gff3 as mod

COLS = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")


def gff(*lines):
    """Entries from tab-less shorthand: 'type start end strand attributes'."""
    out = []
    for line in lines:
        typ, start, end, strand, attrs = line.split(" ", 4)
        out.append(dict(zip(COLS, ("ctg1", "EMBL", typ, int(start), int(end), ".", strand, ".", attrs))))
    return out


def types(entries):
    return [e["type"] for e in entries]


def attr(entry, key):
    m = re.search(r"(?:^|;)" + key + r"=([^;]+)", entry["attributes"])
    return m[1] if m else None


def counted(changed):
    """The counters of normalize_gene_children that are not zero."""
    return {k: v for k, v in changed.items() if v}


def curated(entries):
    return mod.main(entries, "tst1", cur_seqids=True)


def every_feature_has_an_alias(out):
    return all(attr(e, "Alias") for e in out)


CODING = [
    "gene 100 900 + ID=gene-A;Name=A;gbkey=Gene;gene_biotype=protein_coding;locus_tag=A",
    "mRNA 100 900 + ID=rna-A;Parent=gene-A;gbkey=mRNA;locus_tag=A",
    "exon 100 900 + ID=exon-A-1;Parent=rna-A;gbkey=mRNA;locus_tag=A",
    "CDS 100 900 + ID=cds-P1;Parent=rna-A;gbkey=CDS;locus_tag=A;protein_id=P1",
]


def test_a_gff3_in_the_expected_shape_is_left_alone():
    entries = gff(*CODING)
    out, changed = mod.normalize_gene_children(entries)
    assert out == entries
    assert counted(changed) == {}


def test_the_input_list_is_not_modified():
    entries = gff(*CODING, "CDS 100 900 + ID=cds-P2;Parent=gene-A;gbkey=CDS;locus_tag=A;protein_id=P2")
    before = [dict(e) for e in entries]
    mod.normalize_gene_children(entries)
    assert entries == before


def test_a_trna_exon_repeated_under_its_gene_is_removed_and_the_genome_curates():
    # GCA_965637885.1: the extra exon is what failed the genome
    entries = gff(
        *CODING,
        "gene 1000 1070 + ID=gene-T;Name=T;gbkey=Gene;gene_biotype=tRNA;locus_tag=T",
        "tRNA 1000 1070 + ID=rna-T;Parent=gene-T;gbkey=tRNA;locus_tag=T;product=tRNA-Phe",
        "exon 1000 1070 + ID=exon-T-1;Parent=rna-T;gbkey=tRNA;locus_tag=T",
        "exon 1000 1070 + ID=id-T;Parent=gene-T;gbkey=exon;locus_tag=T",
    )
    out, changed = mod.normalize_gene_children(entries)
    assert changed["exons_removed"] == 1
    assert [e for e in out if attr(e, "ID") == "id-T"] == []
    result = curated(entries)
    assert every_feature_has_an_alias(result)
    assert types(result).count("exon") == 2  # the mRNA's and the tRNA's, once each


def test_the_repeated_exon_of_a_coding_gene_is_removed_too_so_exons_are_not_doubled():
    entries = gff(*CODING, "exon 100 900 + ID=id-A;Parent=gene-A;gbkey=exon;locus_tag=A")
    result = curated(entries)
    assert types(result).count("exon") == 1


def test_an_exon_under_a_gene_that_differs_from_the_rna_exons_replaces_the_rnas_whole_gene_exon():
    # the RNA of CODING has one exon over the whole gene, the kind a converter writes as a placeholder
    entries = gff(*CODING, "exon 100 500 + ID=id-A;Parent=gene-A;gbkey=exon;locus_tag=A")
    out, changed = mod.normalize_gene_children(entries)
    assert counted(changed) == {"exons_replaced": 1}
    assert [(e["start"], e["end"]) for e in out if e["type"] == "exon"] == [(100, 500)]


PSEUDOGENE = [
    "pseudogene 100 900 + ID=gene-X;Name=X;gbkey=Gene;gene_biotype=pseudogene;locus_tag=X;pseudo=true",
    "CDS 100 400 + ID=cds-X;Parent=gene-X;gbkey=CDS;locus_tag=X;pseudo=true",
    "CDS 500 900 + ID=cds-X;Parent=gene-X;gbkey=CDS;locus_tag=X;pseudo=true",
    "mRNA 100 900 + ID=rna-X;Parent=gene-X;gbkey=mRNA;locus_tag=X;pseudo=true",
    "exon 100 900 + ID=exon-X-1;Parent=rna-X;gbkey=mRNA;locus_tag=X;pseudo=true",
]


def test_the_cds_of_a_pseudogene_follows_the_pseudogenes_only_rna():
    # GCA_900074735.1, GCA_964275055.1, GCA_965358675.1: the CDS comes before the mRNA and names the gene as parent
    out, changed = mod.normalize_gene_children(gff(*PSEUDOGENE))
    assert counted(changed) == {"cds_reparented": 2}
    order = [(e["type"], attr(e, "Parent")) for e in out]
    assert order.index(("mRNA", "gene-X")) < order.index(("CDS", "rna-X"))
    assert [e["start"] for e in out if e["type"] == "CDS"] == [100, 500]


def test_a_pseudogene_cds_is_kept_and_named_like_the_labs_other_pseudogenes():
    result = curated(gff(*PSEUDOGENE))
    cds = [e for e in result if e["type"] == "CDS"]
    assert len(cds) == 2
    assert every_feature_has_an_alias(result)
    assert {attr(c, "Alias") for c in cds} == {"tst1_pseudogene1"}


def test_a_second_protein_under_the_gene_gets_its_own_mrna():
    # GCA_000182895.1: a CDS for EFI28603.1 hangs off a gene whose mRNA already carries EFI28602.1
    entries = gff(*CODING, "CDS 120 300 + ID=cds-P2;Parent=gene-A;gbkey=CDS;locus_tag=A;protein_id=P2",
                  "CDS 400 600 + ID=cds-P2;Parent=gene-A;gbkey=CDS;locus_tag=A;protein_id=P2")
    out, changed = mod.normalize_gene_children(entries)
    assert counted(changed) == {"cds_reparented": 2, "mrnas_added": 1}
    new = [e for e in out if e["type"] == "mRNA" and attr(e, "ID") != "rna-A"][0]
    assert (new["start"], new["end"], attr(new, "Parent")) == (120, 600, "gene-A")
    result = curated(entries)
    assert every_feature_has_an_alias(result)
    assert {attr(c, "Alias") for c in result if c["type"] == "CDS"} == {"tst1_P1", "tst1_P2"}


def test_two_proteins_on_one_mrna_are_split_onto_two_mrnas():
    # GCA_965358675.1: ORF1 and ORF2 of a retrotransposon on one mRNA
    entries = gff(*CODING, "CDS 120 300 + ID=cds-P2;Parent=rna-A;gbkey=CDS;locus_tag=A;protein_id=P2")
    out, changed = mod.normalize_gene_children(entries)
    assert counted(changed) == {"cds_reparented": 1, "mrnas_added": 1}
    assert types(out).count("mRNA") == 2
    result = curated(entries)
    assert every_feature_has_an_alias(result)
    assert {attr(c, "Alias") for c in result if c["type"] == "CDS"} == {"tst1_P1", "tst1_P2"}


def test_the_segments_of_one_protein_stay_on_one_mrna():
    entries = gff(
        "gene 100 900 + ID=gene-A;Name=A;gbkey=Gene;gene_biotype=protein_coding;locus_tag=A",
        "mRNA 100 900 + ID=rna-A;Parent=gene-A;gbkey=mRNA;locus_tag=A",
        "exon 100 300 + ID=exon-A-1;Parent=rna-A;gbkey=mRNA;locus_tag=A",
        "exon 500 900 + ID=exon-A-2;Parent=rna-A;gbkey=mRNA;locus_tag=A",
        "CDS 100 300 + ID=cds-P1;Parent=rna-A;gbkey=CDS;locus_tag=A;protein_id=P1",
        "CDS 500 900 + ID=cds-P1;Parent=rna-A;gbkey=CDS;locus_tag=A;protein_id=P1",
    )
    out, changed = mod.normalize_gene_children(entries)
    assert out == entries and not any(changed.values())


def test_an_mrna_with_no_cds_beside_one_that_has_a_protein_gets_an_alias():
    # GCA_000976875.3 and GCA_964275055.1: a second, protein-less mRNA of the same gene
    entries = gff(
        *CODING,
        "mRNA 101 900 + ID=rna-A-2;Parent=gene-A;gbkey=mRNA;locus_tag=A",
        "exon 101 900 + ID=exon-A-2-1;Parent=rna-A-2;gbkey=mRNA;locus_tag=A",
    )
    result = curated(entries)
    assert every_feature_has_an_alias(result)
    mrnas = [e for e in result if e["type"] == "mRNA"]
    assert len({attr(m, "Alias") for m in mrnas}) == 2


def test_a_utr_of_a_gene_with_several_mrnas_is_placed_under_the_gene():
    # GCA_905071755.1: the UTR carries the aliases of all of the gene's mRNAs
    entries = gff(
        "gene 100 900 + ID=gene-A;Name=A;gbkey=Gene;gene_biotype=protein_coding;locus_tag=A",
        "mRNA 100 900 + ID=rna-A-1;Parent=gene-A;gbkey=mRNA;locus_tag=A",
        "exon 100 900 + ID=exon-A-1;Parent=rna-A-1;gbkey=mRNA;locus_tag=A",
        "CDS 150 850 + ID=cds-P1;Parent=rna-A-1;gbkey=CDS;locus_tag=A;protein_id=P1",
        "mRNA 100 900 + ID=rna-A-2;Parent=gene-A;gbkey=mRNA;locus_tag=A",
        "exon 100 900 + ID=exon-A-2;Parent=rna-A-2;gbkey=mRNA;locus_tag=A",
        "CDS 160 850 + ID=cds-P2;Parent=rna-A-2;gbkey=CDS;locus_tag=A;protein_id=P2",
        "three_prime_UTR 851 900 + ID=id-utr;Parent=gene-A;gbkey=3'UTR;locus_tag=A",
    )
    result = curated(entries)
    utr = [e for e in result if e["type"] == "three_prime_UTR"][0]
    gene = [e for e in result if e["type"] == "gene"][0]
    assert attr(utr, "Parent") == attr(gene, "ID")


# ---- genes with no RNA, and exons under a gene that are not repeats (11 more genomes of the 825 first fetched) ----

def norna_cds(ident="gene-M", prots=("P1",), exons=False, introns=False, extra=()):
    """A gene with no RNA whose CDSs are its own children (mitochondrial genes, S. pombe, AUGUSTUS output)."""
    lines = [f"gene 100 900 + ID={ident};Name=M;gbkey=Gene;gene_biotype=protein_coding;locus_tag=M"]
    for k, p in enumerate(prots):
        a, b = 100 + 400 * k, 300 + 400 * k
        lines.append(f"CDS {a} {b} + ID=cds-{p};Parent={ident};gbkey=CDS;locus_tag=M;protein_id={p}")
    if exons:
        lines += [f"exon 100 300 + ID=id-M;Parent={ident};gbkey=exon;locus_tag=M",
                  f"exon 500 900 + ID=id-M-2;Parent={ident};gbkey=exon;locus_tag=M"]
    if introns:
        lines.append(f"intron 301 499 + ID=id-M-i;Parent={ident};gbkey=intron;locus_tag=M")
    return list(lines) + list(extra)


def test_a_gene_with_cds_and_introns_and_no_rna_gets_an_mrna_with_exons():
    # GCA_000238015.2 (a mitochondrial coxIII), GCA_000765035.1 (10,000 AUGUSTUS genes)
    entries = gff(*norna_cds(introns=True))
    out, changed = mod.normalize_gene_children(entries)
    assert changed["mrnas_added"] == 1 and changed["exons_added"] == 1 and types(out).count("intron") == 0
    result = curated(entries)
    assert every_feature_has_an_alias(result)
    assert {attr(c, "Alias") for c in result if c["type"] == "CDS"} == {"tst1_P1"}


def test_the_genes_own_exons_become_the_exons_of_its_new_mrna():
    entries = gff(*norna_cds(exons=True))                       # GCA_000313545.1
    out, changed = mod.normalize_gene_children(entries)
    assert changed["mrnas_added"] == 1 and changed["exons_added"] == 0
    assert [e["start"] for e in out if e["type"] == "exon"] == [100, 500]
    assert every_feature_has_an_alias(curated(entries))


def test_two_proteins_on_a_gene_with_no_rna_get_an_mrna_each():
    entries = gff(*norna_cds(prots=("P1", "P2")))               # GCA_000002945.3, S. pombe mei3
    out, changed = mod.normalize_gene_children(entries)
    assert changed["mrnas_added"] == 2
    result = curated(entries)
    assert every_feature_has_an_alias(result)
    assert {attr(c, "Alias") for c in result if c["type"] == "CDS"} == {"tst1_P1", "tst1_P2"}


def test_a_gene_with_one_protein_and_nothing_else_is_left_alone():
    entries = gff(*norna_cds())                                 # thousands per genome already curate as they are
    out, changed = mod.normalize_gene_children(entries)
    assert out == entries and not any(changed.values())


def test_a_pseudogene_with_only_exons_gets_one_mrna_over_them():
    entries = gff("pseudogene 100 900 - ID=gene-X;Name=X;gbkey=Gene;gene_biotype=pseudogene;locus_tag=X;pseudo=true",
                  "exon 100 300 - ID=id-X;Parent=gene-X;gbkey=exon;locus_tag=X;pseudo=true",
                  "exon 400 900 - ID=id-X-2;Parent=gene-X;gbkey=exon;locus_tag=X;pseudo=true")   # GCA_036370985.1
    out, changed = mod.normalize_gene_children(entries)
    assert changed["mrnas_added"] == 1
    assert every_feature_has_an_alias(curated(entries))


def test_a_gene_with_only_introns_loses_them():
    entries = gff("gene 100 900 + ID=gene-C;Name=C;gbkey=Gene;gene_biotype=other;locus_tag=C",
                  "intron 301 499 + ID=id-C;Parent=gene-C;gbkey=intron;locus_tag=C")            # GCA_000313545.1 cox1
    out, changed = mod.normalize_gene_children(entries)
    assert changed["introns_removed"] == 1 and types(out) == ["gene"]
    assert every_feature_has_an_alias(curated(entries))


def test_a_whole_gene_exon_beside_an_rna_with_its_real_exons_is_removed():
    entries = gff(
        "gene 100 900 + ID=gene-A;Name=A;gbkey=Gene;gene_biotype=protein_coding;locus_tag=A",
        "mRNA 100 900 + ID=rna-A;Parent=gene-A;gbkey=mRNA;locus_tag=A",
        "exon 100 300 + ID=exon-A-1;Parent=rna-A;gbkey=mRNA;locus_tag=A",
        "exon 500 900 + ID=exon-A-2;Parent=rna-A;gbkey=mRNA;locus_tag=A",
        "CDS 100 300 + ID=cds-P1;Parent=rna-A;gbkey=CDS;locus_tag=A;protein_id=P1",
        "CDS 500 900 + ID=cds-P1;Parent=rna-A;gbkey=CDS;locus_tag=A;protein_id=P1",
        "exon 100 900 + ID=id-A;Parent=gene-A;gbkey=exon;locus_tag=A")                         # GCA_052576325.1
    out, changed = mod.normalize_gene_children(entries)
    assert changed["exons_removed"] == 1
    assert types(curated(entries)).count("exon") == 2


def test_the_genes_own_exons_replace_a_single_whole_gene_exon_of_its_rna():
    entries = gff(
        "gene 100 900 + ID=gene-A;Name=A;gbkey=Gene;gene_biotype=protein_coding;locus_tag=A",
        "mRNA 100 900 + ID=rna-A;Parent=gene-A;gbkey=mRNA;locus_tag=A",
        "exon 100 900 + ID=exon-A-1;Parent=rna-A;gbkey=mRNA;locus_tag=A",
        "CDS 100 300 + ID=cds-P1;Parent=gene-A;gbkey=CDS;locus_tag=A;protein_id=P1",
        "CDS 500 900 + ID=cds-P1;Parent=gene-A;gbkey=CDS;locus_tag=A;protein_id=P1",
        "exon 100 300 + ID=id-A;Parent=gene-A;gbkey=exon;locus_tag=A",
        "exon 500 900 + ID=id-A-2;Parent=gene-A;gbkey=exon;locus_tag=A")                       # GCA_001373355.1
    out, changed = mod.normalize_gene_children(entries)
    assert changed["exons_replaced"] == 1 and changed["exons_removed"] == 0
    assert sorted(e["start"] for e in out if e["type"] == "exon") == [100, 500]
    assert every_feature_has_an_alias(curated(entries))


def test_a_gene_written_as_several_lines_with_one_id_keeps_its_new_mrna():
    # NCBI splits a multi-part gene into several gene lines; the curator starts the gene over at each one, so the
    # new mRNA must follow the last (GCA_000313545.1 lost 6 proteins before this)
    entries = gff("gene 181 249 - ID=gene-Q;Name=Q;gbkey=Gene;gene_biotype=protein_coding;locus_tag=Q",
                  "gene 117 131 - ID=gene-Q;Name=Q;gbkey=Gene;gene_biotype=protein_coding;locus_tag=Q",
                  "CDS 181 249 - ID=cds-P1;Parent=gene-Q;gbkey=CDS;locus_tag=Q;protein_id=P1",
                  "CDS 117 131 - ID=cds-P1;Parent=gene-Q;gbkey=CDS;locus_tag=Q;protein_id=P1",
                  "exon 117 131 - ID=id-Q;Parent=gene-Q;gbkey=exon;locus_tag=Q",
                  "exon 181 249 - ID=id-Q-2;Parent=gene-Q;gbkey=exon;locus_tag=Q")
    result = curated(entries)
    assert {attr(c, "Alias") for c in result if c["type"] == "CDS"} == {"tst1_P1"}
    assert types(result).count("CDS") == 2
