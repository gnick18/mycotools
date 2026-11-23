# fsgdb branch

Lab integration branch of this fork. It is **not** proposed upstream: nothing here has been
opened as a pull request against xonq/mycotools. CHTC users of the lab database install from a
pinned commit of this branch.

## What it is

upstream `refactor-dev` at `db0e9b1` (draft PR #17, mycotools v2) plus these topic branches,
each also kept separate on this fork so any of them can be proposed upstream later:

| branch | what it fixes |
|---|---|
| `fix/jgi-faa-current-naming` | proteome selection for MycoCosm genomes since about 2023 (`<portal>_FilteredModels1_proteins_<date>.fasta.gz`); the old selector only knew `.aa.fasta` (64 of 100 random public genomes found, 99 with the fix) |
| `fix/jgi-token-check` | the token check used the public search endpoint, which answers 200 to any token; it now uses `/projects/`, which returns 401 for a bad one |
| `feat/jgi-seamless-token` | finds the session token without flags (environment or token file, clipboard, stored credential, prompt), renews it on a 401, and returns unfinished genomes as `deferred` instead of failed when it cannot |
| `fix/ncbi-enumerate-via-datasets-api` | the NCBI half of `mtdb update` listed genomes from `GENOME_REPORTS/eukaryotes.txt`, which NCBI stopped regenerating on 2024-10-10 (newest genome 2024-08-30), so no later genome could ever be found and nothing said so; it now pages the Datasets REST API. Because those reports carry no modification date, an assembly already in the database is an update only when NCBI holds a newer accession VERSION (`.1` to `.2`); see below |

| `fix/gff2list-legacy-rows` | upstream's Feature migration (three commits, 2026-07-25) left `gff2list` returning a `FeatureCol` built without `ingest`, so it raised "No ID obtained for feature" on any GFF3, and none of its 32 consumers reads that form. `predb2mtdb` catches only `IndexError`, so no genome could be curated, and `mtdb accession gff/locus/gbk` failed too. `gff2list` returns dict rows again; the Feature reader is `gff2features` |
| `fix/noninteractive-ncbi-key` | with no credential store the NCBI key was asked for with getpass, which raised `EOFError` in any job without a terminal (CHTC, cron); it is read from `NCBI_API_KEY` and asked for only on a terminal |

`refactor-dev` is upstream's WORK IN PROGRESS (draft PR #17): its unit suite leaves most scripts as skipped
TODO scaffolds, which is how the GFF break went unseen. Treat any command not listed under "Checked" below as unverified.

### Checked (2026-11-22, scratch copy of the 20250430 table, 6,907 genomes)

* `mtdb` link, primary path, `-l`, ome lookup and file paths, `extract` (lineage, source, inverse, by rank, ome list,
  restricted) with counts equal to an independent count of the table, flat file and SQLite (`manage --migrate`) giving
  identical extracts, `files --print`, `configure`
* `accession fa|gff|gbk|locus` byte-identical to mycotools 0.31.17 (acc2fa, acc2gff, acc2gbk, acc2locus)
* curation of a genome from raw GFF3 + FNA (`mtdb predb`): antflo1's 6,877 proteins regenerated identically to the
  database's; a genome released after the 2025-04-30 database (GCA_062053388.1) came out with 5,082 proteins, equal to NCBI's annotation count
* the metadata half of `mtdb update` (NCBI listing, report fetch, redundancy check, MycoCosm catalog and overlap),
  with no genome downloaded; `mtdb update` itself was not run

### Why the NCBI re-version rule changed

Comparing the Datasets release dates with the dates the database stored (the old `Modify Date` rule)
re-versioned 711 of the 4,254 genomes of the 2025-04-30 database that NCBI still lists, because a record
is commonly last modified before it is released. The previous pipeline re-versioned 1 of 6,790 over the
preceding three months. Same-version annotation changes are therefore no longer picked up. `datasets summary
genome taxon` is not used to list genomes: v14.15.0 returned 984 of the 6,055 annotated GenBank fungal assemblies
the API lists.

## Using the JGI session token (it lasts 12 hours and needs an ORCID login with 2FA)

* Interactive: copy the token (Data Portal avatar menu, "Copy My Session Token"), then run the
  command. No flag is needed.
* CHTC or any machine without a terminal: send the token as a file with the job and set
  `JGI_SESSION_TOKEN_FILE` to its name. The file is re-read on every renewal. A token that has
  run out makes the run stop cleanly with the unfinished genomes deferred; copy a new one and
  resubmit. Files already downloaded are kept.

## Not done yet

* request-only mode (submit tape restores and exit, so no job waits on tape holding a token)
* the NCBI listing is checked live for fungi only; plants, animals, bacteria and archaea map to NCBI taxa but were not run
* the update's genome download stages (NCBI `datasets` download, MycoCosm restore and download) were not run for the new listing; only the metadata stages were
* tape restore, assembly download and the interactive prompt have not been run live
* nothing has been run on CHTC
* `mtdb predb` and the update's `predb2mtdb` use the GFF reader on one genome of each kind only; the other readers of `gff2list` (gff/seq/phylo tools) were not run

## Rebuilding

    git checkout -b fsgdb <upstream refactor-dev commit>
    git merge --no-ff fix/jgi-faa-current-naming
    git merge --no-ff feat/jgi-seamless-token    # includes fix/jgi-token-check
    git merge --no-ff fix/ncbi-enumerate-via-datasets-api
    git merge --no-ff fix/gff2list-legacy-rows
    git merge --no-ff fix/noninteractive-ncbi-key

The first two merges conflict only in `test/unit/test_download_jgi.py`, where both sides append tests at
the end; keep both blocks. When upstream takes one of the topic branches, drop it from this list
and rebuild.
