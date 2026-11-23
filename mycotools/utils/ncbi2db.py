#! /usr/bin/env python3


# NEED to refine MTDB to consider GCA and the version in the assembly accession
#   which is important for version parsing and keeping them uniform with
#   biosample formatting in MTDB
#   NEED TO EDIT REDUNDANCY CHECK TO REFERENCE QUERIED ASSEMBLY ACCESSIONS FROM
#   BIOSAMPLES

import logging
import os
import re
import copy
import pandas as pd
from datetime import datetime
from mycotools.download.ncbi import main as ncbi_dwnld

logger = logging.getLogger(__name__)


def accession_version(accession):
    """The version of an assembly accession (GCA_000001405.29 -> 29); 0 when it
    carries none, as the non-NCBI accessions in a database do not."""
    found = re.search(r"\.(\d+)$", str(accession))
    return int(found[1]) if found else 0


def accession_base(accession):
    """An assembly accession without its version (GCA_000001405.29 ->
    GCA_000001405). One that carries no version is returned whole: cutting at
    the first "." instead (`x[: x.find(".")]`) drops the last character of it,
    so GCA_000149645 became GCA_00014964 and could never meet its own NCBI
    record. A tag after the accession (the database's GCA_000002655.v68FungiDB,
    a FungiDB annotation of the Af293 assembly) names the same assembly."""
    text = str(accession)
    found = re.match(r"(GC[AF]_\d+)(\.|$)", text)
    return found[1] if found else re.sub(r"\.\d+$", "", text)


def db_assembly_version(row):
    """The version of the assembly a database row holds. It is the suffix of
    the accession; 490 of the 6,907 genomes of the 2025-04-30 database were
    stored without one (GCA_000149645), and for those the database's `version`
    field holds the assembly version number (1, 2, 4), not a date. Where
    neither says, 1: that is the oldest an assembly can be, so a newer NCBI
    version is never missed, only possibly refetched."""
    suffix = accession_version(row["assembly_acc"])
    if suffix:
        return suffix
    stored = str(row.get("version", "")).strip()
    if re.fullmatch(r"\d{1,3}", stored) and int(stored) > 0:
        return int(stored)
    return 1


def restore_row(ref_db, old_row):
    """Put a database row that redundancy_check removed back into the reference
    database it returned (indexed by ome), as a row of that frame."""
    row = old_row.drop(labels=["index"], errors="ignore")  # the frame's index name
    restored = row.to_frame().T
    restored.index = [row["ome"]]
    restored.index.name = ref_db.index.name
    return pd.concat([ref_db, restored])


def redundancy_check(db, ncbi_df, ass_acc, duplicates={}):
    """should expect that ncbi_df (pd.DataFrame()) is only comprised of biosamples that are
    not overlapping with JGI"""
    db["index"], updates, old_rows = db.loc[:, "ome"], [], {}
    db["versionless"] = [accession_base(x) for x in list(db["assembly_acc"])]
    assembly_accs = set(db["versionless"])
    # prepare a set of assembly accessions without the version tag
    todel_db, todel_ncbi = [], []
    for i, row in ncbi_df.iterrows():  # totally based on assembly accession and version
        organism = "_".join([row["genus"], row["species"], str(row["strain"])])
        assembly_acc = row[ass_acc]
        assembly_acc_check = accession_base(assembly_acc)
        if assembly_acc_check in duplicates:  # special case duplicate
            todel_ncbi.append(i)
        elif assembly_acc_check in assembly_accs:  # assembly accession is in the
            # reference database
            check = db[db["versionless"] == assembly_acc_check]
            if len(check.index) > 1:  # multiple versions of the same accession
                # in reference
                versions = [
                    db_assembly_version(check_row) for _, check_row in check.iterrows()
                ]
                indices = list(check.index)
                max_v = max(versions)
                if versions.count(max_v) > 1:  # multiple max possible version
                    max_indices = [
                        i for i, v in enumerate(indices) if versions[i] == max_v
                    ]
                    new_check = check.iloc[max_indices]
                    dates = list(new_check["version"])
                    date_indices = list(new_check.index)
                    max_date = max(dates)
                    max_i = date_indices[dates.index(max_date)]  # will remove
                    # duplicate versions with duplicate modification dates
                else:
                    max_i = indices[versions.index(max_v)]
                more_dels = [
                    check.loc[v, "ome"] for i, v in enumerate(indices) if v != max_i
                ]
                todel_db.extend(more_dels)  # delete the overlap
                db_row = copy.deepcopy(check.loc[max_i])
            else:
                db_row = copy.deepcopy(check.iloc[0])
            db_ome = db_row["ome"]
            db_organism = "_".join(
                [db_row["genus"], db_row["species"], str(db_row["strain"])]
            )
            # An assembly is an update only when NCBI holds a newer VERSION of
            # it (GCA_x.1 -> GCA_x.2). It used to be decided by comparing the
            # report's "Modify Date" with the date the database stored, but the
            # Datasets reports the genomes are now listed from carry no
            # modification date, and the release dates they do carry are not
            # comparable with the stored ones: a record is commonly last
            # modified before it is released. Compared by date, 711 of the
            # 4,254 genomes of the 2025-04-30 database that NCBI still lists
            # came out as updated, where the previous pipeline re-versioned 1
            # of 6,790 over the preceding three months.
            if accession_version(assembly_acc) > db_assembly_version(db_row):
                updates.append(
                    [
                        db_row["assembly_acc"],
                        db_organism,
                        organism,
                        assembly_acc,
                        db_ome,
                    ]
                )
                todel_db.append(db_ome)
                old_rows[assembly_acc] = db_row
            else:
                todel_ncbi.append(i)
        else:  # new entry
            updates.append([None, None, organism, assembly_acc, None])

    db = db.set_index("index")
    todel_db.sort(reverse=True)
    todel_ncbi.sort(reverse=True)
    for i in set(todel_db):
        db = db.drop(i)
    for i in todel_ncbi:
        ncbi_df = ncbi_df.drop(i)

    return ncbi_df, db, updates, old_rows


def main(
    out_dir,
    ncbi_df=None,
    date=datetime.today().strftime("%Y%m%d"),
    ref_db=None,
    assem=True,
    prot=False,
    gff=True,
    transcript=False,
    rerun=False,
    failed_dict={},
    duplicates={},
    check_MD5=True,
    spacer="\t\t",
    chunk=25,
):

    os.chdir(out_dir)
    todel = []

    if "assembly_acc" in ncbi_df.columns:
        ass_acc = "assembly_acc"
    else:
        ass_acc = "Assembly Accession"
    if not rerun:
        ncbi_df = ncbi_df.set_index(ass_acc)
        if "-" in ncbi_df.keys():
            ncbi_df = ncbi_df.drop("-")
        prev_omes = set(ncbi_df.index)
        for failure in failed_dict:
            if failure in prev_omes:
                if failed_dict[failure]["version"]:
                    version = ncbi_df["version"][failure]
                    try:
                        prev_version = datetime.strptime(
                            failed_dict[failure]["version"], "%Y%m%d"
                        )
                        if not version > prev_version:
                            todel.append(failure)
                    except ValueError:  # if there isnt a date version in the
                        # failed log
                        todel.append(failure)
        for i in todel:
            ncbi_df = ncbi_df.drop(i)
        ncbi_df = ncbi_df.reset_index()

    logger.debug(spacer + "Redundancy check")
    update_check = {}
    if ref_db is not None:
        if len(ref_db) > 0:
            ncbi_df, ref_db, updates, old_rows = redundancy_check(
                ref_db, ncbi_df, ass_acc
            )
            output_str = (
                "old_assembly_acc\tdb_organism\tncbi_organism\tnew_assembly_acc\n"
            )
            for update in updates:
                output_str += "\t".join([str(x) for x in update]) + "\n"
            with open(out_dir + "/ncbiUpdates.tsv", "w") as out:
                out.write(output_str)
            update_check = {i[-2]: i for i in updates if i[0]}
            # dict(update_check) = {assembly_accNEW: [organism, ref organism,
            # old_assembly_acc]}
            logger.debug(spacer + "" + str(len(ncbi_df)) + " genomes to assimilate")

    if len(ncbi_df) > 0:
        logger.debug(spacer + "Initializing NCBI acquisition")
        ncbi_df, failed = ncbi_dwnld(
            assembly=assem,
            proteome=prot,
            gff3=gff,
            ncbi_df=ncbi_df,
            remove=True,
            output_path=out_dir,
            column=ass_acc,
            ncbi_column="assembly",
            check_MD5=check_MD5,
            verbose=True,
            spacer="\t\t\t",
            chunk=chunk,
        )

        logger.debug(
            spacer + "" + str(len(ncbi_df)) + " entries with assemblies and gffs"
        )
        ncbi_df = ncbi_df.rename(
            columns={
                "Release Date": "published",
                "Assembly Accession": "assembly_acc",
                "BioSample Accession": "biosample",
            }
        )
        ncbi_df["source"] = "ncbi"
        for assembly_acc_info in failed:
            if assembly_acc_info[0] in update_check:
                # the new version could not be fetched: the genome it was to
                # replace goes back into the database. (`DataFrame.append`,
                # used here before, no longer exists, so any failed update
                # stopped the whole run with an AttributeError.)
                ref_db = restore_row(ref_db, old_rows[assembly_acc_info[0]])
                del update_check[
                    assembly_acc_info[0]
                ]  # remove it from potential updates

        ref_db = ref_db.set_index("assembly_acc")  # update ome codes
        ncbi_df = ncbi_df.set_index("assembly_acc")
        for assembly_acc, update_d in update_check.items():
            if assembly_acc not in ncbi_df.index:
                # `.at` would add a row for an accession that was neither
                # fetched nor reported failed
                continue
            old_ome = update_d[-1]
            ncbi_df.at[assembly_acc, "ome"] = old_ome

        return ncbi_df.reset_index(), ref_db.reset_index(), failed
    else:
        return ncbi_df, ref_db, []  # , duplicates


if __name__ == "__main__":
    cli()
