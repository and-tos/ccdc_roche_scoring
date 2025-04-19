#!/usr/bin/env python

########################################################################################################################

import shutil
from pathlib import Path

from ccdc import io
from rdkit.Chem import PandasTools
from datetime import datetime

########################################################################################################################


def remove_key(attribute, e):
    for key_string in [
        "Gold.PLP.C",
        "Gold.PLP.p",
        "Gold.PLP.PLP",
        "Gold.Protein",
        "Gold.Chemscore",
        "Gold.Id.",
        "Gold.PLP.S",
        "ABUNDANCE",
    ]:
        if key_string in attribute:
            del e.attributes[attribute]
    return e


def write_concat_sdf(docked_solns):
    # with io.EntryWriter("test_best_docking_solutions.sdf") as w:
    for docked_soln in docked_solns:
        if docked_soln.stat().st_size / (1024 * 1024) > 1:
            docked_soln.unlink()
            continue
        try:
            e = io.EntryReader(str(docked_soln))[0]
        except:
            docked_soln.unlink()
            print("bad file")
            continue
        if "RMSD_to_mcs" not in e.attributes:
            docked_soln.unlink()
            continue
        # if float(e.attributes["RMSD_to_mcs"]) >= 2:
        # continue

        # attributes = list(e.attributes.keys()).copy()
        # for attribute in attributes:
        #     e = remove_key(attribute, e)
        # w.write(e)
    return


def main():
    print("Collecting results...")

    concat_sdf = "concat.sdf"
    with open(concat_sdf, "wb") as outfile:
        rescored_exists = False
        docked_solns = Path(".").glob("docking_job_*/rescored_ligands.sdf")
        # write_concat_sdf(docked_solns)
        for filename in docked_solns:
            rescored_exists = True
            with open(filename, "rb") as infile:
                shutil.copyfileobj(infile, outfile)

        if not rescored_exists:
            docked_solns = Path(".").glob("docking_job_*/*/best_soln_*.sdf")
            # write_concat_sdf(docked_solns)
            for filename in docked_solns:
                with open(filename, "rb") as infile:
                    shutil.copyfileobj(infile, outfile)

    df = PandasTools.LoadSDF(concat_sdf, removeHs=False)
    df = df.astype({"RMSD_to_mcs": float, "Gold.PLP.Fitness": float})
    df = df[df["RMSD_to_mcs"] <= 2]
    df = df.sort_values("Gold.PLP.Fitness", ascending=False).drop_duplicates("ID")
    df = df.drop(columns=[c for c in df.columns if "Gold.PLP.C" in c])
    df = df.drop(columns=[c for c in df.columns if "Gold.PLP.p" in c])
    df = df.drop(columns=[c for c in df.columns if "Gold.PLP.PLP" in c])
    df = df.drop(columns=[c for c in df.columns if "Gold.Protein" in c])
    df = df.drop(columns=[c for c in df.columns if "Gold.Chemscore" in c])
    df = df.drop(columns=[c for c in df.columns if "Gold.Id." in c])
    df = df.drop(columns=[c for c in df.columns if "Gold.PLP.S" in c])
    df = df.drop(columns=[c for c in df.columns if "ABUNDANCE" in c])
    df["Date"] = datetime.today().strftime("%Y-%m-%d")
    PandasTools.WriteSDF(df, "best_docking_solutions.sdf", properties=list(df.columns))

    docking_folders = df["docking_folder"].to_list()
    docked_pockets = [
        str(Path(docking_folder).parent / "best_soln_pocket.mol2")
        for docking_folder in docking_folders
    ]

    rdr = io.EntryReader(docked_pockets)
    with io.EntryWriter("pockets.mol2") as w:
        for docked_pocket in rdr:
            w.write(docked_pocket)


if __name__ == "__main__":
    main()
