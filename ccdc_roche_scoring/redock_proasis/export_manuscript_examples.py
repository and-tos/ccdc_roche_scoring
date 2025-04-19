import pandas as pd
from pathlib import Path
from ccdc import io


def export_artifacts(pdb_rdr, vina_rdr):
    bs_ids = ["7UHI_010", "6G8F_010", "1TVE_010", "1YAT_001", "7LVI_001", "4R68_001"]
    Path("examples/artifacts").mkdir(exist_ok=True)
    for bs_id in bs_ids:
        with io.MoleculeWriter(
            f"examples/artifacts/artifacts_experimental_{bs_id}.mol2"
        ) as w:
            ccdc_mol = pdb_rdr.molecule(bs_id)
            ccdc_mol.identifier = ccdc_mol.identifier + "_experimental"
            w.write(ccdc_mol)

        with io.MoleculeWriter(
            f"examples/artifacts/artifacts_docking_{bs_id}.mol2"
        ) as w:
            try:
                ccdc_mol = vina_rdr.molecule(bs_id)
                ccdc_mol.identifier = ccdc_mol.identifier + "_docking"
                w.write(ccdc_mol)
            except:
                continue

    return


def export_rf_examples(pdb_rdr, vina_rdr):
    bs_ids = ["3PZE_001", "4RCH_001", "6ROT_001", "4X8T_001", "4U43_001", "1ZDP_001"]
    example_dir = Path("examples/rf_examples")
    Path("examples").mkdir(exist_ok=True)
    example_dir.mkdir(exist_ok=True)

    for bs_id in bs_ids:
        with io.MoleculeWriter(example_dir / f"{bs_id}_experimental.mol2") as w:
            ccdc_mol = pdb_rdr.molecule(bs_id)
            ccdc_mol.identifier = ccdc_mol.identifier + "_experimental"
            w.write(ccdc_mol)

    for bs_id in bs_ids:
        with io.MoleculeWriter(example_dir / f"{bs_id}_docking.mol2") as w:
            try:
                ccdc_mol = vina_rdr.molecule(bs_id)
                ccdc_mol.identifier = ccdc_mol.identifier + "_docking"
                w.write(ccdc_mol)
            except:
                continue

    return


def export_low_q_pdbs(df, pdb_rdr, substructure_name):
    df = df[(df["q_pdb"] < 0.01) & df["identifier"].str.endswith("_001")]

    df = df.sort_values("RMSD_vina", ascending=False)
    example_dir = Path("examples/artifacts/")
    example_id = df["identifier"].values[0]

    pdb_out = str(example_dir / f"{example_id}_{substructure_name}_pdb.mol2")

    with io.MoleculeWriter(pdb_out) as w:
        w.write(pdb_rdr.molecule(example_id))
    return


def export_matched_pairs(df, vina_rdr, pdb_rdr, substructure_name):
    df = df[
        (df["q_vina"] < 0.01)
        & (df["q_pdb"] > 0.2)
        & (df["identifier"].str.endswith("_001"))
        & (df["identifier"].str.len() == 8)
        & (df["RMSD_vina"] > 1)
    ]
    if substructure_name == "phenylalanine":
        df = df[(df["TOR1_vina"] > 90) & (df["TOR1_vina"] < 140)]
        df = df[df["TOR2_vina"] < 45]

    df = df.sort_values("RMSD_vina", ascending=False).drop_duplicates("identifier")

    example_dir = Path("examples/torsion_examples/")
    example_dir.mkdir(exist_ok=True)

    for example_id in df["identifier"].values[0:50]:
        vina_out = str(example_dir / f"{substructure_name}_docking_{example_id}.mol2")
        pdb_out = str(
            example_dir / f"{substructure_name}_experimental_{example_id}.mol2"
        )
        with io.MoleculeWriter(vina_out) as w:
            ccdc_mol = vina_rdr.molecule(example_id)
            ccdc_mol.identifier = ccdc_mol.identifier + "_docking"
            w.write(ccdc_mol)
        with io.MoleculeWriter(pdb_out) as w:
            ccdc_mol = pdb_rdr.molecule(example_id)
            ccdc_mol.identifier = ccdc_mol.identifier + "_experimental"
            w.write(ccdc_mol)

    return


def main():
    vina_rdr = io.MoleculeReader("vina_top_docked_complexes_pub.csdsql")
    pdb_rdr = io.MoleculeReader(
        "../protonate3d_2024-12-03/full_p2cq_pub_2024-12-03.csdsql"
    )

    # # Export RF example binding sites
    # export_rf_examples(pdb_rdr, vina_rdr)

    # # Export artifact binding sites
    export_artifacts(pdb_rdr, vina_rdr)

    druglike_df = pd.read_parquet("native_and_docking_pose_check_qmean.gzip")
    druglike_df = druglike_df[(druglike_df["MW"] > 300) & (druglike_df["MW"] < 650)]
    druglike_df = druglike_df[(druglike_df["NumRotatableBonds"] < 10)]
    druglike_df = druglike_df[druglike_df["HBD"] < 4]
    druglike_df = druglike_df[druglike_df["HBA"] < 10]
    druglike_df = druglike_df[druglike_df["NumAromaticRings"] < 5]
    druglike_df = druglike_df[druglike_df["num_pos_charges"] < 2]
    druglike_df = druglike_df[druglike_df["num_neg_charges"] < 2]

    for substructure_name in [
        "butane",
        "phenyl-methyl-pyrrazole",
        "phenyl-fluoro-phenyl",
        "ethylsulfanylethane",
        "phenyl_sulfonyl",
        "5-methyl-1-phenyl-pyrazole",
        "2-methylbutane",
        "phenol",
        "ethoxybenzene",
        "N-benzylacetamide",
        "benzyloxybenzene",
        "phenylalanine",
        "N-phenylformamide",
    ]:
        df = pd.read_csv(Path("q_values") / f"{substructure_name}.csv")
        docking_df = df[df["Database"] == "docking"]
        docking_df = docking_df[docking_df["identifier"].isin(druglike_df["bs_id"])]
        pdb_df = df[df["Database"] == "reference"]
        df = docking_df.join(
            pdb_df.set_index("identifier"),
            on="identifier",
            rsuffix="_pdb",
            lsuffix="_vina",
        )

        # export_low_q_pdbs(df, pdb_rdr, substructure_name)
        export_matched_pairs(df, vina_rdr, pdb_rdr, substructure_name)

    return


if __name__ == "__main__":
    main()
