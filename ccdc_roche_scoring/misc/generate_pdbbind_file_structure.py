########################################################################################################################

import argparse
from pathlib import Path

import generate_metrics_table
import pandas as pd
from ccdc import docking, io
from rf_statistics.los_descriptors import _cut_out_binding_site_by_distance

########################################################################################################################


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="Submit Calculate Rf values on Basel HPC cluster.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument(
        "-np",
        help="Number of parallel processes for multiprocessing.",
        type=int,
        default=24,
    )

    return parser.parse_args()


def write_pdb_files(iteritem):
    cnt, ligand_file = iteritem
    docking_folder = Path(ligand_file).parent
    gold_conf = list(docking_folder.glob("*.conf"))
    if gold_conf:
        gold_conf = str(gold_conf[0])
        docking_settings = docking.Docker.Settings().from_file(gold_conf)
        docking_results = docking.Docker.Results(docking_settings)
        docked_ligand_reader = docking_results.DockedLigandReader(
            ligand_file, docking_settings
        )
        csd_ligand_entry = docked_ligand_reader[0]

        # setup protein-ligand-complex
        protein = docking_results.make_complex(csd_ligand_entry)
        apo_protein = protein.copy()
        apo_protein.remove_ligand(apo_protein.ligands[0].identifier)
        apo_protein.remove_unknown_atoms()

        apo_pocket = protein.copy()
        apo_pocket = _cut_out_binding_site_by_distance(
            apo_pocket, csd_ligand_entry.molecule
        )
        apo_pocket.remove_ligand(apo_pocket.ligands[0].identifier)
        apo_pocket.remove_unknown_atoms()

        strucid = docking_folder.stem.split("_")[0]
        folder = Path(f"pde-10_pdb_bind_format/{strucid}_{cnt}")
        if not folder.is_dir():
            folder.mkdir()
        with io.EntryWriter(folder / f"{strucid}_{cnt}_protein.pdb") as w:
            w.write(apo_protein)
        with io.EntryWriter(folder / f"{strucid}_{cnt}_pocket.pdb") as w:
            w.write(apo_pocket)
        with io.EntryWriter(folder / f"{strucid}_{cnt}_ligand.sdf") as w:
            w.write(csd_ligand_entry)
        with io.EntryWriter(folder / f"{strucid}_{cnt}_ligand.mol2") as w:
            w.write(csd_ligand_entry)


def main():
    args = parse_args()
    best_docked_df = pd.read_csv(
        "docking_template_date/docked_rf_count_best_docked.csv"
    )
    folder = Path(f"pde-10_pdb_bind_format")
    if not folder.is_dir():
        folder.mkdir()

    # nproc = args.np
    # if nproc == 1:
    #     for iteritem in best_docked_df['ligand_file'].items():
    #         print(iteritem)
    #         write_pdb_files(iteritem)
    # else:
    #     parallel_return_df = partial(write_pdb_files)
    #     pool = ProcessingPool(nproc)
    #     pool.map(parallel_return_df, best_docked_df['ligand_file'].items())

    PDB_codes = []
    dicts = []
    for d in Path(f"pde-10_pdb_bind_format_labeled").glob("*_*"):
        if d.is_dir():
            ligand_file = list(d.glob("*.sdf"))
            if ligand_file:
                ligand_file = ligand_file[0]
                ligand_entry = io.EntryReader(str(ligand_file))[0]
                dicts.append(ligand_entry.attributes)
                PDB_codes.append(d.stem)
    attributes_dict = {
        k: [d[k] for d in dicts]
        for k in dicts[0]
        if k
        not in [
            "Gold.PLP.SBar",
            "Gold.Protein.RotatedWaterAtoms",
            "NETCHARGE",
            "Proasis ID",
            "series",
        ]
        and "Gold.Protein" not in k
    }
    attributes_dict["PDB_code"] = PDB_codes
    attributes_dict["release_year"] = attributes_dict[
        "PDE10_FULL_LEN_CGMP_SPA_IC50_h-PDE10A(14-779)-E.Coli-c: AP005128;Min;EXP Test Date"
    ]
    del attributes_dict[
        "PDE10_FULL_LEN_CGMP_SPA_IC50_h-PDE10A(14-779)-E.Coli-c: AP005128;Min;EXP Test Date"
    ]
    df = pd.DataFrame(attributes_dict)
    series_df = pd.read_csv("series_assignment.csv")
    series_df = series_df[["SRN", "series"]]
    df = df.join(series_df.set_index("SRN"), on="SRN")
    # df.to_csv(Path(f'pde-10_pdb_bind_format') / 'pde-10_data.csv', index=False)

    pdb_proasis_df = pd.read_csv("proasis_pdb.csv").rename(
        columns={"IDs": "PDB_code", "Label": "strucid"}
    )
    pdb_proasis_df["PDB_code"] = pdb_proasis_df["PDB_code"].str.lower()
    proasis_date_df = pd.read_csv("proasis_date.csv")[["STRUCID", "DATESOLVED"]]
    pdb_proasis_df = pdb_proasis_df.join(
        proasis_date_df.set_index("STRUCID"), on="strucid"
    )
    pdb_proasis_df = pdb_proasis_df.drop(columns="strucid")
    si_df = df[
        [
            "SRN",
            "pde10_h_pic50",
            "PDB_code",
            "release_year",
            "series",
            "template_strucid",
        ]
    ]
    si_df["compound_id"] = si_df["PDB_code"].str.split("_", expand=True)[1]
    si_df = si_df.drop(columns="PDB_code")
    si_df = si_df.rename(
        columns={
            "release_year": "pic50_date",
            "series": "binding_mode_class",
            "pde10_h_pic50": "pic50",
        }
    )
    si_df = si_df.join(pdb_proasis_df.set_index("PDB_code"), on="template_strucid")
    si_df["docking_folder"] = si_df["template_strucid"] + "_" + si_df["compound_id"]
    smiles_df = pd.read_csv("dgl_models_random/attentive_fp/2d_dataset_df.csv")
    smiles_df = smiles_df[["SRN", "SMILES"]]
    si_df = si_df.join(smiles_df.set_index("SRN"), on="SRN")
    si_df = si_df[
        [
            "SRN",
            "compound_id",
            "SMILES",
            "binding_mode_class",
            "pic50",
            "pic50_date",
            "DATESOLVED",
            "docking_folder",
        ]
    ]
    si_df = si_df.sort_values("compound_id")
    scenarios = Path("..").glob("dgl_models_*")
    for scenario in scenarios:
        scenario_name = str.join("_", scenario.stem.split("_")[2:])
        scenario_df = pd.read_csv(scenario / "dataset_df.csv")
        scenario_df = scenario_df[["SRN", "split"]]
        scenario_df = scenario_df.rename(columns={"split": scenario_name + "_split"})
        si_df = si_df.join(scenario_df.set_index("SRN"), on="SRN")

    scenarios = [
        "random",
        "temporal_2011",
        "temporal_2012",
        "temporal_2013",
        "aminohetaryl_c1_amide",
        "c1_hetaryl_alkyl_c2_hetaryl",
        "aryl_c1_amide_c2_hetaryl",
    ]

    table_maker = generate_metrics_table.TableMaker()
    for scenario in scenarios:
        prediction_df = table_maker.get_predictions_df(
            dgl_path=Path(f"dgl_models_{scenario}"),
            rf_score_path=Path(
                f"docking_template_date/first_iteration_{scenario}_split"
            ),
        )

        si_df["SRN"] = si_df["SRN"].str.split("-", expand=True)[0]
        prediction_df = prediction_df.join(si_df.set_index("SRN"), on="SRN")
        prediction_df = prediction_df.drop(columns="SRN")
        prediction_df.to_csv(f"raw_predictions_{scenario}.csv", index=False, sep="\t")

    si_df = si_df.drop(columns="SRN")

    # si_df.to_csv("table_s1.csv", index=False)


if __name__ == "__main__":
    main()
