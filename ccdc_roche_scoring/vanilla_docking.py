#!/usr/bin/env python

########################################################################################################################

import argparse
from pathlib import Path

from ccdc import io, protein
from ccdc.docking import Docker

########################################################################################################################


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="""
        Execute Line of sight contact scripts.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument(
        "--input_ligands", help="SDF file with aligned ligands", default=False
    )

    parser.add_argument(
        "--native_ligand", help="SDF file with native ligand", default=False
    )

    parser.add_argument(
        "--docking_dir", help="Directory with docked ligands", default="."
    )

    parser.add_argument("-t", "--target", help="Target name.", default=None)

    return parser.parse_args()


def _dock(
    dry_receptor_file,
    target_protein: protein.Protein,
    native_ligand,
    ligand_filename,
    docking_folder,
    water_paths,
    reference_ligand_file=False,
):
    docker = Docker()
    settings = docker.settings
    settings.fitness_function = "plp"
    settings.save_lone_pairs = False

    settings.add_protein_file(str(dry_receptor_file.resolve()))
    native_ligand.add_hydrogens(mode="missing")
    settings.binding_site = settings.BindingSiteFromLigand(
        target_protein, native_ligand, 10.0
    )
    for water_path in water_paths:
        settings.add_water_file(str(water_path))
    ndocks = 3

    if reference_ligand_file:
        settings.reference_ligand_file = str(reference_ligand_file.resolve())

    settings.add_ligand_file(str(ligand_filename.resolve()), ndocks=ndocks)
    settings.output_file = str(docking_folder / Path(f"./docked_ligands.sdf"))
    settings.write_options = [
        "NO_LINK_FILES",
        "NO_RNK_FILES",
        "NO_PLP_MOL2_FILES",
        "NO_BESTRANKING_LST_FILE",
        "NO_GOLD_LIGAND_MOL2_FILE",
        "NO_LOG_FILES",
        "NO_FIT_PTS_FILES",
        "NO_SEED_LOG_FILE",
        "NO_GOLD_LIGAND_MOL2_FILE",
        "NO_GOLD_PROTEIN_MOL2_FILE",
        "NO_PID_FILE",
    ]
    gold_conf = docking_folder / Path(f"gold.conf")
    settings.make_absolute_file_names(str(gold_conf.resolve()))
    settings.write(str(gold_conf.resolve()))

    docker = Docker()
    docker.settings = docker.settings.from_file(str(gold_conf.resolve()))
    docker.dock(file_name=str(gold_conf.resolve()))
    print("GOLD docking finished.")

    return


def main():
    args = parse_args()
    docking_home = args.docking_dir
    ccdc_protein = protein.Protein.from_file(args.target)
    ccdc_protein.assign_bond_types()
    native_ligand = io.MoleculeReader(args.native_ligand)[0]
    water_paths = Path(docking_home).glob("water*.mol2")
    docking_dir = Path(docking_home)
    if not docking_dir.is_dir():
        docking_dir.mkdir()
    _dock(
        Path(args.target),
        ccdc_protein,
        native_ligand,
        Path(args.input_ligands),
        docking_dir,
        water_paths,
    )

    return


if __name__ == "__main__":
    main()
