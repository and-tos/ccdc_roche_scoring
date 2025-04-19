import argparse
from pathlib import Path

from ccdc import io, protein
from ccdc.docking import Docker

from ccdc_roche_scoring.docking import _dock


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="""
        Template docking standalone script.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument(
        "-l",
        "--ligands",
        help="SDF with ligands",
        required=True,
    )

    parser.add_argument(
        "-t",
        "--template",
        help="SDF with template",
        required=True,
    )

    parser.add_argument(
        "-p",
        "--protein",
        help="MOL2 with protein",
        required=True,
    )

    parser.add_argument(
        "-d",
        "--docking_folder",
        help="Folder where docking poses are saved",
        default="dock_dir",
        required=False,
    )

    parser.add_argument(
        "-n",
        "--ndocks",
        help="Number of poses generated per ligand.",
        default=3,
        type=int,
        required=False,
    )

    parser.add_argument(
        "-fr",
        "--flexible_residues",
        nargs="+",
        help="Residues that should be treates as flexible.",
        default=[],
        required=False,
    )

    return parser.parse_args()


def main():
    args = parse_args()
    ccdc_protein = protein.Protein.from_entry(io.EntryReader(args.protein)[0])
    native_ligand = io.MoleculeReader(args.template)[0]
    ligand_filename = Path(args.ligands)
    docking_folder = "."
    strucid = ccdc_protein.identifier
    ligand_id = "ligand"
    water_paths = None
    ccdc_template = io.MoleculeReader(args.template)[0]

    _dock(
        Docker(),
        Path(args.protein),
        ccdc_protein,
        native_ligand,
        ligand_filename,
        docking_folder,
        strucid,
        ligand_id,
        water_paths,
        ccdc_template,
        ccdc_template.copy(),
        False,
        True,
        ndocks=args.ndocks,
    )
    return


if __name__ == "__main__":
    main()
