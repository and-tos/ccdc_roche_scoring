from ccdc import io, protein
from ccdc.docking import Docker
from ccdc_roche_scoring.docking import _dock

import click
from pathlib import Path


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-l",
    "--ligands",
    help="SDF with ligands",
    required=True,
)
@click.option(
    "-p",
    "--protein_file",
    help="MOL2 with protein",
    required=True,
)
@click.option(
    "-nl",
    "--native_ligand_file",
    help="SDF with native ligand. Required for binding site definition.",
    required=True,
)
@click.option(
    "-t",
    "--template",
    is_flag=True,
    help="If specified, native ligand will be used for template guided docking. Otherwise, docking without template constraints will be performed.",
    required=False,
)
@click.option(
    "-nd",
    "--ndocks",
    help="Number of poses generated per ligand.",
    default=3,
    type=int,
    required=False,
)
@click.option(
    "-fr",
    "--flexible_residues",
    multiple=True,
    default=[],
    help="Residues that should be treated as flexible. Currently not exposed.",
    required=False,
)
@click.option(
    "-gc",
    "--gold_conf",
    default=None,
    help="Path to GOLD config file.",
    required=False,
)
def gold_dock(
    ligands,
    protein_file,
    native_ligand_file,
    ndocks,
    template=True,
    flexible_residues=[],
    gold_conf="",
):
    ccdc_protein = protein.Protein.from_entry(io.EntryReader(protein_file)[0])
    native_ligand = io.MoleculeReader(native_ligand_file)[0]
    ligand_filename = Path(ligands)
    docking_folder = "."
    strucid = ccdc_protein.identifier
    if template:
        scaffold = io.MoleculeReader(native_ligand_file)[0]
        strict_scaffold = io.MoleculeReader(native_ligand_file)[0]
    else:
        scaffold = False
        strict_scaffold = False
    kwargs = {"flexible_residues": flexible_residues, "gold_conf": gold_conf}
    _dock(
        Docker(),
        Path(protein_file),
        ccdc_protein,
        native_ligand,
        ligand_filename,
        docking_folder,
        strucid,
        ligand_identifier="ligand",
        water_paths=None,
        scaffold=scaffold,
        strict_scaffold=strict_scaffold,
        ndocks=ndocks,
        **kwargs
    )
    return


if __name__ == "__main__":
    gold_dock()
