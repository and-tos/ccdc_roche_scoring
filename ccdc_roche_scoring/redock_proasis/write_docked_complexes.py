from ccdc import io, protein, molecule
from pathlib import Path
from rf_statistics.database_utils import rdkit_library
from rdkit import Chem
import multiprocessing
from functools import partial
import click

# Pickle molecules contain properties
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


def combine_protein_ligand(ccdc_ligand, ccdc_protein):
    for a in ccdc_ligand.atoms:
        a.label = "_Z" + a.label

    ccdc_protein.add_ligand(ccdc_ligand)
    return ccdc_protein


def run_multiprocessing(nproc, iterable, function, **kwargs):
    parallel_func = partial(function, **kwargs)
    pool = multiprocessing.Pool(nproc)
    output = pool.map(parallel_func, iterable)
    return output


def parse_gold_pose(docking_dir):
    ligand_files = Path(docking_dir).glob("gold_soln_*.sdf")
    ligand_entries = io.EntryReader([str(f) for f in ligand_files])
    if ligand_entries:
        top_pose = max(
            ligand_entries,
            key=lambda ligand_entry: float(ligand_entry.attributes["Gold.PLP.Fitness"]),
        )
        ccdc_gold_entry = top_pose
        ccdc_protein = protein.Protein.from_file(str(docking_dir / "protein.mol2"))
        ccdc_protein = combine_protein_ligand(ccdc_gold_entry.molecule, ccdc_protein)
        ccdc_protein.attributes["docking_dir"] = docking_dir
        ccdc_gold_entry.identifier = ccdc_protein.identifier
        ccdc_gold_entry.attributes["docking_dir"] = docking_dir
        return ccdc_protein.to_string(), ccdc_gold_entry.to_string()
    else:
        return (None, None)


def parse_vina_pose(docking_dir):
    ligand_file = Path(docking_dir) / "vina_results.sdf"
    if ligand_file.is_file():
        ligand_entries = io.EntryReader(str(ligand_file))
        if ligand_entries:
            top_pose = ligand_entries[0]
            ccdc_vina_entry = top_pose
            ccdc_protein = protein.Protein.from_file(str(docking_dir / "protein.mol2"))
            ccdc_protein = combine_protein_ligand(
                ccdc_vina_entry.molecule, ccdc_protein
            )
            ccdc_vina_entry.identifier = ccdc_protein.identifier
            ccdc_vina_entry.attributes["docking_dir"] = docking_dir
            ccdc_protein.attributes["docking_dir"] = docking_dir
            return (ccdc_protein.to_string(), ccdc_vina_entry.to_string())
    else:
        return (None, None)


def parse_reference_pose(docking_dir):
    ligand_file = Path(docking_dir) / "ligand_native.sdf"
    if ligand_file.is_file():
        ligand = io.EntryReader(str(ligand_file))[0]
        ccdc_protein = protein.Protein.from_file(str(docking_dir / "protein.mol2"))
        ligand.identifier = ccdc_protein.identifier
        return ligand.to_string()


def write_gold_poses(docking_dirs, db):
    rdkit_gold_mols = []
    ccdc_gold_proteins = []
    ccdc_gold_ligands = []

    output = run_multiprocessing(32, docking_dirs, parse_gold_pose)
    for ccdc_gold_protein_string, ccdc_gold_ligand_string in output:
        if ccdc_gold_protein_string and ccdc_gold_ligand_string:
            try:
                ccdc_gold_ligand = molecule.Molecule.from_string(
                    ccdc_gold_ligand_string
                )
                ccdc_gold_ligands.append(ccdc_gold_ligand)

                ccdc_gold_protein = protein.Protein.from_string(
                    ccdc_gold_protein_string
                )
                ccdc_gold_proteins.append(ccdc_gold_protein)

                rdkit_gold_mols.append(
                    rdkit_library.rdkit_mol_from_ccdc_mol(ccdc_gold_protein)
                )
            except:
                continue

    with io.MoleculeWriter(f"gold_top_docked_complexes_{db}.csdsql") as w:
        for ccdc_gold_protein in ccdc_gold_proteins:
            w.write(ccdc_gold_protein)

    with io.MoleculeWriter(f"gold_top_docked_ligands_{db}.csdsql") as w:
        for ccdc_gold_ligand in ccdc_gold_ligands:
            w.write(ccdc_gold_ligand)

    rdkit_library.write_rdkit_library(
        rdkit_gold_mols, f"gold_top_docked_complexes_{db}_rdkit.p"
    )

    return


def write_vina_poses(docking_dirs, db):
    rdkit_vina_mols = []
    ccdc_vina_proteins = []
    ccdc_vina_ligands = []

    output = run_multiprocessing(32, docking_dirs, parse_vina_pose)
    for ccdc_vina_protein_string, ccdc_vina_ligand_string in output:
        if ccdc_vina_protein_string and ccdc_vina_ligand_string:
            ccdc_vina_ligand = molecule.Molecule.from_string(ccdc_vina_ligand_string)
            print(ccdc_vina_ligand.identifier)
            try:
                ccdc_vina_ligands.append(ccdc_vina_ligand)

                ccdc_vina_protein = protein.Protein.from_string(
                    ccdc_vina_protein_string
                )
                ccdc_vina_proteins.append(ccdc_vina_protein)

                rdkit_vina_mols.append(
                    rdkit_library.rdkit_mol_from_ccdc_mol(ccdc_vina_protein)
                )
            except:
                continue

    with io.MoleculeWriter(f"vina_top_docked_complexes_{db}.csdsql") as w:
        for ccdc_vina_protein in ccdc_vina_proteins:
            w.write(ccdc_vina_protein)

    with io.MoleculeWriter(f"vina_top_docked_ligands_{db}.csdsql") as w:
        for ccdc_vina_ligand in ccdc_vina_ligands:
            w.write(ccdc_vina_ligand)

    rdkit_library.write_rdkit_library(
        rdkit_vina_mols, f"vina_top_docked_complexes_{db}_rdkit.p"
    )
    return


def write_reference_poses(docking_dirs, db):
    output = run_multiprocessing(32, docking_dirs, parse_reference_pose)
    with io.MoleculeWriter(f"reference_ligands_{db}.csdsql") as w:
        for refernce_ligand_string in output:
            if refernce_ligand_string:
                refernce_ligand = molecule.Molecule.from_string(refernce_ligand_string)
                w.write(refernce_ligand)


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-i",
    "--input",
    help="Which type of files to collect.",
    type=click.Choice(["gold", "vina", "reference"], case_sensitive=False),
    required=True,
)
@click.option(
    "-db",
    "--db",
    help="Which type of files to collect.",
    type=click.Choice(["pub", "roche"], case_sensitive=False),
    required=True,
)
def main(input, db):
    docking_dir = "redocking_gold"
    if db == "roche":
        docking_dirs = Path(docking_dir).glob("?????_???")
    elif db == "pub":
        docking_dirs = Path(docking_dir).glob("????_???")
    if input == "gold":
        write_gold_poses(docking_dirs, db)
    elif input == "vina":
        write_vina_poses(docking_dirs, db)
    elif input == "reference":
        write_reference_poses(docking_dirs, db)
    return


if __name__ == "__main__":
    main()
