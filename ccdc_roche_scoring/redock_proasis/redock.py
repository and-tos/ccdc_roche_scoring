from pathlib import Path
import subprocess as sp
import click
import pandas as pd
import numpy as np
from ccdc import io, docking, protein
from vina import Vina

from ccdc_roche_scoring.redock_proasis.score_native_pose import (
    smallest_enclosing_cuboid,
)
from rf_statistics.structure_quality import StructureQuality


def is_water_atom(atom):
    if atom.label.startswith("HOH"):
        return True
    if atom.label.startswith("H") and not atom.label.startswith("HOH"):
        if atom.neighbours[0].label.startswith("HOH"):
            return True
    return False


def return_atom_type(atom):
    if is_water_atom(atom):
        return "water"
    if atom.label.startswith("_Z"):
        return "ligand"
    else:
        return "protein"


def remove_waters(ccdc_protein):
    contacts = ccdc_protein.contacts(path_length_range=(-1, 0))
    contact_dict = {
        "atom0": [],
        "atom1": [],
        "atom0_label": [],
        "atom1_label": [],
        "atom0_type": [],
        "atom1_type": [],
        "atom0_index": [],
        "atom1_index": [],
        "distance": [],
    }
    for contact in contacts:
        if contact.is_in_line_of_sight:
            contact_dict[f"distance"].append(contact.length)
            for cnt, contact_atom in enumerate(contact.atoms):
                contact_dict[f"atom{cnt}"].append(contact_atom)
                atom_type = return_atom_type(contact_atom)
                contact_dict[f"atom{cnt}_type"].append(atom_type)

                if contact_atom.label.startswith(
                    "H"
                ) and not contact_atom.label.startswith("HOH"):
                    contact_atom = contact_atom.neighbours[0]
                contact_dict[f"atom{cnt}_label"].append(contact_atom.label)
                contact_dict[f"atom{cnt}_index"].append(contact_atom.index)
    contact_df = pd.DataFrame(contact_dict)
    water_contact_df = contact_df[
        (contact_df["atom0_type"] == "water") | (contact_df["atom1_type"] == "water")
    ]
    water_contact_df = water_contact_df[
        ~(
            (water_contact_df["atom0_type"] == "water")
            & (water_contact_df["atom1_type"] == "water")
        )
    ]

    atom0_water = water_contact_df[water_contact_df["atom0_type"] == "water"]
    water_contact_df.loc[atom0_water.index, "water_index"] = atom0_water["atom0_index"]

    atom1_water = water_contact_df[water_contact_df["atom1_type"] == "water"]
    water_contact_df.loc[atom1_water.index, "water_index"] = atom1_water["atom1_index"]

    keep_atoms = []
    for water_index, group_df in water_contact_df.groupby("water_index"):
        atom_types = group_df["atom0_type"].to_list() + group_df["atom1_type"].to_list()
        protein_count = atom_types.count("protein")
        ligand_count = atom_types.count("ligand")
        if (ligand_count > 0 and protein_count > 0) or protein_count > 2:
            keep_atoms.append(ccdc_protein.atoms[int(water_index)])
    remove_atoms = []
    for a in ccdc_protein.heavy_atoms:
        if a.label.startswith("HOH"):
            if a not in keep_atoms:
                remove_atoms.append(a)
                remove_atoms.extend(list(a.neighbours))
    ccdc_protein.remove_atoms(remove_atoms)
    return ccdc_protein


def get_low_quality_ids():
    quality_df = pd.read_parquet(
        "/LoS/protonate3d_2024-12-03/pub_roche_project_annotation_2024-12-03.gz"
    )
    quality_df = quality_df.rename(columns={"bs_id": "identifier"})
    sq = StructureQuality(quality_df)
    low_q_bs_ids = sq.low_quality_bs_ids
    return low_q_bs_ids


def randomize_position(ligand_file, output_file):
    ligand = io.MoleculeReader(ligand_file)[0]
    ligand.rotate((1, 2, 3), 45)
    ligand.translate((0, 1, 2))
    for bond in ligand.bonds:
        if bond.is_rotatable and bond.sybyl_type != "am":
            ba0 = bond.atoms[0]
            ba1 = bond.atoms[1]
            na0 = [a for a in ba0.neighbours if a.atomic_symbol != "H" and a != ba1]
            na1 = [a for a in ba1.neighbours if a.atomic_symbol != "H" and a != ba0]
            if na0 and na1:
                ta0 = na0[0]
                ta1 = na1[0]
            else:
                continue
            ligand.set_torsion_angle(ta0, ba0, ba1, ta1, np.random.uniform(0, 180))
    with io.MoleculeWriter(output_file) as mol_writer:
        mol_writer.write(ligand)

    return ligand


def get_covalend_ids(df):
    return list(df[df["ligand_atom_is_warhead"]]["identifier"].unique())


class Redocker(object):
    def __init__(self, db) -> None:
        self.low_quality_ids = get_low_quality_ids()
        contacts_df = pd.read_parquet(
            f"/LoS/protonate3d_2024-12-03/{Path(db).stem}_rf_contacts.gzip"
        )
        self.covalent_ids = get_covalend_ids(contacts_df)

        self.cofactor_list = protein.Protein.known_cofactor_codes() + [
            "AMP",
            "ADP",
            "ATP",
            "GMP",
            "GDP",
            "GTP",
        ]

    def gold_docking(self, central_ligand):
        docker = docking.Docker()
        settings = docker.settings

        settings.make_absolute_file_names
        settings.add_ligand_file(self.random_conformer_file, ndocks=20)
        settings.add_protein_file(self.protein_file)
        settings.fitness_function = "plp"
        settings.autoscale = 5.0
        settings.early_termination = False
        settings.binding_site = settings.BindingSiteFromLigand(
            settings.proteins[0], central_ligand, 6
        )
        settings.reference_ligand_file = self.ligand_file
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
        docker.dock(str(self.dock_dir / "api_gold.conf"))

    def vina_docking(self):
        # Setup ligand for Vina
        lig_pdbqt = (self.dock_dir / "ligand.pdbqt").absolute()

        prepare_ligand_args = f"-i {self.random_conformer_file} -o {lig_pdbqt}".split()
        sp.check_output(
            [
                "ccdc_roche_scoring/bin/mk_prepare_ligand.py"
            ]
            + prepare_ligand_args
        )

        # Setup receptor for Vina
        rec_pdbqt = (self.dock_dir / "protein.pdbqt").absolute()
        prepare_receptor_args = (
            f"-r {self.protein_file} -o {rec_pdbqt} -U nphs_lps_waters"
        ).split()
        sp.check_output(
            ["/home/tosstora/ADFRsuite/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor"]
            + prepare_receptor_args
        )

        # Docking with Vina
        n_poses = 20
        v = Vina()
        v.set_receptor(str(rec_pdbqt))
        v.set_ligand_from_file(str(lig_pdbqt))

        # specify docking box
        ccdc_ligand = io.MoleculeReader(str(self.dock_dir / "ligand_native.sdf"))[0]
        center, side_length = smallest_enclosing_cuboid(ccdc_ligand)
        ligand_native_pdbqt = str(self.dock_dir / "ligand_native.pdbqt")
        prepare_ligand_args = f"-i {self.ligand_file} -o {ligand_native_pdbqt}".split()
        sp.check_output(
            [
                "ccdc_roche_scoring/bin/mk_prepare_ligand.py"
            ]
            + prepare_ligand_args
        )
        v.set_ligand_from_file(str(self.dock_dir / "ligand_native.pdbqt"))
        v.compute_vina_maps(center, side_length)

        v.randomize()
        v.dock(exhaustiveness=16)
        results_pdbqt = str((self.dock_dir / "vina_results.pdbqt").absolute())
        v.write_poses(results_pdbqt, n_poses, overwrite=True)
        results_sdf = str((self.dock_dir / "vina_results.sdf").absolute())
        sp.check_output(
            [
                "ccdc_roche_scoring/bin/mk_export.py",
                results_pdbqt,
                "-o",
                results_sdf,
            ]
        )

        vina_entries = io.EntryReader(results_sdf)
        new_vina_entries = []
        for cnt, ve in enumerate(vina_entries):
            ve.attributes["vina_score"] = v.energies(n_poses)[cnt][0]
            new_vina_entries.append(ve)
        with io.EntryWriter(results_sdf) as w:
            for ve in new_vina_entries:
                w.write(ve)

    def _redock_proasis_entry(self, entry_index, entry_name, proasis_csdsql):
        if not (entry_index or entry_name):
            raise click.UsageError("Must specify --entry_index or --entry_name")
        if entry_index:
            ccdc_entry = io.EntryReader(proasis_csdsql)[int(entry_index)]
        if entry_name:
            ccdc_entry = io.EntryReader(proasis_csdsql).entry(entry_name)
        entry_id = ccdc_entry.identifier

        if entry_id in self.low_quality_ids:
            print(entry_id + " is low quality.")
            return

        if entry_id in self.covalent_ids:
            print(entry_id + " is covalent ligand.")
            return

        ccdc_protein = protein.Protein.from_entry(ccdc_entry)

        central_ligand = [
            c for c in ccdc_protein.components if "_Z" in c.atoms[0].label
        ][0]

        ligand_name = central_ligand.atoms[0].residue_label[:3]
        if ligand_name in self.cofactor_list:
            print("ligand is cofactor")
            return

        non_standard_ligand_atoms = [
            a.atomic_symbol
            for a in central_ligand.heavy_atoms
            if a.atomic_symbol not in ["C", "N", "O", "S", "I", "Cl", "F", "Br", "P"]
        ]
        if non_standard_ligand_atoms:
            print("Ligand contains non-standard atoms: ", non_standard_ligand_atoms)
            return
        if 7 > len(central_ligand.heavy_atoms) or len(central_ligand.heavy_atoms) > 80:
            print(
                entry_id,
                " not docked. Heavy atom count:",
                len(central_ligand.heavy_atoms),
            )
            return

        self.dock_dir = Path(ccdc_protein.identifier)
        self.dock_dir.mkdir(exist_ok=True)
        self.protein_file = str((self.dock_dir / "protein.mol2").absolute())
        self.ligand_file = str((self.dock_dir / "ligand_native.sdf").absolute())
        self.random_conformer_file = str(
            (self.dock_dir / "ligand_random.sdf").absolute()
        )

        ccdc_protein = remove_waters(ccdc_protein)
        remove_atoms = [ccdc_protein.atom(a.label) for a in central_ligand.heavy_atoms]

        neighbours = []
        for a in remove_atoms:
            for n in a.neighbours:
                neighbours.append(n)

        remove_atoms = list(set(remove_atoms + neighbours))
        ccdc_protein.remove_atoms(remove_atoms)

        with io.MoleculeWriter(self.protein_file) as w:
            w.write(ccdc_protein)

        with io.MoleculeWriter(self.ligand_file) as w:
            w.write(central_ligand)
        random_conformer = randomize_position(
            self.ligand_file, str((self.dock_dir / "ligand_random.sdf").absolute())
        )
        # self.vina_docking()
        self.gold_docking(central_ligand)


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-i",
    "--entry_index",
    type=int,
    help="Entry index in the csdsql.",
    required=False,
)
@click.option(
    "-n",
    "--entry_name",
    type=str,
    help="Entry name in the csdsql.",
    required=False,
)
@click.option(
    "--proasis_csdsql",
    type=str,
    default="full_p2cq_pub_jul2023.csdsql",
    help="Entry index in the csdsql.",
    required=False,
)
def redock_proasis_entry(proasis_csdsql, entry_name=False, entry_index=False):
    if entry_index:
        print("Entry index: ", entry_index)
        redocker = Redocker(proasis_csdsql)
    elif entry_name:
        print("Entry name: ", entry_name)
        redocker = Redocker(proasis_csdsql)
    redocker._redock_proasis_entry(entry_index, entry_name, proasis_csdsql)
    return


def main():
    redock_proasis_entry()


if __name__ == "__main__":
    main()
