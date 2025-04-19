#!/usr/bin/env python

########################################################################################################################
import warnings

warnings.filterwarnings("ignore", ".*deprecated.*", DeprecationWarning, ".*", 0)

import argparse
import os
import subprocess as sp
import sys
from pathlib import Path
from collections import defaultdict

import pandas as pd

# import this before csd api
from rdkit import Chem

from ccdc_roche_scoring import template_docking_mcs

from ccdc import descriptors, entry, io, protein, search
from ccdc.docking import Docker
from rf_statistics.los_descriptors import _cut_out_binding_site_by_distance


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
        "--docking_dir", help="Directory with docked ligands", default=False
    )

    parser.add_argument("-t", "--target", help="Target name.", default="default")
    parser.add_argument("-gc", "--gold_conf", help="GOLD config file.", default=False)

    parser.add_argument(
        "-fr",
        "--flexible_residues",
        nargs="+",
        help="Residues that should be treates as flexible.",
        default=[],
    )

    return parser.parse_args()


class GeometryFlag(object):
    def __init__(
        self, smarts, atom_index_1, atom_index_2, min_distance, has_to_be_los=True
    ):
        self.smarts = smarts
        self.substructure = search.SMARTSSubstructure(smarts)
        self.atom_index_1 = atom_index_1
        self.atom_index_2 = atom_index_2
        self.min_distance = min_distance
        self.has_to_be_los = has_to_be_los


class ConformerPruner(object):
    def __init__(self):
        self.flagged_geometries = self._flagged_geometries()

    def _flagged_geometries(self):
        flagged_geometries = []
        flagged_geometries.append(
            GeometryFlag("N(H)-C(=O)-!@[#6]@[#6]-!@N(H)", 1, 7, 2.2)
        )
        flagged_geometries.append(
            GeometryFlag("N(H)-C(=O)-!@[#6](~[#7D2H0])@[#6]-!@N(H)", 1, 8, 2.5)
        )
        flagged_geometries.append(
            GeometryFlag("N(H)-C(=O)-!@[#6](~[#7D2H0])@[#6]-!@N(H)", 3, 6, 2.9)
        )
        flagged_geometries.append(
            GeometryFlag("N(H)-C(=O)-!@[#6](~[#7D2H0])@[#6]-!@N(H)", 1, 7, 3.0)
        )
        flagged_geometries.append(
            GeometryFlag("N(H)-C(=O)-!@[#6]@[#6]-!@N(H)", 0, 6, 3.0)
        )
        flagged_geometries.append(
            GeometryFlag("N(H)-C(=O)-!@[#6](~[#6H1])~[#6]~[#7D2H0]", 3, 7, 3.8)
        )
        flagged_geometries.append(
            GeometryFlag("N([H])-C(=O)-!@[#6]@[#7]([H])", 1, 6, 2.1)
        )
        flagged_geometries.append(
            GeometryFlag("C(=O)-!@[#6]@[#6]-!@N-C(=O)", 1, 6, 2.5)
        )
        flagged_geometries.append(
            GeometryFlag("[#6](=O)-!@[c]([cD2])@[c]([cD2])-O-H", 1, 6, 3)
        )
        flagged_geometries.append(
            GeometryFlag("[#6](=O)-[c]([cD3][#6])@[c]([cD2])-O-H", 1, 7, 3)
        )
        flagged_geometries.append(
            GeometryFlag("[NH1][C](=O)!@[#6]~[#6][N-]", 2, 5, 3.5)
        )
        flagged_geometries.append(
            GeometryFlag("[#7D2H0]~[#6](~[#6])-!@[#6](~[#6])~[#7D2H0]", 0, 5, 3)
        )
        flagged_geometries.append(
            GeometryFlag("[#7D2H0]~[#6](~[#6][H])-!@[#6](~[#6][H])~[#7D2H0]", 3, 6, 2.5)
        )
        flagged_geometries.append(GeometryFlag("C(=O)-[#6]@[#7D2H0]", 1, 3, 3.0, False))
        flagged_geometries.append(GeometryFlag("C(=O)-N-[#6]@[#7D2H0]", 1, 4, 3.0))

        return flagged_geometries

    def is_bad_conformer(self, ligand):
        for flagged_geometry in self.flagged_geometries:
            searcher = search.SubstructureSearch()
            searcher.add_substructure(flagged_geometry.substructure)
            hits = searcher.search(ligand)
            for hit in hits:
                match_atoms = hit.match_atoms()
                is_los = match_atoms[flagged_geometry.atom_index_1].is_in_line_of_sight(
                    match_atoms[flagged_geometry.atom_index_2]
                )
                distance = descriptors.MolecularDescriptors.atom_distance(
                    match_atoms[flagged_geometry.atom_index_1],
                    match_atoms[flagged_geometry.atom_index_2],
                )
                if distance < flagged_geometry.min_distance:
                    if flagged_geometry.has_to_be_los:
                        if is_los:
                            return True
                    else:
                        return True
        return False


def update_gold_conf(gold_conf, water_paths=False, fixed_bonds=""):
    f = open(gold_conf, "r")
    newdata = f.read()
    f.close()

    newdata = newdata.replace(
        "internal_ligand_h_bonds = 0", "internal_ligand_h_bonds = 1"
    )
    newdata = newdata.replace("rms_tolerance = 1.5", "rms_tolerance = 0.5")
    newdata = newdata.replace("save_lone_pairs = 1", "save_lone_pairs = 0")

    newdata = newdata.replace("pt_crosswt = 0", "pt_crosswt = 95")
    newdata = newdata.replace("allele_mutatewt = 0", "allele_mutatewt = 95")
    newdata = newdata.replace("migratewt = 0", "migratewt = 10")

    # tor_lib_2020.tordist should be the default starting with CSDS 2023.3 update
    # newdata = newdata.replace("tor_lib_2020.tordist", "gold.tordist")

    newdata = newdata.replace("solvate_all = 1", f"solvate_all = 1\n" + fixed_bonds)
    newdata = newdata.replace(
        "param_file = DEFAULT",
        f"param_file = gold_custom.params",
    )
    newdata = newdata.replace(
        "rescore_param_file = DEFAULT",
        f"rescore_param_file = gold_custom.params",
    )

    if water_paths:
        newdata = newdata + "  WATER DATA\n"
        for water_file in water_paths:
            newdata = newdata + f"water 1 toggle trans_spin 0.5 {water_file}\n"

    f = open(gold_conf, "w")
    f.write(newdata)
    f.close()


def _mcs_metrics(docked_ligand_entry, scaffold):
    docked_ligand = docked_ligand_entry.molecule
    docked_ligand.set_formal_charges()
    rdkit_ligand = Chem.MolFromMolBlock(docked_ligand.to_string("sdf"), removeHs=False)
    docked_ligand = docked_ligand.from_string(Chem.MolToMolBlock(rdkit_ligand))

    docked_ligand.remove_unknown_atoms()
    mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    mcs_searcher.settings.check_bond_type = False
    mcs_searcher.settings.ignore_hydrogens = True
    mcs_searcher.settings.check_bond_count = False
    temp_docked_ligand = docked_ligand.copy()
    temp_scaffold = scaffold.copy()
    mcs_atoms_docked_template = mcs_searcher.search(
        temp_docked_ligand, temp_scaffold, search_step_limit=1000000
    )[0]

    # Keep only MCS part to allow RMSD calculation to consider symmetry
    ligand_mcs_atoms = [mcs_pair[0] for mcs_pair in mcs_atoms_docked_template]
    remove_atoms = [a for a in temp_docked_ligand.atoms if a not in ligand_mcs_atoms]
    temp_docked_ligand.remove_atoms(remove_atoms)
    temp_docked_ligand.remove_hydrogens()

    scaffold_mcs_atoms = [mcs_pair[1] for mcs_pair in mcs_atoms_docked_template]
    remove_atoms = [a for a in temp_scaffold.atoms if a not in scaffold_mcs_atoms]
    temp_scaffold.remove_atoms(remove_atoms)
    temp_scaffold.remove_hydrogens()

    rmsd = descriptors.MolecularDescriptors.rmsd(
        temp_docked_ligand,
        temp_scaffold,
        # atoms=mcs_atoms_docked_template,
        exclude_hydrogens=True,
        with_symmetry=True,
    )
    docked_ligand_entry.attributes["RMSD_to_mcs"] = rmsd
    abs_mcs_size = len(scaffold.heavy_atoms)
    rel_mcs_size = abs_mcs_size / len(docked_ligand.heavy_atoms)
    docked_ligand_entry.attributes["rel_mcs_size"] = rel_mcs_size
    docked_ligand_entry.attributes["abs_mcs_size"] = abs_mcs_size
    return docked_ligand_entry, docked_ligand


def _select_best_soln(
    docking_folder, scaffold, row, srn, pdb_id, attributes, template_filename
):
    for docking_subjob_dir in [dir for dir in docking_folder.glob("*") if dir.is_dir()]:
        docked_ligand_files = docking_subjob_dir.glob(f"gold_soln_input_*_*_*.sdf")

        best_soln = entry.Entry()
        best_soln_fitness = -999
        for docked_ligand_file in docked_ligand_files:
            docked_ligand_entries = []
            docked_ligand_file = docked_ligand_file.resolve()

            with io.EntryReader(str(docked_ligand_file)) as rdr:
                for docked_ligand_entry in rdr:
                    docked_ligand_entry, docked_ligand = _mcs_metrics(
                        docked_ligand_entry,
                        scaffold,
                    )
                    docked_ligand_entry.attributes["template_strucid"] = row[
                        "template_strucid"
                    ]
                    docked_ligand_entry.attributes[
                        "rel_mcs_size_to_native_ligand"
                    ] = row["rel_mcs_size_to_native_ligand"]
                    docked_ligand_entry.attributes[
                        "tanimoto_similiarity_to_native_ligand"
                    ] = row["similarity"]
                    docked_ligand_entry.attributes["is_decoy"] = False
                    docked_ligand_entry.attributes.update(attributes)
                    docked_ligand_entries.append(docked_ligand_entry)
                    gold_rescore_fitness = float(
                        docked_ligand_entry.attributes["Gold.PLP.Fitness"]
                    )
                    if (
                        gold_rescore_fitness
                        > best_soln_fitness
                        # and not ConformerPruner().is_bad_conformer(
                        #     docked_ligand_entry.molecule
                        # )
                    ):
                        best_soln = docked_ligand_entry
                        best_soln_file = docked_ligand_file
                        best_soln_fitness = gold_rescore_fitness

            with io.EntryWriter(str(docked_ligand_file)) as wr:
                for docked_ligand_entry in docked_ligand_entries:
                    wr.write(docked_ligand_entry)

        shape_similarity = sp.check_output(
            [
                sys.executable,
                Path(__file__).parent / "shape_similarity.py",
                "--template_ligand",
                template_filename,
                "--docked_ligand",
                str(best_soln_file),
            ],
            env=os.environ,
        )

        shape_similarity = float(shape_similarity.decode("utf-8").split("\n")[0])
        print("shape", shape_similarity)

        best_soln.attributes["is_decoy"] = False
        best_soln.attributes["TanimotoCombo"] = shape_similarity
        best_soln.attributes.update(attributes)

        best_soln_mol = best_soln.molecule
        best_soln_mol.set_formal_charges()
        rdkit_ligand = Chem.MolFromMolBlock(
            best_soln_mol.to_string("sdf"), removeHs=False
        )
        best_soln_mol = best_soln_mol.from_string(Chem.MolToMolBlock(rdkit_ligand))
        best_soln_mol.remove_unknown_atoms()

        new_best_soln = entry.Entry.from_molecule(best_soln_mol)
        new_best_soln.attributes = best_soln.attributes
        new_best_soln.attributes["docking_folder"] = str(docking_subjob_dir.absolute())

        best_soln_file = str(
            Path(docking_folder) / Path(f"best_soln_{pdb_id}_{srn}.sdf")
        )
        new_best_soln.identifier = new_best_soln.attributes["SRN"]
        with io.EntryWriter(best_soln_file) as w:
            w.write(new_best_soln)

        gold_conf = str(list(Path(docking_subjob_dir).glob("*.conf"))[0])
        _write_best_soln_pocket(gold_conf, best_soln_file)
    return True


def _select_best_decoy(
    docking_folder, mcs_searcher, scaffold, srn, pdb_id, reference_ligand_file
):
    docked_ligand_files = docking_folder.glob(f"gold_soln_docking_input_*_*_*.mol2")

    best_decoy = entry.Entry()
    best_decoy_rmsd = 0

    for docked_ligand_file in docked_ligand_files:
        docked_ligand_entries = []
        docked_ligand_file = docked_ligand_file.resolve()
        with io.EntryReader(str(docked_ligand_file)) as rdr:
            for docked_ligand_entry in rdr:
                docked_ligand_entry, docked_ligand = _mcs_metrics(
                    docked_ligand_entry, scaffold
                )
                rmsd_to_mcs = docked_ligand_entry.attributes["RMSD_to_mcs"]
                docked_ligand_entry.attributes["template_strucid"] = pdb_id
                docked_ligand_entries.append(docked_ligand_entry)

                # if docked_ligand_entry.attributes['RMSD_to_mcs'] > best_decoy_rmsd:
                #     # best_decoy = docked_ligand_entry
                #     best_decoy_rmsd = docked_ligand_entry.attributes['RMSD_to_mcs']
                #     # best_decoy_file = docked_ligand_file

        shape_similarity = sp.check_output(
            [
                sys.executable,
                Path(__file__).parent / "shape_similarity.py",
                "--template_ligand",
                str(reference_ligand_file),
                "--docked_ligand",
                str(docked_ligand_file),
            ],
            env=os.environ,
        )

        shape_similarity = float(shape_similarity.decode("utf-8").split("\n")[0])
        print("shape", shape_similarity)

        # set attributes for decoy
        with io.EntryReader(str(reference_ligand_file)) as rdr:
            attributes = rdr[0].attributes
        docked_ligand_entry.attributes.update(attributes)
        docked_ligand_entry.attributes["RMSD_to_mcs"] = rmsd_to_mcs
        docked_ligand_entry.attributes["TanimotoCombo_to_reference"] = shape_similarity
        docked_ligand_entry.attributes["is_decoy"] = True

        with io.EntryWriter(docked_ligand_file) as wr:
            for docked_ligand_entry in docked_ligand_entries:
                wr.write(docked_ligand_entry)

    # shape_similarity = sp.check_output(
    #     [sys.executable, '/rf_scoring/shape_similarity.py',
    #      '--template_ligand', str(reference_ligand_file),
    #      '--docked_ligand', str(best_decoy_file)], env=os.environ)
    #
    # shape_similarity = float(shape_similarity.decode('utf-8').split('\n')[0])
    # print('shape', shape_similarity)
    #
    # # set attributes for decoy
    # with io.EntryReader(str(reference_ligand_file)) as rdr:
    #     attributes = rdr[0].attributes
    # best_decoy.attributes.update(attributes)
    # best_decoy.attributes['TanimotoCombo_to_reference'] = shape_similarity
    # best_decoy.attributes['is_decoy'] = True
    #
    # # write molecule with RDKit because CSD messes up formal charges of carboxylate groups from GOLD output
    # best_decoy_mol = best_decoy.molecule
    # best_decoy_mol.kekulize()
    # best_decoy_mol.remove_unknown_atoms()
    #
    # new_best_decoy = entry.Entry.from_molecule(best_decoy_mol)
    # new_best_decoy.attributes = best_decoy.attributes
    #
    # with io.EntryWriter(str(Path(docking_folder) / Path(f'best_decoy_{pdb_id}_{srn}.sdf'))) as w:
    #     w.write(new_best_decoy)
    return True


def _protein_preparation(pdb, dry_receptor_file, pdb_id, target):
    protein_res_labels = [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLU",
        "GLN",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HOH",
        "",
    ]

    target_protein = protein.Protein.from_file(str(pdb))

    water_atoms = [
        a
        for a in target_protein.atoms
        if a.protein_atom_type == "Water"
        and is_important_water(a, target_protein, target) == False
    ]
    ligand_atoms = [
        a
        for a in target_protein.atoms
        if ("LIG1" in a.residue_label or a.residue_label[:3] not in protein_res_labels)
    ]

    target_protein.remove_atoms(ligand_atoms + water_atoms)
    target_protein.standardise_aromatic_bonds()
    target_protein.standardise_delocalised_bonds()

    if not dry_receptor_file.parent.is_dir():
        dry_receptor_file.parent.mkdir()

    water_paths = list(dry_receptor_file.parent.glob("*water*.mol2"))

    if len(target_protein.waters) > 0 and len(water_paths) == 0:
        # water_paths = []
        for cnt, water in enumerate(target_protein.waters):
            water_path = (
                dry_receptor_file.parent / Path(f"{pdb_id}_water_{cnt}.mol2")
            ).resolve()
            water_paths.append(water_path)
            with io.MoleculeWriter(water_path) as w_wr:
                water.add_hydrogens()
                w_wr.write(water)
    target_protein.remove_all_waters()

    # save as pdb first to ensure correct format
    temp_pdb = str(dry_receptor_file.parent / Path(dry_receptor_file.stem + ".pdb"))
    with io.MoleculeWriter(temp_pdb) as p_wr:
        p_wr.write(target_protein)

    target_protein = protein.Protein.from_file(temp_pdb)
    target_protein.add_hydrogens()
    with io.MoleculeWriter(str(dry_receptor_file)) as p_wr:
        p_wr.write(target_protein)

    return water_paths, target_protein


def is_important_water(atom, protein, target="pde-10"):
    if target == "pde-10" :
        distance_searcher = descriptors.MolecularDescriptors.AtomDistanceSearch(protein)
        close_atoms = distance_searcher.atoms_within_range(atom.coordinates, 3.7)

    else:
        return False

    if target == "pde-10":
        tyr_contact = False
        gln_contact = False
        trp_contact = False
        asp1_contact = False
        asp2_contact = False

        for at in close_atoms:
            res_label = at.residue_label[0:3]
            atom_label = at.label
            if res_label == "TYR" and atom_label == "OH":
                tyr_contact = True
                continue

            if res_label == "TRP" and atom_label == "NE1":
                trp_contact = True
                continue

            if res_label == "GLN" and atom_label == "OE1":
                gln_contact = True
                continue

            if res_label == "ASP" and atom_label == "O":
                asp1_contact = True
                continue

            if res_label == "ASP" and atom_label == "OD1":
                asp2_contact = True
                continue

        if gln_contact and trp_contact and tyr_contact:
            return True
        elif tyr_contact and asp1_contact and asp2_contact:
            return True

        else:
            return False


def _fix_rotatable_bond(mcs_scaffold_bond, strict_scaffold):
    """
    Returns True if rotatable bond should be fixed in docking.
    :param mcs_scaffold_bond:
    :param strict_scaffold:
    :return:
    """
    if template_docking_mcs.atom_has_open_valency(
        mcs_scaffold_bond.atoms[0]
    ) or template_docking_mcs.atom_has_open_valency(mcs_scaffold_bond.atoms[1]):
        return False
    if template_docking_mcs.is_in_ring_with_open_valency(
        mcs_scaffold_bond.atoms[0], strict_scaffold
    ):
        return False
    if template_docking_mcs.is_in_ring_with_open_valency(
        mcs_scaffold_bond.atoms[1], strict_scaffold
    ):
        return False
    return True


class GOLDSetup(object):
    def __init__(
        self,
        dry_receptor_file,
        native_ligand,
        target_protein,
        reference_ligand_file=None,
    ) -> None:
        self.dry_receptor_file = dry_receptor_file
        self.native_ligand = native_ligand
        self.target_protein = target_protein
        self.reference_ligand_file = reference_ligand_file

    def default_gold_settings(
        self,
        docker,
        docking_folder,
        flexible_residues,
        lig_cnt,
        ligand_identifier,
        ndocks,
        pdb_id,
        ligand_mol,
    ):
        settings = docker.Settings()
        settings.fitness_function = "plp"

        settings.diverse_solutions = True, 1, 0.5

        settings.add_protein_file(str(self.dry_receptor_file.resolve()))
        self.native_ligand.add_hydrogens(mode="missing")
        settings.binding_site = settings.BindingSiteFromLigand(
            self.target_protein, self.native_ligand, 10.0
        )

        for flexible_residue in flexible_residues:
            res = [
                r
                for r in settings.proteins[0].residues
                if flexible_residue in r.identifier
            ][0]
            flexible_residues.append(res)
        for flexible_residue in flexible_residues:
            rl = settings.RotamerLibrary(settings.protein_files[0], flexible_residue)
            rl.add_default_rotamers()
            settings.add_rotamer_library(settings.proteins[0], rl)

        settings.autoscale = 1  # in percent, lower is faster

        if self.reference_ligand_file:
            settings.reference_ligand_file = str(self.reference_ligand_file.resolve())

        settings.output_directory = str(".")
        settings.output_file = f"{pdb_id}_docked_ligands_{ligand_identifier}.sdf"

        # Setup ligand
        docking_subjob_dir = Path(docking_folder) / f"dock_{lig_cnt}"
        docking_subjob_dir.mkdir(exist_ok=True)
        input_ligand_filename = docking_subjob_dir / "input_ligand.sdf"
        with io.MoleculeWriter(str(input_ligand_filename)) as w:
            w.write(ligand_mol)
        settings.clear_ligand_files()
        settings.add_ligand_file(str(input_ligand_filename.resolve()), ndocks=ndocks)

        self.gold_conf_file = docking_subjob_dir / Path(f"gold.conf")
        settings.make_absolute_file_names(str(self.gold_conf_file.resolve()))
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
        settings._settings.set_flip_planar_N(False)
        return settings

    def from_conf_file(
        self, gold_conf_file, ligand_mol, docking_folder, ndocks, lig_cnt=0
    ):
        docker = Docker()
        settings = docker.settings.from_file(str(gold_conf_file))
        docking_subjob_dir = Path(docking_folder) / f"dock_{lig_cnt}"
        docking_subjob_dir.mkdir(exist_ok=True)
        input_ligand_filename = docking_subjob_dir / "input_ligand.sdf"
        with io.MoleculeWriter(str(input_ligand_filename)) as w:
            w.write(ligand_mol)
        settings.clear_ligand_files()
        settings.add_ligand_file(str(input_ligand_filename.resolve()), ndocks=ndocks)
        settings.output_directory = str(".")
        settings.output_file = f"docking_soln.sdf"
        self.gold_conf_file = docking_subjob_dir / Path(f"gold.conf")
        return settings

    def add_template_docking_settings(
        self,
        settings,
        ligand_mol,
        water_paths,
        scaffold=None,
        strict_scaffold=None,
    ):
        """SETTINGS for template docking"""
        if scaffold:
            scaffold.add_hydrogens(mode="missing")
            scaffold.assign_bond_types()
            settings.add_constraint(
                settings.TemplateSimilarityConstraint("all", scaffold, weight=40)
            )

        if strict_scaffold:
            _strict_scaffold = strict_scaffold.copy()
            settings.add_constraint(
                settings.ScaffoldMatchConstraint(_strict_scaffold, weight=1000)
            )

            mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
            mcs_searcher.settings.ignore_hydrogens = True
            mcs_searcher.settings.check_bond_type = False
            if ligand_mol.largest_ring_size > 20:
                mcs_atoms, mcs_bonds = template_docking_mcs.macrocycle_mcs(
                    ligand_mol, _strict_scaffold, mcs_searcher
                )
            else:
                mcs_atoms, mcs_bonds = mcs_searcher.search(ligand_mol, _strict_scaffold)

            # fix bonds
            for mcs_ligand_bond, mcs_scaffold_bond in mcs_bonds:
                if mcs_ligand_bond.is_rotatable:
                    if _fix_rotatable_bond(mcs_scaffold_bond, _strict_scaffold):
                        settings.add_specific_fixed_rotatable_bond(mcs_ligand_bond)

            _strict_scaffold.remove_hydrogens()

        settings.write(str(self.gold_conf_file.resolve()))

        update_gold_conf(self.gold_conf_file, water_paths)

        return self.gold_conf_file


def _dock(
    docker,
    dry_receptor_file,
    target_protein: protein.Protein,
    native_ligand,
    ligand_filename,
    docking_folder,
    pdb_id,
    ligand_identifier,
    water_paths,
    scaffold=False,
    strict_scaffold=False,
    reference_ligand_file=False,
    ndocks=3,
    gold_conf=False,
    flexible_residues: list = [],
    **kwargs,
) -> None:
    for lig_cnt, ligand_mol in enumerate(io.MoleculeReader(str(ligand_filename))):
        gold_setup = GOLDSetup(
            dry_receptor_file, native_ligand, target_protein, reference_ligand_file
        )

        if gold_conf:
            gold_settings = gold_setup.from_conf_file(
                gold_conf, ligand_mol, docking_folder, ndocks
            )
        else:
            gold_settings = gold_setup.default_gold_settings(
                docker,
                docking_folder,
                flexible_residues,
                lig_cnt,
                ligand_identifier,
                ndocks,
                pdb_id,
                ligand_mol,
            )

        gold_conf = gold_setup.add_template_docking_settings(
            gold_settings,
            ligand_mol,
            water_paths,
            scaffold,
            strict_scaffold,
        )

        docker = Docker()
        docker.settings = docker.settings.from_file(str(gold_conf.resolve()))
        docker.dock(file_name=str(gold_conf.resolve()))
        print("GOLD docking finished.")

    return


def _mcs_templates_df(native_ligand_entries, ligand_mol, series_template_strucids=None):
    """
    Find maximum common substructure. Eliminate ring and hybridization mismatch atoms.
    :param native_ligand_entries:
    :param mcs_searcher:
    :param ligand_mol:
    :return:
    """
    templates_df = pd.DataFrame()
    templates = defaultdict(list)
    for cnt, native_ligand_entry in enumerate(native_ligand_entries):
        if (
            series_template_strucids is not None
            and native_ligand_entry.attributes["STRUCID"]
            not in series_template_strucids
        ):
            continue

        mcs = template_docking_mcs.MCS(ligand_mol, native_ligand_entry.molecule)
        # scaffold = mcs.mcs_scaffold
        scaffold = mcs.return_mcs_scaffold(
            partial_ring_matches_allowed=False, ignore_hydrogens=True
        )[0]
        mcs_atoms = mcs.mcs_atoms
        mcs_atom_labels = [
            (mcs_pair[0].label, mcs_pair[1].label)
            for mcs_pair in mcs_atoms
            if mcs_pair[0] in scaffold.atoms
        ]

        mcs_size = len(scaffold.heavy_atoms)
        templates["template_strucid"].append(native_ligand_entry.attributes["STRUCID"])
        templates["abs_mcs_size"].append(mcs_size)
        temp_ligand = ligand_mol.copy()
        templates["rel_mcs_size_to_ligand"].append(
            mcs_size / len(temp_ligand.heavy_atoms)
        )
        templates["rel_mcs_size_to_native_ligand"].append(
            mcs_size / len(native_ligand_entry.molecule.heavy_atoms)
        )
        templates["native_ligand"].append(native_ligand_entry)
        templates["scaffold"].append(scaffold)
        templates["mcs_object"].append(mcs)
        templates["mcs_atom_labels"].append(mcs_atom_labels)
        templates_df = (
            pd.DataFrame(templates)
            .sort_values(by="abs_mcs_size", ascending=False)
            .reset_index(drop=True)
        )

    # Tanimoto similarity for highest MCS compounds
    searcher = search.SimilaritySearch(ligand_mol)
    for index, row in templates_df.iterrows():
        similarity = searcher.search_molecule(row["native_ligand"].molecule).similarity
        templates_df.loc[index, "similarity"] = similarity

    templates_df.loc[templates_df["rel_mcs_size_to_ligand"] < 0.35, "similarity"] = 0
    templates_df = templates_df.sort_values(
        by=["abs_mcs_size", "similarity"], ascending=False
    ).reset_index(drop=True)
    print(templates_df)
    return templates_df


def return_mcs_atoms_of_identical_mols(mol1, mol2):
    mol1_noh = mol1.copy()
    mol1_noh.remove_hydrogens()
    substructure = search.MoleculeSubstructure(mol1_noh)
    searcher = search.SubstructureSearch()
    searcher.add_substructure(substructure)
    hit = searcher.search(mol2, max_hit_structures=1)[0]
    mol2_match_atoms = hit.match_atoms()

    hit = searcher.search(mol1, max_hit_structures=1)[0]
    mol1_match_atoms = hit.match_atoms()
    mcs_atoms = tuple(zip(mol1_match_atoms, mol2_match_atoms))
    return mcs_atoms


def _write_starting_ligand_macrocycle(
    ligand_filename,
    docking_folder,
    template_filename,
    mcs,
):
    "Open macrocycle for faster conformer generation, then close again."
    ring_mcs_bonds = [b for b in mcs.mcs_bonds if b[0].is_cyclic]
    overall_largest_ring_mcs_bond = ring_mcs_bonds[0]
    overall_largest_ring = max(
        overall_largest_ring_mcs_bond[0].rings, key=lambda ring: len(ring.atoms)
    )
    for ring_bond_pair in ring_mcs_bonds:
        ligand_ring_bond = ring_bond_pair[0]
        largest_ring = max(ligand_ring_bond.rings, key=lambda ring: len(ring.atoms))
        if len(largest_ring.atoms) > len(overall_largest_ring.atoms):
            overall_largest_ring_mcs_bond = ring_bond_pair
            overall_largest_ring = largest_ring

    ligand_bond_atoms = overall_largest_ring_mcs_bond[1].atoms
    ligand_bond_type = overall_largest_ring_mcs_bond[1].bond_type
    mcs.ligand.remove_bond(overall_largest_ring_mcs_bond[1])

    template_bond_atoms = overall_largest_ring_mcs_bond[0].atoms
    template_bond_type = overall_largest_ring_mcs_bond[0].bond_type
    mcs.mcs_scaffold.remove_bond(overall_largest_ring_mcs_bond[0])
    mcs.mcs_scaffold.add_hydrogens(mode="missing")
    rdkit_scaffold = Chem.MolFromMol2Block(mcs.mcs_scaffold.to_string(), removeHs=False)
    w = Chem.SDWriter(str(template_filename))
    w.write(rdkit_scaffold)
    w.close()

    tmp_ligand_file = "tmp_ligand_mol.sdf"
    tmp_ligand_file = docking_folder / tmp_ligand_file
    mcs.ligand.add_hydrogens(mode="missing")
    mcs.ligand.normalise_labels()
    rdkit_ligand = Chem.MolFromMol2Block(mcs.ligand.to_string("mol2"), removeHs=False)
    Chem.MolToMolFile(rdkit_ligand, str(tmp_ligand_file))

    sp.check_output(
        [
            sys.executable,
            Path(__file__).parent / "template_alignment.py",
            "-l",
            str(Path(tmp_ligand_file).resolve()),
            "-t",
            str(Path(template_filename).resolve()),
            "-o",
            str(Path(ligand_filename).resolve()),
        ],
        env=os.environ,
    )

    # Ensure stereo centers were preserved by omega:
    rdkit_ligand_stereo_dict = (
        template_docking_mcs._return_stereocenter_hybridization_dict(rdkit_ligand)[0]
    )

    ccdc_input_ligand = io.MoleculeReader(str(Path(ligand_filename).resolve()))[0]
    ccdc_input_ligand.kekulize()
    input_ligand = Chem.MolFromMol2Block(
        ccdc_input_ligand.to_string("mol2"), removeHs=False
    )
    input_ligand_stereo_dict = (
        template_docking_mcs._return_stereocenter_hybridization_dict(input_ligand)[0]
    )

    ccdc_input_ligand.standardise_aromatic_bonds()
    mcs.ligand.standardise_aromatic_bonds()

    mcs_atoms = return_mcs_atoms_of_identical_mols(mcs.ligand, ccdc_input_ligand)
    # mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    # mcs_searcher.settings.ignore_hydrogens = True
    # mcs_atoms = mcs_searcher.search(
    #     mcs.ligand, ccdc_input_ligand, search_step_limit=1000000
    # )[0]

    stereo_center_pairs = [
        a for a in mcs_atoms if a[0].label in rdkit_ligand_stereo_dict.keys()
    ]
    for stereo_center_pair in stereo_center_pairs:
        if (
            rdkit_ligand_stereo_dict[stereo_center_pair[0].label]
            != input_ligand_stereo_dict[stereo_center_pair[1].label]
        ):
            print("Omega returned wrong stereo center")
            return False

    # stitch macrocycle back together
    a0 = [a[1] for a in mcs_atoms if a[0] == ligand_bond_atoms[0]][0]
    a1 = [a[1] for a in mcs_atoms if a[0] == ligand_bond_atoms[1]][0]
    ccdc_input_ligand.add_bond(ligand_bond_type, a0, a1)
    ccdc_input_ligand.remove_hydrogens()
    ccdc_input_ligand.add_hydrogens()
    with io.MoleculeWriter(ligand_filename) as w:
        w.write(ccdc_input_ligand)
    return True


def _write_starting_ligand(
    ligand_mol, ligand_filename, docking_folder, template, template_filename
):
    rdkit_template = Chem.MolFromMol2Block(template.to_string(), removeHs=False)
    w = Chem.SDWriter(str(template_filename))
    w.write(rdkit_template)
    w.close()

    tmp_ligand_file = "tmp_ligand_mol.sdf"
    tmp_ligand_file = docking_folder / tmp_ligand_file
    ligand_mol.add_hydrogens(mode="missing")
    ligand_mol.normalise_labels()
    rdkit_ligand = Chem.MolFromMol2Block(ligand_mol.to_string("mol2"), removeHs=False)
    Chem.MolToMolFile(rdkit_ligand, str(tmp_ligand_file))

    # template_alignment._omega_overlay(
    #     str(Path(tmp_ligand_file).resolve()),
    #     str(Path(template_filename).resolve()),
    #     str(Path(ligand_filename).resolve()),
    # )
    sp.check_output(
        [
            sys.executable,
            Path(__file__).parent / "template_alignment.py",
            "-l",
            str(Path(tmp_ligand_file).resolve()),
            "-t",
            str(Path(template_filename).resolve()),
            "-o",
            str(Path(ligand_filename).resolve()),
        ],
        env=os.environ,
    )

    # Ensure stereo centers were preserved by omega:
    rdkit_ligand_stereo_dict = (
        template_docking_mcs._return_stereocenter_hybridization_dict(rdkit_ligand)[0]
    )

    ccdc_input_ligand = io.MoleculeReader(str(Path(ligand_filename).resolve()))[0]
    ccdc_input_ligand.kekulize()
    input_ligand = Chem.MolFromMol2Block(
        ccdc_input_ligand.to_string("mol2"), removeHs=False
    )
    input_ligand_stereo_dict = (
        template_docking_mcs._return_stereocenter_hybridization_dict(input_ligand)[0]
    )

    ccdc_input_ligand.standardise_aromatic_bonds()
    ligand_mol.standardise_aromatic_bonds()
    mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    mcs_searcher.settings.ignore_hydrogens = True
    mcs_atoms = mcs_searcher.search(
        ligand_mol, ccdc_input_ligand, search_step_limit=1000000
    )[0]

    stereo_center_pairs = [
        a for a in mcs_atoms if a[0].label in rdkit_ligand_stereo_dict.keys()
    ]
    for stereo_center_pair in stereo_center_pairs:
        if (
            rdkit_ligand_stereo_dict[stereo_center_pair[0].label]
            != input_ligand_stereo_dict[stereo_center_pair[1].label]
        ):
            print("Omega returned wrong stereo center")
            return False
    return True


def gold_decoy_docking(reference_docking_job, args=None):
    from ccdc_roche_scoring import join_docked_rf_counts

    mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
    mcs_searcher.settings.ignore_hydrogens = True

    df = pd.read_csv(Path(reference_docking_job) / Path("docked_rf_count_df.csv"))
    best_solns = join_docked_rf_counts.get_best_docking_solutions(df, args.target)[0]
    best_solns = [
        Path(ligand_file).parent.stem
        for ligand_file in best_solns["ligand_file"].values
    ]
    reference_docking_dirs = [
        reference_dir
        for reference_dir in Path(reference_docking_job).glob("*_*")
        if reference_dir.stem in best_solns
    ]
    for reference_docking_dir in reference_docking_dirs:
        if reference_docking_dir.is_dir():
            docking_folder = Path(reference_docking_dir.name)
            if not docking_folder.is_dir():
                docking_folder.mkdir()
            pdb_id = docking_folder.name.split("_")[0].lower()
            ligand_identifier = docking_folder.name.split("_")[1]

            reference_ligand_file = list(reference_docking_dir.glob("best_soln_*.sdf"))[
                0
            ].resolve()

            # load scaffold file
            scaffold_file = list(reference_docking_dir.glob("scaffold_*.sdf"))[
                0
            ].resolve()
            with io.MoleculeReader(str(scaffold_file)) as rdr:
                scaffold = rdr[0]

            # load cavity file
            cavity_ligand_filepath = list(reference_docking_dir.glob("cavity_*.mol2"))[
                0
            ].resolve()
            with io.MoleculeReader(str(cavity_ligand_filepath)) as rdr:
                cavity_ligand = rdr[0]

            # load target filesrn
            dry_receptor_file = (
                reference_docking_dir / Path("gold_protein.mol2")
            ).resolve()
            with io.EntryReader(str(dry_receptor_file)) as rdr:
                target_protein = protein.Protein.from_entry(rdr[0])
                # target_protein = rdr[0]

            # run docking
            docker = Docker()
            ligand_filepath = list(reference_docking_dir.glob("docking_input_*.mol2"))[
                0
            ].resolve()
            _dock(
                docker,
                dry_receptor_file,
                target_protein,
                cavity_ligand,
                ligand_filepath,
                docking_folder,
                pdb_id,
                srn,
                water_paths=False,
                scaffold=False,
                reference_ligand_file=reference_ligand_file,
                autoscale=25,
                args=args,
            )

            _select_best_decoy(
                docking_folder,
                mcs_searcher,
                scaffold,
                srn,
                pdb_id,
                reference_ligand_file,
            )


def _setup_docking(
    docking_folder,
    targets,
    pdb_id,
    project_home,
    ligand_mol,
    srn,
    scaffold,
    template_filename,
    native_ligand,
    strict_scaffold,
    mcs,
    args,
):
    if not docking_folder.is_dir():
        docking_folder.mkdir()

    # get target structure and water
    dry_receptor_file = (
        Path(project_home) / Path("apo_proteins") / Path(f"{pdb_id}_dry.mol2")
    )
    if not Path(str(dry_receptor_file)).is_file():
        pdb = [p for p in targets if pdb_id.lower() in str(p).lower()][0]
        print("Preparing protein MOL2 files...")
        _protein_preparation(pdb, dry_receptor_file, pdb_id, target=args.target)
    with io.EntryReader(str(dry_receptor_file)) as t_rdr:
        target_protein = protein.Protein.from_entry(t_rdr[0])
        water_paths = [
            str(water_path.resolve())
            for water_path in dry_receptor_file.parent.glob(f"{pdb_id}_water*.mol2")
        ]
        if len(water_paths) == 0:
            print("No water molecules selected for docking...")
            water_paths = False
    target_protein.identifier = pdb_id

    ligand_filename = Path(docking_folder) / Path(
        f"docking_input_{srn}.mol2"
    )  # mol2 to preserve atom types

    if mcs.ligand.rings:
        largest_ring = max(mcs.ligand.rings, key=lambda ring: len(ring.atoms))

    ready_to_dock = False
    if (
        len(mcs.ligand.heavy_atoms) > 80
        and mcs.ligand.rings
        and len(largest_ring.atoms) > 15
    ):
        if _write_starting_ligand_macrocycle(
            ligand_filename,
            docking_folder,
            template_filename,
            mcs,
        ):
            ready_to_dock = True
    else:
        if _write_starting_ligand(
            ligand_mol, ligand_filename, docking_folder, scaffold, template_filename
        ):
            ready_to_dock = True

    if ready_to_dock:
        docker = Docker()
        _dock(
            docker,
            dry_receptor_file,
            target_protein,
            native_ligand,
            ligand_filename,
            docking_folder,
            pdb_id,
            srn,
            water_paths,
            scaffold=scaffold,
            strict_scaffold=strict_scaffold,
            **vars(args),
        )


def _write_best_soln_pocket(gold_conf, ligand_file):
    docking_settings = Docker.Settings().from_file(gold_conf)
    docking_results = Docker.Results(docking_settings)
    csd_ligand_entry = docking_results.DockedLigandReader(
        ligand_file, docking_settings
    )[0]

    ligand_identifier = csd_ligand_entry.identifier
    # setup protein-ligand-complex
    protein = docking_results.make_complex(csd_ligand_entry)
    apo_protein = protein.copy()
    apo_protein.remove_ligand(":" + ligand_identifier)
    apo_protein.remove_unknown_atoms()
    apo_protein.identifier = ligand_identifier

    apo_pocket = protein.copy()

    # workaround: gold somehow messes up ligand sometimes
    apo_pocket.remove_ligand(":" + ligand_identifier)
    apo_pocket.remove_unknown_atoms()

    csd_ligand_entry.identifier = "LIG000"
    apo_pocket.add_ligand(csd_ligand_entry.molecule)

    apo_pocket = _cut_out_binding_site_by_distance(
        apo_pocket, csd_ligand_entry.molecule, cutoff=5.0
    )
    apo_pocket.remove_unknown_atoms()
    apo_pocket.identifier = ligand_identifier

    pocket_file = Path(ligand_file).parent / "best_soln_pocket.mol2"
    protein_soln_file = Path(ligand_file).parent / "best_soln_protein.pdb"
    with io.EntryWriter(pocket_file) as w:
        w.write(apo_pocket)

    apo_protein.kekulize()
    with io.EntryWriter(protein_soln_file) as w:
        w.write(apo_protein)
    return


def target_setup(target):
    targets = list(
        Path(f"{target}_bindingsites/").glob(
            "*prot3d.pdb"
        )
    )
    template_ligands_sdf = list(
        Path(f"{target}_bindingsites/").glob(
            "*prot3d.sdf"
        )
    )
    pdb_ligand_df_file = Path(
        f"{target}_bindingsites/Structure.info"
    )
    if pdb_ligand_df_file.is_file():
        pdb_ligand_df = pd.read_csv(pdb_ligand_df_file, sep="\t")
        pdb_ligand_df["LIGCLASHES"] = pdb_ligand_df["LIGCLASHES"].str.replace("-", "0")
        pdb_ligand_df = pdb_ligand_df[pdb_ligand_df["LIGSYMCONTACT"].isna() == False]
        pdb_ligand_df = pdb_ligand_df.astype(
            {"RESOLUTION": float, "LIGSYMCONTACT": int, "LIGCLASHES": int}
        )
        pdb_ligand_df = pdb_ligand_df[
            (pdb_ligand_df["RESOLUTION"] <= 2.5)
            & (pdb_ligand_df["LIGSYMCONTACT"] == 0)
            & (pdb_ligand_df["LIGCLASHES"] == 0)
        ]
        if "EXCLUDE_FROM_DOCKING" in pdb_ligand_df.columns:
            pdb_ligand_df = pdb_ligand_df[
                pdb_ligand_df["EXCLUDE_FROM_DOCKING"] == False
            ]
        high_q_strucids = pdb_ligand_df["STRUCID"].to_list()
        targets = [
            target_path
            for target_path in targets
            if target_path.name.split("_")[0] in high_q_strucids
        ]
        template_ligands_sdf = [
            lig
            for lig in template_ligands_sdf
            if lig.name.split("_")[0] in high_q_strucids
        ]
    template_ligands_sdf = [str(lig) for lig in template_ligands_sdf]

    return template_ligands_sdf, targets


def template_ligand_setup(args, pdb_ligand_df_file, template_ligands_sdf):
    native_ligand_entries = []
    with io.EntryReader(template_ligands_sdf) as sdf_rdr:
        for cnt, native_ligand_entry in enumerate(sdf_rdr):
            if "strucid" in native_ligand_entry.attributes:
                strucid = native_ligand_entry.attributes["strucid"]
            elif args.target != "default":
                strucid = native_ligand_entry.identifier.split(".")[0]
            elif pdb_ligand_df_file.is_file() and args.target == "default":
                pdb_ligand_df = pd.read_csv(pdb_ligand_df_file)
                ligand_file = str(
                    Path(native_ligand_entry.attributes["$File"]).resolve()
                )
                strucid = pdb_ligand_df[pdb_ligand_df["template_file"] == ligand_file][
                    "strucid"
                ]
                if strucid.shape[0] == 1:
                    strucid = strucid.to_list()[0]
                else:
                    strucid = False
            else:
                strucid = native_ligand_entry.identifier.split(".")[0]

            native_ligand = native_ligand_entry.molecule
            if args.target != "default":
                max_comp_size = 0
                max_comp = None
                for native_ligand_comp in native_ligand.components:
                    comp_size = len(native_ligand_comp.heavy_atoms)
                    if comp_size > max_comp_size:
                        max_comp_size = comp_size
                        max_comp = native_ligand_comp
                native_ligand = max_comp
                native_ligand.add_hydrogens(mode="missing")

            native_ligand.normalise_labels()
            new_native_ligand_entry = entry.Entry.from_molecule(native_ligand)
            new_native_ligand_entry.attributes = native_ligand_entry.attributes
            native_ligand_entry = new_native_ligand_entry
            if strucid:
                native_ligand_entry.attributes["STRUCID"] = strucid
                native_ligand_entries.append(native_ligand_entry)
    return native_ligand_entries


def gold_scaffold_docking(ligands_for_alignment, project_home, args=None):
    rel_mcs_size_to_ligand_threshold = 0.25
    rel_mcs_size_to_native_ligand_threshold = 0.30
    similarity_threshold = 0.7
    pdb_ligand_df_file = Path(project_home, "pdb_ligand.csv")
    print("Target: ", args.target)
    if args.target == "default":
        targets = list(
            ((Path(project_home) / Path("tmp_aligned_for_MOE_sanitized")).glob("*.pdb"))
        )
        template_ligands_sdf = str(
            Path(project_home) / "tmp_aligned_3d_sdf_sanitized/native_ligands.sdf"
        )
        if pdb_ligand_df_file.is_file():
            pdb_ligand_df = pd.read_csv(pdb_ligand_df_file)
            pdb_ligand_df = pdb_ligand_df[pdb_ligand_df["is_high_quality"] == True]
            pdb_ligand_df = pdb_ligand_df[
                pdb_ligand_df["template_file"].isna() == False
            ]
    else:
        template_ligands_sdf, targets = target_setup(args.target)

    series_df = Path("../../series_assignment.csv")
    if series_df.is_file():
        series_df = pd.read_csv(series_df)
        series_df = series_df[
            [c for c in series_df.columns if c in ["Proasis ID", "SRN", "series"]]
        ].rename({"Proasis ID": "strucid"}, axis=1)
        pdb_ligand_df = pdb_ligand_df.join(series_df.set_index("strucid"), on="strucid")
    else:
        series_df = pd.DataFrame()

    native_ligand_entries = template_ligand_setup(
        args, pdb_ligand_df_file, template_ligands_sdf
    )

    ligands_for_alignment_sdf = Path(project_home) / Path(ligands_for_alignment)
    with io.EntryReader(str(ligands_for_alignment_sdf)) as lig_sd_rdr:
        for lig_cnt, ligand in enumerate(lig_sd_rdr):
            # try:
            ligand_entry = lig_sd_rdr[lig_cnt]
            if "SRN" in ligand_entry.attributes:
                srn = ligand_entry.attributes["SRN"].replace(" ", "")
            else:
                srn = ligand_entry.identifier.replace(" ", "")
                ligand_entry.attributes["SRN"] = srn

            ###
            # if 'RO6089776' not in srn:
            #     continue
            # if lig_cnt == 0:
            #     continue
            ###

            print(srn)
            series_template_strucids = None
            if series_df.shape[0] > 1:
                series = series_df[series_df["SRN"] == srn]["series"].values[0]
                if series == series:
                    series_template_strucids = pdb_ligand_df[
                        pdb_ligand_df["series"] == series
                    ]["strucid"].values
                else:
                    continue

            ligand_mol = ligand_entry.molecule.components[0]

            # kekulize with RDKit because it handles aromaticity better
            rdkit_ligand = Chem.MolFromMolBlock(
                ligand_mol.to_string("sdf"), removeHs=False
            )
            ligand_mol = ligand_mol.from_string(Chem.MolToMolBlock(rdkit_ligand))

            ligand_mol.normalise_labels()
            templates_df = _mcs_templates_df(
                native_ligand_entries, ligand_mol, series_template_strucids
            )

            # Remove template structures that were solved after the assay date
            if args.target == "pde-10":
                proasis_df = pd.read_csv("../../proasis_date.csv").astype(
                    {"DATESOLVED": "datetime64"}
                )
                templates_df = templates_df.join(
                    proasis_df.set_index("STRUCID"), on="template_strucid"
                )
                assay_date = ligand_entry.attributes[
                    "PDE10_FULL_LEN_CGMP_SPA_IC50_h-PDE10A(14-779)-E.Coli-c: AP005128;Min;EXP Test Date"
                ]
                assay_date = pd.to_datetime(assay_date)
                templates_df = templates_df[templates_df["DATESOLVED"] < assay_date]

            dock = False
            for index, row in templates_df.iterrows():
                criterion1 = (
                    row["rel_mcs_size_to_ligand"] >= rel_mcs_size_to_ligand_threshold
                    and row["rel_mcs_size_to_native_ligand"]
                    >= rel_mcs_size_to_native_ligand_threshold
                )
                criterion2 = (
                    row["abs_mcs_size"] >= 10
                    and row["similarity"] >= similarity_threshold
                )
                if index < 1 and row["abs_mcs_size"] > 5 and (criterion1 or criterion2):
                    pdb_id = row["template_strucid"].split(".")[0].lower()
                    native_ligand = row["native_ligand"].molecule
                    docking_folder = Path(f"{pdb_id}_{srn}_{lig_cnt}")
                    scaffold = row["scaffold"]
                    mcs = row["mcs_object"]
                    template_filename = docking_folder / Path(f"scaffold_{srn}.sdf")
                    strict_scaffold = mcs.return_mcs_scaffold(
                        partial_ring_matches_allowed=False, ignore_hydrogens=False
                    )[0]

                    _setup_docking(
                        docking_folder,
                        targets,
                        pdb_id,
                        project_home,
                        ligand_mol,
                        srn,
                        scaffold,
                        template_filename,
                        native_ligand,
                        strict_scaffold,
                        mcs,
                        args,
                    )
                    _select_best_soln(
                        docking_folder,
                        scaffold,
                        row,
                        srn,
                        pdb_id,
                        ligand_entry.attributes,
                        template_filename,
                    )
                    dock = True
            if not dock:
                print("No suitable scaffold for ", srn)

            # except Exception as e:
            #     print('EXCEPTION')
            #     print(e)
            #     continue


def gold_rescoring(ligand_file: str, protein_file: str):
    print("rescoring......")

    with io.MoleculeReader(ligand_file) as rdr:
        ligand = rdr[0]
    ligand.add_hydrogens(mode="missing")
    target_protein = protein.Protein.from_file(protein_file)

    docker = Docker()
    settings = docker.settings
    settings.fitness_function = None
    settings.rescore_function = "plp"
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
    settings._settings.set_rescore_with_simplex(False)
    settings._settings.set_fix_protein_rotatable_bonds(True)
    settings.add_protein_file(str(protein_file))

    settings.binding_site = settings.BindingSiteFromLigand(target_protein, ligand, 10.0)
    settings.add_ligand_file(ligand_file)
    settings.output_file = str(Path(f"./rescored_ligand.sdf"))

    gold_conf = "gold_rescoring.conf"
    settings.make_absolute_file_names(gold_conf)
    settings.write(gold_conf)

    docker = Docker()
    docker.settings = docker.settings.from_file(gold_conf)
    docker.dock(file_name=gold_conf)

    ligand_attributes = {}
    if Path("rescored_ligand.sdf").is_file():
        ligand_entry = io.EntryReader("rescored_ligand.sdf")[0]
        ligand_attributes = ligand_entry.attributes
        ligand_attributes["Gold.PLP.Chemscore.Protein.Energy"] = 0
        ligand_entry.attributes = ligand_attributes
        with io.EntryWriter("rescored_ligand.sdf") as w:
            w.write(ligand_entry)

    return ligand_attributes


def main():
    args = parse_args()
    args.flexible_residues = [s.upper() for s in args.flexible_residues if s.strip()]
    project_home = "../.."

    print("Project home:", Path(project_home).resolve())
    if args.input_ligands:
        gold_scaffold_docking(args.input_ligands, project_home=project_home, args=args)
        print("finished scaffold docking")
    if args.docking_dir:
        gold_decoy_docking(args.docking_dir, args=args)
        print("finished decoy docking")


if __name__ == "__main__":
    main()
