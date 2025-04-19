#!/usr/bin/env python

"""
Post-processing of docking poses.
"""


import argparse
from functools import partial
from pathlib import Path
from collections import defaultdict

import pandas as pd
from ccdc import descriptors, io, search
from pathos.multiprocessing import ProcessingPool

from ccdc_roche_scoring import q_values


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="""
        Post-process docking poses. Count donor-donor, acceptor-acceptor contacts and steric clashes.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument(
        "-np",
        "--nproc",
        type=int,
        help="Number of parallel processes for multiprocessing.",
        default=24,
    )

    parser.add_argument(
        "-isdf", help="SD file with ligand poses", default="ligands.sdf"
    )

    parser.add_argument("--outdir", "-o", help="Output directory", default=".")

    return parser.parse_args()


def smarts_to_hits(
    smarts,
    db,
    atom_property=False,
    max_hits_per_structure=None,
    torsional_constraint=None,
):
    searcher = search.SubstructureSearch()
    searcher.Settings(max_hits_per_structure=max_hits_per_structure)
    smarts = search.SMARTSSubstructure(smarts)
    if atom_property:
        smarts.atoms[0].label_match = "^_Z"
    searcher.add_substructure(smarts)

    if torsional_constraint:
        searcher.add_torsion_angle_constraint(*torsional_constraint)
    return searcher.search(db, max_hits_per_structure=max_hits_per_structure)


def bond_is_rotatable(bond):
    if not bond.is_rotatable:
        return False
    if bond.sybyl_type == "am":
        return False

    aniline_smarts = search.SMARTSSubstructure("[N$([NH2]c)]")
    sulfonyl_smarts = search.SMARTSSubstructure("[S$(S(~O)(~O)[*r])]")
    for smarts in [aniline_smarts, sulfonyl_smarts]:
        for atom in bond.atoms:
            if smarts.match_atom(atom):
                return False

    return True


def not_lone_pair_clash(contact):
    nitrile_smarts = search.SMARTSSubstructure("[N$(N#C)]")
    contact_atoms = contact.atoms
    for atom in contact_atoms:
        if nitrile_smarts.match_atom(atom):
            other_acceptor = [a for a in contact_atoms if a != atom][0]
            contact_angle = descriptors.MolecularDescriptors.atom_angle(
                atom.neighbours[0], atom, other_acceptor
            )
            if -170 < contact_angle or 170 > contact_angle:
                return True
    return False


class IntramolecularContactChecker(object):
    def __init__(self, contact_atom_0, contact_atom_1, vdw_distance, ccdc_mol, contact):
        self.contact_atom_0 = contact_atom_0
        self.contact_atom_1 = contact_atom_1
        self.vdw_distance = vdw_distance
        self.ccdc_mol = ccdc_mol
        self.contact = contact
        self.substructure_dict = self.get_substructure_dict()
        return

    def get_substructure_dict(self):
        smarts_dict = {
            "aromatic_oxygen": "[O$(O1-*=*-*=*1)]",
            "aromatic_sulfur": "[S$([S]1@-[*]@=[*]@-[*]@=[*]1)]",
            "amidinium": "[#7$([#7]=,:[#6]-!@[#7]),$([#7]-,:[#6]=!@[#7]),$([#7]=!@[#6]-,:[#7]),$([#7]-!@[#6]=,:[#7])]",
            "aniline": "[N$(N(-[CX4])-[#6r]),$(N(-[#6])-@[#6])]",
            "urea": "[N$(N-C(=[O,N])-[N])]",
        }
        substructure_dict = {}
        for sub in smarts_dict:
            substructure_dict[sub] = search.SMARTSSubstructure(smarts_dict[sub])
        return substructure_dict

    def contact_is_sigma_hole(self):
        if self.substructure_dict["aromatic_sulfur"].match_atom(self.contact_atom_0):
            plane = descriptors.MolecularDescriptors.ring_plane(
                self.contact_atom_0.rings[0]
            )
            if plane.point_distance(self.contact_atom_1.coordinates) < 1:
                return True
        elif self.substructure_dict["aromatic_sulfur"].match_atom(self.contact_atom_1):
            plane = descriptors.MolecularDescriptors.ring_plane(
                self.contact_atom_1.rings[0]
            )
            if plane.point_distance(self.contact_atom_0.coordinates) < 1:
                return True
        return False

    def contact_is_with_weakest_acceptor_of_neighbours(self, neighbours):
        acceptor_strengths = {"N": 4, "O": 3, "S": 2, "Cl": 1, "Br": 1, "I": 1}
        weakest_acceptor_neigbour = {"weakest_acceptor_neigbour": None, "strength": 0}
        for neighbour in neighbours:
            if neighbour != self.contact_atom_0 and neighbour != self.contact_atom_1:
                if neighbour.is_acceptor:
                    strength = acceptor_strengths[neighbour.atomic_symbol]
                else:
                    strength = 0
                if weakest_acceptor_neigbour["weakest_acceptor_neigbour"] is None:
                    weakest_acceptor_neigbour["weakest_acceptor_neigbour"] = neighbour
                    weakest_acceptor_neigbour["strength"] = strength
            else:
                contact_atom_strength = acceptor_strengths[neighbour.atomic_symbol]

        if contact_atom_strength <= weakest_acceptor_neigbour["strength"]:
            return True

        return False

    def is_acc_acc_clash(self):
        """
        Returns True if two acceptors are clashing.
        """

        if (
            self.contact_atom_0.is_acceptor
            and self.contact_atom_1.is_acceptor
            and self.vdw_distance < 0.2
        ):
            if self.substructure_dict["aromatic_oxygen"].match_atom(
                self.contact_atom_0
            ) or self.substructure_dict["aromatic_oxygen"].match_atom(
                self.contact_atom_1
            ):
                return False

            if self.contact_is_sigma_hole():
                return False

            shortest_path_bonds = self.ccdc_mol.shortest_path_bonds(
                self.contact_atom_0, self.contact_atom_1
            )
            rotatable_bonds_on_path = [
                bond for bond in shortest_path_bonds if bond_is_rotatable(bond)
            ]

            num_rotatable_bonds_on_path = len(rotatable_bonds_on_path)

            if num_rotatable_bonds_on_path >= 5:
                return True

            elif num_rotatable_bonds_on_path == 0:
                return False

            if not_lone_pair_clash(self.contact):
                return False

            else:
                # dimethoxy use case
                if (
                    shortest_path_bonds[0].is_rotatable
                    or shortest_path_bonds[-1].is_rotatable
                ):
                    return False

                rotatable_bonds_on_path_atoms = []
                if rotatable_bonds_on_path:
                    rotatable_bonds_on_path_atoms = rotatable_bonds_on_path[0].atoms
                for rot_atom in rotatable_bonds_on_path_atoms:
                    rot_atom_neighbours = [
                        a
                        for a in rot_atom.neighbours
                        if a not in rotatable_bonds_on_path_atoms
                    ]
                    if (
                        self.contact_atom_0 in rot_atom_neighbours
                        or self.contact_atom_1 in rot_atom_neighbours
                    ):
                        if self.contact_is_with_weakest_acceptor_of_neighbours(
                            rot_atom_neighbours
                        ):
                            return False
                    else:
                        continue
                return True
        return False

    def is_don_don_clash(self):
        if (
            self.contact_atom_0.atomic_symbol == "H"
            and self.contact_atom_1.atomic_symbol == "H"
            and self.vdw_distance < 0.2
            and self.contact.is_in_line_of_sight
        ) or (
            self.contact_atom_0.atomic_symbol == "H"
            and self.contact_atom_1.atomic_symbol == "H"
            and self.vdw_distance < -0.1
        ):
            donor_atom_0 = self.contact_atom_0.neighbours[0]
            donor_atom_1 = self.contact_atom_1.neighbours[0]
            if donor_atom_0.is_donor and donor_atom_1.is_donor:
                bonds_on_path = self.ccdc_mol.shortest_path_bonds(
                    donor_atom_0, donor_atom_1
                )
                rotatable_bonds_on_path = [
                    bond for bond in bonds_on_path if bond_is_rotatable(bond)
                ]

                if len(rotatable_bonds_on_path) > 0:
                    if self.substructure_dict["aniline"].match_atom(
                        donor_atom_0
                    ) and self.substructure_dict["aniline"].match_atom(donor_atom_1):
                        return False

                    if self.substructure_dict["urea"].match_atom(
                        donor_atom_0
                    ) and self.substructure_dict["urea"].match_atom(donor_atom_1):
                        return False

                    # check amidinium
                    if self.substructure_dict["amidinium"].match_atom(
                        donor_atom_0
                    ) and self.substructure_dict["amidinium"].match_atom(donor_atom_1):
                        return False

                    return True
        return False


def count_cis_amides(mol):
    cis_amide_smarts = "[CX4]-!@N(H)-!@C(=O)-C"
    torsion_constrain_args = (
        "cis_amide",
        (0, 2),
        (0, 1),
        (0, 3),
        (0, 4),
        (-20.0, 20.0),
    )
    hits = smarts_to_hits(
        cis_amide_smarts, mol, torsional_constraint=torsion_constrain_args
    )

    return len(hits)


def get_primary_contacts(contacts):
    primary_contacts_dict = {}
    for contact in contacts:
        for atom in contact.atoms:
            if atom not in primary_contacts_dict.keys():
                primary_contacts_dict[atom] = contact
            else:
                if contact.length < primary_contacts_dict[atom].length:
                    primary_contacts_dict[atom] = contact
    return primary_contacts_dict


def count_bad_contacts(ccdc_mol):
    """Return dictionary with intramolecular clashes for ccdc molecule.

    Args:
        ccdc_mol: ccdc_mol

    Returns:
        Dictionary with bad contact counts.
    """

    ccdc_mol.remove_unknown_atoms()
    ccdc_mol.standardise_aromatic_bonds()
    bad_contacts_dict = {
        "acc_acc_contact": 0,
        "don_don_contact": 0,
        "intramolecular_steric_clashes": 0,
        "cis_amide_count": 0,
    }

    # contacts = ccdc_mol.contacts(path_length_range=(3, 999))
    # primary_contacts = get_primary_contacts(contacts)
    # hbond_atoms = [(hb.atoms[0], hb.atoms[2]) for hb in ccdc_mol.hbonds()]
    # for contact in contacts:
    #     contact_atom_0 = contact.atoms[0]
    #     contact_atom_1 = contact.atoms[1]

    #     if (contact_atom_0, contact_atom_1) in hbond_atoms or (
    #         contact_atom_1,
    #         contact_atom_0,
    #     ) in hbond_atoms:
    #         continue

    #     vdw_distance = (
    #         contact.length - contact_atom_0.vdw_radius - contact_atom_1.vdw_radius
    #     )

    #     if (
    #         vdw_distance < -0.5
    #         and contact_atom_0.atomic_symbol != "H"
    #         and contact_atom_1.atomic_symbol != "H"
    #     ):
    #         shortest_path_bonds = ccdc_mol.shortest_path_bonds(
    #             contact_atom_0, contact_atom_1
    #         )
    #         rotatable_bonds_on_path = len(
    #             [bond for bond in shortest_path_bonds if bond.is_rotatable]
    #         )
    #         if rotatable_bonds_on_path > 0:
    #             threshold = -0.6
    #             if [b for b in shortest_path_bonds if b.sybyl_type != "1"]:
    #                 threshold = -0.62
    #             if len(shortest_path_bonds) == 3 and vdw_distance < threshold:
    #                 bad_contacts_dict["intramolecular_steric_clashes"] += 1
    #             elif len(shortest_path_bonds) > 3:
    #                 bad_contacts_dict["intramolecular_steric_clashes"] += 1
    #         continue

    #     contact_checker = IntramolecularContactChecker(
    #         contact_atom_0, contact_atom_1, vdw_distance, ccdc_mol, contact
    #     )
    #     if contact_checker.is_acc_acc_clash():
    #         bad_contacts_dict["acc_acc_contact"] += 1

    #     if contact_checker.is_don_don_clash():
    #         bad_contacts_dict["don_don_contact"] += 1

    # bad_contacts_dict["cis_amide_count"] = count_cis_amides(ccdc_mol)
    bad_contacts_dict["q2_gmean"] = q_values.QValue().get_q2_gmean(ccdc_mol)

    return bad_contacts_dict


def return_clash_dict(entry_index, sdf):
    # try:
    mol = io.MoleculeReader(sdf)[entry_index]
    clash_dict = count_bad_contacts(mol)
    clash_dict["mol_id"] = mol.identifier
    return clash_dict

    # except:
    #     return None


def single_process(args):
    rdr = io.EntryReader(args.isdf)
    combined_clash_dict = defaultdict(list)
    with io.EntryWriter(str(Path(args.outdir) / "clash_check.sdf")) as w:
        for i in range(len(rdr)):
            clash_dict = return_clash_dict(i, args.isdf)
            if clash_dict:
                e = rdr.entry(clash_dict["mol_id"])
                e.attributes.update(clash_dict)
                combined_clash_dict["mol_id"].append(clash_dict["mol_id"])
                combined_clash_dict["acc_acc_contact"].append(
                    clash_dict["acc_acc_contact"]
                )
                combined_clash_dict["don_don_contact"].append(
                    clash_dict["don_don_contact"]
                )
                combined_clash_dict["intramolecular_steric_clashes"].append(
                    clash_dict["intramolecular_steric_clashes"]
                )

                combined_clash_dict["cis_amide_count"].append(
                    clash_dict["cis_amide_count"]
                )
                combined_clash_dict["q2_gmean"].append(clash_dict["q2_gmean"])
                w.write(e)
    return combined_clash_dict


def multi_process(args):
    rdr = io.EntryReader(args.isdf)
    parallel_func = partial(return_clash_dict, sdf=args.isdf)
    pool = ProcessingPool(args.nproc)
    it = range(len(rdr))
    output = pool.map(parallel_func, it)
    combined_clash_dict = defaultdict(list)
    with io.EntryWriter(str(Path(args.outdir) / "clash_check.sdf")) as w:
        for clash_dict in output:
            if clash_dict:
                e = rdr.entry(clash_dict["mol_id"])
                e.attributes.update(clash_dict)
                combined_clash_dict["mol_id"].append(clash_dict["mol_id"])
                combined_clash_dict["acc_acc_contact"].append(
                    clash_dict["acc_acc_contact"]
                )
                combined_clash_dict["don_don_contact"].append(
                    clash_dict["don_don_contact"]
                )
                combined_clash_dict["intramolecular_steric_clashes"].append(
                    clash_dict["intramolecular_steric_clashes"]
                )

                combined_clash_dict["cis_amide_count"].append(
                    clash_dict["cis_amide_count"]
                )
                combined_clash_dict["q2_gmean"].append(clash_dict["q2_gmean"])
                w.write(e)

    return combined_clash_dict


def main():
    args = parse_args()

    if args.nproc == 1:
        combined_clash_dict = single_process(args)
    else:
        combined_clash_dict = multi_process(args)
    df = pd.DataFrame(combined_clash_dict)
    df.to_csv(str(Path(args.outdir) / "clash_check.csv"), index=False, sep="\t")
    return


if __name__ == "__main__":
    main()
