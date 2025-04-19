import warnings

warnings.filterwarnings("ignore", ".*deprecated.*", DeprecationWarning, ".*", 0)
import itertools
import time

from ccdc import descriptors, search
from rdkit import Chem


class MCS(object):
    """Calculate Maximum Common Substructure between ligand and template."""

    def __init__(self, ligand, template):
        """

        :param ligand:
        :param template:
        """
        self.ligand = ligand
        self.template = template
        (
            self.default_mcs_atoms,
            self.default_mcs_bonds,
        ) = self.return_default_mcs_scaffold()
        self.mcs_scaffold, self.mcs_atoms, self.mcs_bonds = self.return_mcs_scaffold()

    def return_mcs_searcher(self, ignore_hydrogens=True):
        mcs_searcher = descriptors.MolecularDescriptors.MaximumCommonSubstructure()
        # mcs_searcher.settings.ignore_hydrogens = ignore_hydrogens
        mcs_searcher.settings.ignore_hydrogens = True
        mcs_searcher.settings.check_hydrogen_count = not ignore_hydrogens
        mcs_searcher.settings.check_bond_type = False
        mcs_searcher.settings.check_bond_count = True
        return mcs_searcher

    def return_default_mcs_scaffold(self, ignore_hydrogens=True):
        mcs_searcher = self.return_mcs_searcher(ignore_hydrogens)
        print("Get MCS...")
        t1 = time.time()

        if self.ligand.largest_ring_size > 20:
            mcs_atoms, mcs_bonds = macrocycle_mcs(
                self.template, self.ligand, mcs_searcher
            )
        else:
            mcs_atoms, mcs_bonds = mcs_searcher.search(
                self.template, self.ligand, search_step_limit=1000000
            )
        print("MCS done: ", time.time() - t1)
        return mcs_atoms, mcs_bonds

    def return_mcs_scaffold(
        self, partial_ring_matches_allowed=True, ignore_hydrogens=True
    ):
        mcs_searcher = self.return_mcs_searcher(ignore_hydrogens)

        # assign mcs atoms to copy of template
        scaffold = self.template.copy()
        if self.ligand.largest_ring_size > 20:
            mcs_atoms, mcs_bonds = macrocycle_mcs(scaffold, self.ligand, mcs_searcher)
        else:
            mcs_atoms, mcs_bonds = mcs_searcher.search(
                scaffold, self.ligand, search_step_limit=1000000
            )

        if partial_ring_matches_allowed:
            mcs_bonds = [
                mcs_pair
                for mcs_pair in mcs_bonds
                if mcs_pair[0].is_cyclic == mcs_pair[1].is_cyclic
            ]
            mcs_template_atoms_from_bonds = list(
                itertools.chain.from_iterable([b[0].atoms for b in mcs_bonds])
            )
            mcs_template_atoms_from_bonds = [
                a.label for a in mcs_template_atoms_from_bonds
            ]
            mcs_atoms = [
                mcs_pair
                for mcs_pair in mcs_atoms
                if mcs_pair[0].label in mcs_template_atoms_from_bonds
            ]
        else:
            mcs_atoms = [
                mcs_pair
                for mcs_pair in mcs_atoms
                if mcs_pair[0].is_cyclic == mcs_pair[1].is_cyclic
            ]

        remove_atoms_template = []

        for ta, la in mcs_atoms:
            if ta.sybyl_type in ["C.ar", "C.2", "N.3", "N.pl"]:
                ta_h = [a for a in ta.neighbours if a.atomic_symbol == "H"]
                la_h = [a for a in la.neighbours if a.atomic_symbol == "H"]
                if len(ta_h) > len(la_h):
                    remove_atoms_template.extend(ta_h)

        mcs_template_atoms = [a[0] for a in mcs_atoms]

        template_atoms = scaffold.heavy_atoms
        remove_atoms_template.extend(
            [a for a in template_atoms if a not in mcs_template_atoms]
        )

        # if ignore_hydrogens:
        #     remove_atoms_template = [
        #         a
        #         for a in template_atoms
        #         if a not in mcs_template_atoms and a.atomic_symbol != "H"
        #     ]
        # else:
        # remove_atoms_template = [
        #     a for a in template_atoms if a not in mcs_template_atoms
        # ]

        if partial_ring_matches_allowed:
            remove_atoms_template = (
                remove_atoms_template
                + self._mcs_partial_ring_match_from_bond(mcs_bonds)
            )
        else:
            remove_atoms_template = (
                remove_atoms_template + self._mcs_partial_ring_match(mcs_atoms)
            )
            remove_atoms_template = (
                remove_atoms_template
                + SubstitutentComparer().compare(self.ligand, mcs_atoms)
            )

        remove_atoms_template = remove_atoms_template + _compare_mcs_stereo_chemistry(
            scaffold, self.ligand, mcs_atoms
        )
        remove_atoms_template = list(set(remove_atoms_template))

        # remove partial ring matches
        while remove_atoms_template:
            mcs_searcher.settings.check_bond_count = False
            scaffold_atoms = scaffold.atoms
            remove_atoms_template = [
                a for a in remove_atoms_template if a in scaffold_atoms
            ]
            scaffold.remove_atoms(remove_atoms_template)
            scaffold = scaffold.heaviest_component

            if (
                self.ligand.largest_ring_size > 20
                and scaffold.rings
                and scaffold.largest_ring_size > 20
            ):
                mcs_atoms, mcs_bonds = macrocycle_mcs(
                    scaffold, self.ligand, mcs_searcher
                )
            else:
                mcs_atoms, mcs_bonds = mcs_searcher.search(
                    scaffold, self.ligand, search_step_limit=1000000
                )

            # mcs_atoms_ = mcs_searcher.search(scaffold, self.ligand, search_step_limit=1000000)[0]
            # remove_atoms_template = [mcs_pair[0] for mcs_pair in mcs_atoms if mcs_pair not in mcs_atoms_]
            remove_atoms_template = []

            if partial_ring_matches_allowed:
                remove_atoms_template = (
                    remove_atoms_template
                    + self._mcs_partial_ring_match_from_bond(mcs_bonds)
                )
            else:
                remove_atoms_template = (
                    self._mcs_partial_ring_match(mcs_atoms) + remove_atoms_template
                )
                remove_atoms_template = (
                    remove_atoms_template
                    + SubstitutentComparer().compare(self.ligand, mcs_atoms)
                )
            remove_atoms_template = (
                remove_atoms_template
                + _compare_mcs_stereo_chemistry(scaffold, self.ligand, mcs_atoms)
            )
            remove_atoms_template = list(set(remove_atoms_template))
        return scaffold, mcs_atoms, mcs_bonds

    def _mcs_partial_ring_match(self, mcs_atoms):
        """
        :param mcs_atoms:
        :return: List of atoms that have different ring attributes.
        """
        remove_atoms_template = []
        for mcs_pair in mcs_atoms:
            atom0 = mcs_pair[0]
            atom1 = mcs_pair[1]
            if (
                atom0.is_cyclic != atom1.is_cyclic
            ):  # and len(atom1.rings[0]) < 8: # allow partial ring matches for macrocycles
                remove_atoms_template.append(atom0)
            elif atom0.is_spiro != atom1.is_spiro:
                remove_atoms_template.append(atom0)
            elif atom0.is_cyclic:
                if len(atom0.rings) != len(atom1.rings):
                    remove_atoms_template.append(atom0)

        return remove_atoms_template

    def _mcs_partial_ring_match_from_bond(self, mcs_bonds, strict=False):
        remove_atoms_template = []
        safe_atoms = []
        for mcs_pair in mcs_bonds:
            bond0 = mcs_pair[0]
            bond1 = mcs_pair[1]
            if bond0.is_cyclic != bond1.is_cyclic:
                if strict:
                    atom0 = mcs_pair[0].atoms[0]
                    if atom_has_open_valency(atom0):
                        remove_atoms_template.append(atom0)
                    atom1 = mcs_pair[0].atoms[1]
                    if atom_has_open_valency(atom1):
                        remove_atoms_template.append(atom1)
                else:
                    remove_atoms_template.append(bond0.atoms[0])
                    remove_atoms_template.append(bond0.atoms[1])
                    if [a for a in bond0.atoms if a.atomic_symbol == "H"]:
                        safe_atoms.append(bond0.atoms[0].label)
                        safe_atoms.append(bond0.atoms[1].label)
        remove_atoms_template = [
            a for a in remove_atoms_template if a.label not in safe_atoms
        ]

        for mcs_pair in mcs_bonds:
            if not mcs_pair[0].is_cyclic and "ar" in mcs_pair[0].sybyl_type:
                atom0 = mcs_pair[0].atoms[0]
                if atom_has_open_valency(atom0):
                    remove_atoms_template.append(atom0)
                atom1 = mcs_pair[0].atoms[1]
                if atom_has_open_valency(atom1):
                    remove_atoms_template.append(atom1)
        return remove_atoms_template


def is_in_ring_with_open_valency(atom, mol):
    """
    Check if an atom belongs to a "flat" ring system that has an open valency
    :param atom:
    :param mol:
    :return:
    """

    if atom not in mol.atoms:
        return Exception("Atom is not part of Molecule.")
    if atom.is_cyclic:
        for ring in mol.rings:
            # Don't rotate saturated rings larger than 5, as their conformation is not sampled
            if ring.is_aromatic or ring.is_fully_conjugated or len(ring.atoms) == 5:
                if atom in ring.atoms:
                    for ring_atom in ring.atoms:
                        if atom_has_open_valency(ring_atom):
                            print("open valency")
                            return True
                        if ring.is_fused:
                            for ring_atom in ring.atoms:
                                fused_rings = [
                                    fused_ring
                                    for fused_ring in ring_atom.rings
                                    if fused_ring != ring
                                ]
                                if len(fused_rings) > 1:
                                    for fused_ring in fused_rings:
                                        for fused_ring_atom in fused_ring.atoms:
                                            if atom_has_open_valency(fused_ring_atom):
                                                print("Fused ring open valency")
                                                return True
    return False


def atom_has_open_valency(atom):
    num_neighbours = len(atom.neighbours)
    if atom.sybyl_type == "C.3" and num_neighbours < 4:
        return True
    if atom.sybyl_type == "C.2" and num_neighbours < 3:
        return True
    if atom.sybyl_type == "N.pl3" and num_neighbours < 3:
        return True
    if atom.sybyl_type == "N.3" and atom.formal_charge == 0 and num_neighbours < 3:
        return True
    if atom.sybyl_type == "N.am" and atom.formal_charge == 0 and num_neighbours < 3:
        return True
    if atom.sybyl_type == "N.3" and atom.formal_charge == 1 and num_neighbours < 4:
        return True
    if atom.sybyl_type == "N.4" and num_neighbours < 4:
        return True
    return False


class SubstitutentComparer(object):
    def __init__(self):
        self.smarts_dict = {"amide": "C(=O)N"}

    def compare(self, ligand, mcs_atoms):
        """
        Return atoms that should be removed from strict MCS because their number of heavy atom neighbours is
        different from ligand. Also removes partially matched amide carbonyls.
        :param ligand:
        :param mcs_atoms:
        :return: list of atoms
        """
        remove_atoms_template = []
        for substructure in self.smarts_dict:
            smarts = self.smarts_dict[substructure]
            searcher = search.SubstructureSearch()
            searcher.add_substructure(search.SMARTSSubstructure(smarts))
            ligand_hits = searcher.search(ligand)

            for hit in ligand_hits:
                ligand_match_atoms = hit.match_atoms()
                ligand_mcs_pairs = [a for a in mcs_atoms if a[1] in ligand_match_atoms]
                for mcs_pair in ligand_mcs_pairs:
                    scaffold_atom = mcs_pair[0]
                    heavy_scaffold_neighbours = [
                        a for a in scaffold_atom.neighbours if a.atomic_symbol != "H"
                    ]
                    heavy_ligand_neighbours = [
                        a for a in mcs_pair[1].neighbours if a.atomic_symbol != "H"
                    ]

                    if scaffold_atom.atomic_symbol == "N":
                        # tolerance for macrocycles
                        if scaffold_atom.rings and len(scaffold_atom.rings[0]) > 10:
                            continue

                        # for partial amide matches remove entire amide, to enable its rotation
                        if len(heavy_scaffold_neighbours) != len(
                            heavy_ligand_neighbours
                        ):
                            remove_atoms_template.append(scaffold_atom)
                            if len(ligand_mcs_pairs) == len(ligand_match_atoms):
                                remove_atoms_template.append(ligand_mcs_pairs[1][0])
                        # N(-C)-C will lead to symmetric matches. If both C have open valency, strict scaffold match can be off
                        elif len(heavy_scaffold_neighbours) == 3:
                            open_valency_neigbours = [
                                n
                                for n in heavy_scaffold_neighbours
                                if atom_has_open_valency(n)
                            ]
                            if len(open_valency_neigbours) == 2:
                                remove_atoms_template.extend(open_valency_neigbours)
        return remove_atoms_template


def _compare_mcs_stereo_chemistry(scaffold, ligand, mcs_atoms) -> bool:
    """
    :return: True if MCS have the same stereochemistry, else return False
    """

    ligand = ligand.copy()

    # mcs_atoms, mcs_bonds = mcs_searcher.search(
    #     scaffold, ligand, search_step_limit=1000000
    # )
    if not mcs_atoms:
        return []
    ligand_mcs_atoms = [a[1] for a in mcs_atoms]
    ligand_remove_atoms = [a for a in ligand.atoms if a not in ligand_mcs_atoms]
    ligand.remove_atoms(ligand_remove_atoms)

    # rdkit_scaffold
    rdkit_scaffold = Chem.MolFromMol2Block(scaffold.to_string("mol2"), removeHs=False)
    params = Chem.RemoveHsParameters()
    params.removeDegreeZero = True
    rdkit_scaffold = Chem.RemoveHs(rdkit_scaffold, params)

    # rdkit_ligand
    rdkit_ligand = Chem.MolFromMol2Block(ligand.to_string("mol2"))

    Chem.rdmolops.FindPotentialStereo(rdkit_scaffold)
    Chem.rdmolops.FindPotentialStereo(rdkit_ligand)

    (
        scaffold_atom_stereo,
        scaffold_atom_hybridization,
    ) = _return_stereocenter_hybridization_dict(rdkit_scaffold)
    (
        ligand_atom_stereo,
        ligand_atom_hybridization,
    ) = _return_stereocenter_hybridization_dict(rdkit_ligand)

    template_atoms_stereo_mismatch = []
    hybridization_mismatch = []
    for mcs_atom_pair in mcs_atoms:
        if mcs_atom_pair[0].atomic_symbol == "H":
            continue
        ligand_stereo = None
        scaffold_stereo = None
        scaffold_atom = mcs_atom_pair[0]
        ligand_atom = mcs_atom_pair[1]
        scaffold_atom_label = scaffold_atom.label
        ligand_atom_label = ligand_atom.label

        if scaffold_atom_label in scaffold_atom_stereo.keys():
            scaffold_stereo = scaffold_atom_stereo[scaffold_atom_label]
        if ligand_atom_label in ligand_atom_stereo.keys():
            ligand_stereo = ligand_atom_stereo[ligand_atom_label]
        if (
            scaffold_stereo is not None
            and ligand_stereo is not None
            and scaffold_stereo != ligand_stereo
        ):
            ligand_neighbours = sorted(
                a.atomic_symbol
                for a in list(
                    itertools.chain.from_iterable(
                        n.neighbours for n in ligand_atom.neighbours
                    )
                )
            )
            scaffold_neighbours = sorted(
                a.atomic_symbol
                for a in list(
                    itertools.chain.from_iterable(
                        n.neighbours for n in scaffold_atom.neighbours
                    )
                )
            )
            if len(ligand_neighbours) > len(scaffold_neighbours):
                template_atoms_stereo_mismatch.append(scaffold_atom)
            elif ligand_neighbours == scaffold_neighbours:
                template_atoms_stereo_mismatch.append(scaffold_atom)
        if (
            ligand_atom_hybridization[ligand_atom_label]
            != scaffold_atom_hybridization[scaffold_atom_label]
        ):
            hybridization_mismatch.append(scaffold.atom(scaffold_atom_label))
    return list(set(template_atoms_stereo_mismatch + hybridization_mismatch))


def _return_stereocenter_hybridization_dict(rdkit_mol):
    stereo_center_dict = {}
    hybridization_dict = {}
    for a in rdkit_mol.GetAtoms():
        hybridization_dict[a.GetProp("_TriposAtomName")] = a.GetHybridization()
        if "_CIPCode" in a.GetPropsAsDict().keys():
            stereo_center_dict[a.GetProp("_TriposAtomName")] = a.GetProp("_CIPCode")
    return stereo_center_dict, hybridization_dict


def _return_cis_trans_bicycle(atom):
    atom_hydrogen = [a for a in atom.neighbours if a.is_cyclic == False][0]
    for neigh in atom.neighbours:
        if len(neigh.rings) == len(atom.rings):
            if len(neigh.neighbours) < 4:
                return False
            neigh_hydrogen = [a for a in neigh.neighbours if a.is_cyclic == False][0]
            if (
                descriptors.MolecularDescriptors.atom_distance(
                    neigh_hydrogen, atom_hydrogen
                )
                < 2.8
            ):
                return "cis"
            else:
                return "trans"


def _bicycle_cis_trans_mismatch(atom1, atom2):
    if atom1.sybyl_type.split(".")[1] in ["2", "ar"]:
        return False
    if atom2.sybyl_type.split(".")[1] in ["2", "ar"]:
        return False
    if atom1.atomic_symbol != "C":
        return False
    if atom2.atomic_symbol != "C":
        return False
    atom1_stereo = _return_cis_trans_bicycle(atom1)
    atom2_stereo = _return_cis_trans_bicycle(atom2)

    if atom1_stereo != atom2_stereo:
        return True
    else:
        return False


def _bicyclic_cis_trans_chirality(atom1, atom2):
    """Check if the bicyclic system has cis-trans chirality."""
    if (
        1 < len(atom1.rings) == len(atom2.rings)
        and atom1.is_spiro == False
        and atom2.is_spiro == False
    ):
        atom1_second_bridge_atom = [n for n in atom1.neighbours if len(n.rings) > 1][0]
        atom2_second_bridge_atom = [n for n in atom2.neighbours if len(n.rings) > 1][0]
        if (
            len(atom1_second_bridge_atom.neighbours)
            == len(atom2_second_bridge_atom.neighbours)
            == 4
        ):
            return True
        else:
            return False
    else:
        return False


def assign_mcs_results_to_other_object(new_molecule, old_mcs_atoms, old_mcs_bonds):
    mcs_atoms = []
    scaffold_atom_labels = [a.label for a in new_molecule.atoms]
    for sa, qa in old_mcs_atoms:
        if sa.label in scaffold_atom_labels:
            mcs_atoms.append((new_molecule.atom(sa.label), qa))

    mcs_bonds = []
    for sb, qb in old_mcs_bonds:
        if (
            sb.atoms[0].label in scaffold_atom_labels
            and sb.atoms[1].label in scaffold_atom_labels
        ):
            _b = [
                b
                for b in new_molecule.bonds
                if b.atoms[0].label == sb.atoms[0].label
                and b.atoms[1].label == sb.atoms[1].label
            ][0]
            mcs_bonds.append((_b, qb))
    return mcs_atoms, mcs_bonds


def get_macrocycle_peptide_bonds(ccdc_molecule):
    scaffold_ring = max(ccdc_molecule.rings, key=lambda x: len(x.atoms))
    amides_smarts = search.SMARTSSubstructure("NC=O")
    searcher = search.SubstructureSearch()
    searcher.add_substructure(amides_smarts)
    hits = searcher.search(ccdc_molecule)
    scaffold_ring_atoms = scaffold_ring.atoms
    scaffold_ring_bonds = scaffold_ring.bonds
    peptide_bonds = []
    for hit in hits:
        match_atoms = hit.match_atoms()[:2]
        if (
            match_atoms[0] in scaffold_ring_atoms
            and match_atoms[1] in scaffold_ring_atoms
        ):
            for b in scaffold_ring_bonds:
                bond_atoms = b.atoms
                if match_atoms[0] in bond_atoms and match_atoms[1] in bond_atoms:
                    peptide_bonds.append(b)
    return peptide_bonds


def macrocycle_mcs(template_mol, query_mol, mcs_searcher):
    mcs_searcher.settings.check_hydrogen_count = True
    template_peptide_bonds = get_macrocycle_peptide_bonds(template_mol)

    open_templates = []
    for peptide_bond in template_peptide_bonds:
        open_template = template_mol.copy()
        _peptide_bond = [
            b
            for b in open_template.bonds
            if b.atoms[0].label == peptide_bond.atoms[0].label
            and b.atoms[1].label == peptide_bond.atoms[1].label
        ][0]
        open_template.remove_bond(_peptide_bond)
        open_templates.append(open_template)

    ligand_peptide_bonds = get_macrocycle_peptide_bonds(query_mol)
    remove_bond = ligand_peptide_bonds[0]
    remove_bond_atoms = remove_bond.atoms
    remove_bond_type = remove_bond.bond_type
    query_mol.remove_bond(remove_bond)

    mcs_atoms = []
    mcs_bonds = []
    for open_template in open_templates:
        t1 = time.time()
        _mcs_atoms, _mcs_bonds = mcs_searcher.search(
            open_template, query_mol, search_step_limit=500000
        )
        t2 = time.time()
        if len(_mcs_atoms) > len(mcs_atoms):
            mcs_atoms, mcs_bonds = _mcs_atoms, _mcs_bonds
            print(len(mcs_atoms), t2 - t1)
    query_mol.add_bond(remove_bond_type, *remove_bond_atoms)
    mcs_atoms, mcs_bonds = assign_mcs_results_to_other_object(
        template_mol, mcs_atoms, mcs_bonds
    )
    return mcs_atoms, mcs_bonds
