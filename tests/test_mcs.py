from pathlib import Path

import pytest
from ccdc import io

from ccdc_roche_scoring.docking import _mcs_templates_df
from ccdc_roche_scoring.template_docking_mcs import MCS
from tests import template_docking_mcs_testdata


@pytest.mark.parametrize(
    "ligand_file, native_ligand_file,template_file",
    [
        ("test_ligand_1.sdf", "1faoj.sdf", "test_strict_template_1.sdf"),
        ("test_ligand_2.sdf", "8s9a.sdf", "test_strict_template_2.sdf"),
    ],
)
def test_template(ligand_file, native_ligand_file, template_file):
    ligand_file = str(Path(template_docking_mcs_testdata.__file__).parent / ligand_file)
    native_ligand_file = str(
        Path(template_docking_mcs_testdata.__file__).parent / native_ligand_file
    )
    template_file = str(
        Path(template_docking_mcs_testdata.__file__).parent / template_file
    )

    ligand_mol = io.MoleculeReader(ligand_file)[0]
    native_ligand_mol = io.MoleculeReader(native_ligand_file)[0]
    test_template_mol = io.MoleculeReader(template_file)[0]
    mcs = MCS(ligand_mol, native_ligand_mol)
    template = mcs.return_mcs_scaffold(
        partial_ring_matches_allowed=False, ignore_hydrogens=False
    )[0]
    assert len(template.heavy_atoms) == len(test_template_mol.heavy_atoms)


@pytest.mark.parametrize(
    "ligand_file, native_ligand_file, template_file",
    [
        ("test_ligand_1.sdf", "1faoj.sdf", "test_template_df.sdf"),
        ("test_ligand_2.sdf", "8s9a.sdf", "test_strict_template_2.sdf"),
    ],
)
def test_template_df(ligand_file, native_ligand_file, template_file):
    ligand_file = str(Path(template_docking_mcs_testdata.__file__).parent / ligand_file)
    native_ligand_file = str(
        Path(template_docking_mcs_testdata.__file__).parent / native_ligand_file
    )
    template_file = str(
        Path(template_docking_mcs_testdata.__file__).parent / template_file
    )
    test_template_mol = io.MoleculeReader(template_file)[0]

    native_entries = io.EntryReader(native_ligand_file)
    ligand_mol = io.MoleculeReader(ligand_file)[0]

    templates_df = _mcs_templates_df(native_entries, ligand_mol)

    template = templates_df["scaffold"].values[0]
    assert len(template.heavy_atoms) == len(test_template_mol.heavy_atoms)


def main():
    test_template("test_ligand_2.sdf", "8s9a.sdf", "test_strict_template_2.sdf")
    test_template_df("test_ligand_2.sdf", "8s9a.sdf", "test_strict_template_2.sdf")


if __name__ == "__main__":
    main()
