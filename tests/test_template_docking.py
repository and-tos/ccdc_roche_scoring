from pathlib import Path
import pytest
from ccdc import io

from tests import template_docking_testdata

from pathlib import Path

from ccdc import io, protein
from ccdc.docking import Docker

from ccdc_roche_scoring.misc.template_dock_standalone import _dock
from ccdc_roche_scoring.docking import _mcs_metrics
import tempfile


@pytest.mark.parametrize(
    "protein_file, ligand_file, template_file",
    [("7prm_apo.mol2", "7prm_native_ligand.sdf", "7prm_template.sdf")],
)
def test_dock(protein_file, ligand_file, template_file):
    protein_file = str(Path(template_docking_testdata.__file__).parent / protein_file)
    ligand_file = str(Path(template_docking_testdata.__file__).parent / ligand_file)
    template_file = str(Path(template_docking_testdata.__file__).parent / template_file)

    ccdc_protein = protein.Protein.from_entry(io.EntryReader(protein_file)[0])
    native_ligand = io.MoleculeReader(ligand_file)[0]
    ligand_filename = Path(ligand_file)
    docking_folder = "."
    strucid = ccdc_protein.identifier
    ligand_id = "ligand"
    water_paths = None
    ccdc_template = io.MoleculeReader(template_file)[0]

    docking_folder = tempfile.TemporaryDirectory()
    _dock(
        Docker(),
        Path(protein_file),
        ccdc_protein,
        native_ligand,
        ligand_filename,
        docking_folder.name,
        strucid,
        ligand_id,
        water_paths,
        ccdc_template,
        ccdc_template.copy(),
        False,
        True,
    )

    docking_solns = list(
        (Path(docking_folder.name) / "dock_0").glob("gold_soln_input_ligand_m1_*.sdf")
    )
    scaffold = io.MoleculeReader(template_file)[0]
    rmsds = []
    for docking_soln in docking_solns:
        dockin_soln_entry = io.EntryReader(str(docking_soln))[0]
        docked_ligand_entry, docked_ligand = _mcs_metrics(dockin_soln_entry, scaffold)
        rmsds.append(docked_ligand_entry.attributes["RMSD_to_mcs"])
    assert min(rmsds) < 1
    assert len(docking_solns) == 3
    docking_folder.cleanup()
