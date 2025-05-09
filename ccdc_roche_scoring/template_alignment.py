#!/usr/bin/env python

import argparse
import time

from openeye import oechem, oeomega, oequacpac


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="""
        Execute Line of sight contact scripts.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument(
        "-l", "--ligand", help="Path to template_ligand.", default="tmp_ligand_mol.sdf"
    )

    parser.add_argument(
        "-t",
        "--template",
        help="Path to docking_ligand.",
        default="tmp_scaffold_mol.sdf",
    )

    parser.add_argument(
        "-o", "--output", help="Path to docking_ligand.", default="tmp_oe_ligand.mol2"
    )

    return parser.parse_args()


def generate_fraglib(mol, fname):
    flofs = oechem.oemolostream()
    if not flofs.open(fname):
        oechem.OEThrow.Fatal(f"Unable to open {fname} for writing fraglib")
    makefraglib = oeomega.OEMakeFragLib()
    makefraglib.ClearFragLibs()

    for frag in makefraglib.GetMissingFrags(mol):
        oechem.OEWriteMolecule(flofs, frag)

    flofs.close()


def gen_confs_from_mcs(
    mol,
    mcsmatch=None,
    maxintconfs=1,
    torlib="guba",
    timeout=10.0,
    energy_window=20,
    fixrms=0.15,
    atomexpr=oechem.OEExprOpts_AtomicNumber | oechem.OEExprOpts_Aromaticity,
    bondexpr=oechem.OEExprOpts_DefaultBonds,
    forcefield="mmff94smod_NoEstat",
):
    omegaOpts = oeomega.OEOmegaOptions()
    if torlib == "guba":
        torlib = oeomega.OETorLib(oeomega.OETorLibType_GubaV21)
    elif torlib == "original":
        torlib = oeomega.OETorLib(oeomega.OETorLibType_Original)
    else:
        print("Invalid selection for torsion library: " + str(torlib))
        exit(13)

    omegaOpts.SetTorLib(torlib)
    omegaOpts.SetStrictStereo(False)  # DBX make input option
    omegaOpts.SetMaxConfs(maxintconfs)
    omegaOpts.SetMaxSearchTime(timeout)
    omega = oeomega.OEOmega(omegaOpts)
    omega.SetMaxConfs(maxintconfs)
    omega.SetEnergyWindow(energy_window)
    omega.SetSearchForceField(forcefield)

    if mcsmatch:
        fixopts = oeomega.OEConfFixOptions()
        fixopts.SetFixMol(mcsmatch)
        fixopts.SetFixDeleteH(True)
        fixopts.SetAtomExpr(atomexpr)
        fixopts.SetBondExpr(bondexpr)
        omegaOpts.SetConfFixOptions(fixopts)
    omega = oeomega.OEOmega(omegaOpts)
    mol_with_confs = None

    ret_code = omega.Build(mol)
    if ret_code == oeomega.OEOmegaReturnCode_Success:
        mol_with_confs = mol

    # Try different fragment library formats.
    else:
        generate_fraglib(mcsmatch, "fraglib.oeb.gz")
        omega.AddFragLib("fraglib.oeb.gz")
        mol_with_confs = None
        ret_code = omega.Build(mol)
        if ret_code == oeomega.OEOmegaReturnCode_Success:
            mol_with_confs = mol
        else:
            # DBX TO DO ERROR CHECKING
            # Write out enantiomer for debugging
            oechem.OEThrow.Warning(
                "Enantiomer throws error %s" % (oeomega.OEGetOmegaError(ret_code))
            )
            pass

    if mol_with_confs:
        for dp in oechem.OEGetSDDataPairs(mol):
            oechem.OESetSDData(mol_with_confs, dp.GetTag(), dp.GetValue())
        if mcsmatch:
            oechem.OESetSDData(
                mol_with_confs, "MCS_Query", oechem.OEGetSDData(mcsmatch, "MCS_Query")
            )
    del omega
    if mol_with_confs:
        return mol_with_confs.GetConfs()


def _omega_overlay(
    tmp_ligand_file="tmp_ligand_mol.mol2",
    tmp_scaffold_file="tmp_scaffold_mol.mol2",
    output="tmp_oe_ligand.sdf",
):
    lig_ifs = oechem.oemolistream()
    if lig_ifs.open(tmp_ligand_file):
        oe_ligand = next(lig_ifs.GetOEMols())

    scaf_ifs = oechem.oemolistream()
    if scaf_ifs.open(tmp_scaffold_file):
        oe_scaffold = next(scaf_ifs.GetOEMols())

    confs = gen_confs_from_mcs(oe_ligand, oe_scaffold, maxintconfs=1)

    ofs = oechem.oemolostream()
    if not ofs.open(output):
        oechem.OEThrow.Fatal("unable to open output file tmp_oe_ligand.sdf")
    for conf in confs:
        oechem.OESuppressHydrogens(conf)
        oechem.OEAddExplicitHydrogens(conf)
        chargeEngine = oequacpac.OEFormalToPartialCharges()
        oequacpac.OEAssignCharges(conf, chargeEngine)
        oechem.OEWriteMolecule(ofs, conf)
        break
    ofs.close()
    return True


def main():
    args = parse_args()
    _omega_overlay(args.ligand, args.template, args.output)


if __name__ == "__main__":
    print("Running Omega overlay...")
    t1 = time.time()
    main()
    t2 = time.time()
    print("main took ", t2 - t1, " s.")
