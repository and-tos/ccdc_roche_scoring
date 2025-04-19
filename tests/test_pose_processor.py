from pathlib import Path

import pandas as pd
import pytest
from ccdc import io

from ccdc_roche_scoring.pose_processor import count_bad_contacts
from tests import pose_processor_testdata


@pytest.mark.parametrize(
    "file, acc, don, steric, cis_amide_count",
    [
        ("acc_arom_oxygen_no_clash.sdf", 0, 0, 0, 0),
        ("acc_acc_no_clash.sdf", 0, 0, 0, 0),
        ("don_don_clash.sdf", 0, 1, 0, 0),
        ("no_clash.sdf", 0, 0, 0, 0),
        ("acc_weak_acc_O_no_clash.sdf", 0, 0, 0, 0),
        ("acc_weak_acc_S_no_clash.sdf", 0, 0, 0, 0),
        ("acc_weak_acc_S_no_clash_2.sdf", 0, 0, 0, 0),
        ("acc_strong_acc_clash.sdf", 1, 0, 0, 0),
        ("acc_acc_no_clash_rotatable_bond_edge.sdf", 0, 0, 0, 0),
        ("intram_hbond_no_clash.sdf", 0, 0, 0, 0),
        ("dimethoxybenzen_no_acc_clash.sdf", 0, 0, 0, 0),
        ("thiophene_no_acc_clash.sdf", 0, 0, 0, 0),
        ("amide_no_steric_clash.sdf", 0, 0, 0, 0),
        ("LOMRUD_gauche_oxygen.sdf", 0, 0, 0, 0),
        ("LESXAJ.sdf", 0, 0, 0, 1),
        ("0atx0_001.sdf", 0, 0, 0, 0),
        ("1bwit_004.sdf", 0, 0, 1, 0),
        ("1bwdl_001.sdf", 0, 0, 1, 0),
        ("19PK1_001.sdf", 0, 0, 0, 0),
        ("1AVBP_008.sdf", 0, 0, 0, 0),
        ("3nlr_004.sdf", 0, 0, 0, 0),
        ("4f3c_011.sdf", 0, 0, 0, 0),
        ("1ck0_001.sdf", 2, 1, 1, 0),
        ("1vwbv_001.sdf", 0, 0, 0, 0),
        ("1lkrv_001.sdf", 1, 0, 2, 0),
    ],
)
def test_count_bad_contacts(file, acc, don, steric, cis_amide_count):
    testfile = str(Path(pose_processor_testdata.__file__).parent / file)
    # testfile = "tests/pose_processor_testdata" + file
    rdr = io.EntryReader(testfile)
    bad_contacts = count_bad_contacts(rdr[0].molecule)
    assert bad_contacts == {
        "acc_acc_contact": acc,
        "don_don_contact": don,
        "intramolecular_steric_clashes": steric,
        "cis_amide_count": cis_amide_count,
    }


df = pd.read_csv(str(Path(pose_processor_testdata.__file__).parent / "csd_test.csv"))


@pytest.mark.parametrize(
    ",".join(df.columns.values), [tuple(row) for _, row in df.iterrows()]
)
def test_count_bad_contacts_csd(csd_code, acc, don, steric, cis_amide_count):
    rdr = io.MoleculeReader("CSD")
    bad_contacts = count_bad_contacts(rdr.molecule(csd_code))
    assert bad_contacts == {
        "acc_acc_contact": acc,
        "don_don_contact": don,
        "intramolecular_steric_clashes": steric,
        "cis_amide_count": cis_amide_count,
    }


def main():
    # test_dict = {"HAQFAK": [0, 0, 0, 0]}
    # for csd_code in test_dict:
    #     test_count_bad_contacts_csd(csd_code, *test_dict[csd_code])

    test_dict = {
        "1lkrv_001.sdf": [0, 0, 0, 0],
    }
    for file in test_dict:
        test_count_bad_contacts(file, *test_dict[file])

    return


if __name__ == "__main__":
    main()
