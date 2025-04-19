from ccdc import io
from rf_statistics import los_descriptors
from ccdc_roche_scoring import pose_processor
import argparse
from pathlib import Path


def parse_args():
    """Define and parse the arguments to the script."""
    parser = argparse.ArgumentParser(
        description="""
        Pose check standalone script. Writes SDF with additional attributes: RF bins, intramolecular and intermolecular clash count, cis-amide torsion count.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,  # To display default values in help message.
    )

    parser.add_argument(
        "-l",
        "--ligands",
        help="SDF with ligands",
        required=True,
    )

    parser.add_argument(
        "-g",
        "--gold_conf",
        help="Path to gold.conf file.",
        default="dock_dir",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Path to output file.",
        default="",
        required=False,
    )

    return parser.parse_args()


def add_pose_metrics(docked_ligand_file, gold_conf):
    with io.EntryReader(docked_ligand_file) as rdr:
        assert len(rdr) == 1
        ccdc_entry = rdr[0]
        describer = los_descriptors.CsdDescriptorsFromGold(
            docked_ligand_file,
            gold_conf=gold_conf,
            only_binding_site=False,
        )
        contact_df = describer.contact_df()
        rf_count_df = los_descriptors.rf_count_df(contact_df, ccdc_entry.molecule)
        rf_count_df = rf_count_df.iloc[:, 20:].drop(columns="smiles")
        rf_count_df.columns = [str(c) for c in rf_count_df.columns]
        rf_count_df = rf_count_df.rename({"clash_count": "P-L_steric_clashes"})
        clash_dict = pose_processor.count_bad_contacts(ccdc_entry.molecule)
        ccdc_entry.attributes.update(clash_dict)
        ccdc_entry.attributes.update(rf_count_df.iloc[0].to_dict())
    return ccdc_entry


def main():
    args = parse_args()
    ccdc_entry = add_pose_metrics(args.ligands, args.gold_conf)
    if args.output:
        output_path = args.output
    else:
        output_path = Path(args.ligands).stem + "_pose_check.sdf"
    with io.EntryWriter(output_path) as w:
        w.write(ccdc_entry)


if __name__ == "__main__":
    main()
