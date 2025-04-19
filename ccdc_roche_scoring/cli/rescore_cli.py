from ccdc_roche_scoring.docking import gold_rescoring
from ccdc_roche_scoring.cli import pose_check_cli

import click
from pathlib import Path


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-l",
    "--ligand_file",
    default="ligand.sdf",
    help="SDF with one ligand structure.",
    required=False,
)
@click.option(
    "-p",
    "--protein_file",
    help="Protein file.",
    default="protein.mol2",
    required=True,
    type=str,
)
def rescore(ligand_file, protein_file):
    gold_rescoring(ligand_file, protein_file)
    pose_check_cli.pose_check(
        [
            "--gold_conf",
            "gold_rescoring.conf",
            "--output",
            "rescored_pose_checked_ligand.sdf",
            "--ligands",
            "rescored_ligand.sdf",
        ]
    )
    return


def main():
    rescore()


if __name__ == "__main__":
    main()
