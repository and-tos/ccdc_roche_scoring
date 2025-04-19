from pathlib import Path
import click
from ccdc import io
from ccdc_roche_scoring import q_values
import pandas as pd
from scipy.stats.mstats import gmean
import numpy as np
from collections import defaultdict


def binned_q_values(_q_values=[]):
    bins = [0, 0.1, 0.3, 0.5, 0.8, 1.01]
    binned_contacts = (
        pd.cut(_q_values, bins)
        .value_counts()
        .to_frame()
        .sort_index()
        .transpose()
        .reset_index(drop=True)
    )
    return binned_contacts


def _q_value_gmean(mol_rdr):
    q_gmean_dict = defaultdict(list)
    bins = []
    for ccdc_mol in mol_rdr:
        df = q_values.QValue().get_q_values(ccdc_mol)
        q_gmean = np.nan
        if df.shape[0] > 0:
            _q_values = []
            for c in df.columns:
                if c.startswith("q") and not c.startswith("q_"):
                    _q_values.extend(df[df[c].isna() == False][c].values)
            q_gmean = gmean(_q_values)
            bins.append(binned_q_values(_q_values))
        else:
            bins.append(binned_q_values())
        q_gmean_dict["identifier"].append(ccdc_mol.identifier)
        q_gmean_dict["q_gmean"].append(q_gmean)
    df = pd.DataFrame(q_gmean_dict)
    bins = pd.concat(bins, ignore_index=True)
    df = pd.concat([df, bins], axis=1)
    return df


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-isdf",
    "--isdf",
    help="SDF with one ligand structure.",
    required=True,
)
@click.option(
    "-o",
    "--output",
    help="Path to output directory.",
    default="q_values.csv",
    required=False,
)
@click.option(
    "--gmean",
    help="Calculate geometric mean only.",
    is_flag=True,
    default=False,
    required=False,
)
def q_values_csv(isdf: str, output: str, gmean: bool = False):
    mol_rdr = io.MoleculeReader(isdf)
    if gmean:
        df = _q_value_gmean(mol_rdr)
    else:
        q_values_dfs = []
        for mol_cnt, ccdc_mol in enumerate(mol_rdr):
            df = q_values.QValue().get_q_values(ccdc_mol)
            if df.shape[0] > 0:
                dfs = []
                if "q2" in df.columns:
                    _df = df[df["q2"].isna() == False][
                        ["TOR1_atom1", "TOR1_atom2", "q2", "substructure_name", "TOR1"]
                    ]
                    _df["torsion_name"] = "TOR1"
                    dfs.append(_df)

                    _df = df[df["q2"].isna() == False][
                        ["TOR2_atom1", "TOR2_atom2", "q2", "substructure_name", "TOR2"]
                    ]
                    _df["torsion_name"] = "TOR2"
                    dfs.append(_df)

                if "q1" in df.columns:
                    _df = df[df["q1"].isna() == False][
                        ["TOR1_atom1", "TOR1_atom2", "q1", "substructure_name", "TOR1"]
                    ]
                    _df["torsion_name"] = "TOR1"
                    dfs.append(_df)

                for _df in dfs:
                    _df.columns = [
                        "TOR_atom1",
                        "TOR_atom2",
                        "q",
                        "substructure_name",
                        "TOR",
                        "torsion_name",
                    ]
                df = pd.concat(dfs, ignore_index=True)
                df[["TOR_atom1", "TOR_atom2"]] = [
                    sorted(t) for t in (zip(df["TOR_atom1"], df["TOR_atom2"]))
                ]
                df = df.drop_duplicates(["TOR_atom1", "TOR_atom2"])
            df["identifier"] = ccdc_mol.identifier
            df["mol_cnt"] = mol_cnt
            q_values_dfs.append(df)
        df = pd.concat(q_values_dfs, ignore_index=True)
    output = Path(output)
    df.to_csv(output, index=False, sep="\t")
    print(df)
    print("Output: ", output)
    return


if __name__ == "__main__":
    q_values_csv()
