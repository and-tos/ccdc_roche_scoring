import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib.colors import LogNorm
import matplotlib.colors as mcolors
from collections import defaultdict
from pathlib import Path
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit import Chem
from ccdc import io, descriptors
from ccdc_roche_utils.plot_utils import roche_branding
from ccdc_roche_utils.plot_utils.seaborn_plots import plot_unity
from rf_statistics import structure_quality
from rf_statistics.utils import run_multiprocessing

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Roche Sans Light"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Roche Sans Light",
        "mathtext.it": "Roche Sans Light:italic",
        "mathtext.bf": "Roche Sans Light:bold",
    }
)


def relative_rf_error(df):
    """
    Calculate error of relative RF through taylor expansion.

    f(x)=x/rf_pdb
    f'(x)=1/rf_pdb

    f(y)=rf_docking/y
    f'(y)=-rf_docking/(y^2)

    d_rf_rel = |1/rf_pdb|*d_rf_docking + |-rf_docking/(rf_pdb^2)|*d_rf_pdb
    """
    d_rf_docking = (df["rf_high_vina"] - df["rf_low_vina"]) / 2
    d_rf_pdb = (df["rf_high_pdb"] - df["rf_low_pdb"]) / 2
    d_rf_rel = (
        1 / df["rf_pdb"] * d_rf_docking + df["rf_vina"] / (df["rf_pdb"] ** 2) * d_rf_pdb
    )
    return d_rf_rel


def plot_distributions(
    df,
    ax=None,
    outfile="cumulative_distribution.tif",
    hue="Docking method",
):
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())

    f = sns.ecdfplot(
        data=df,
        x="RMSD [Å]",
        legend=True,
        hue=hue,
        palette=[
            roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
            roche_branding.RocheColours().roche_colours_rgb_dict["light_red"],
            roche_branding.RocheColours().roche_colours_rgb_dict["dark_red"],
            roche_branding.RocheColours().roche_colours_rgb_dict["orange"],
            roche_branding.RocheColours().roche_colours_rgb_dict["purple"],
        ],
        ax=ax,
    )
    f.set_xlim(-0.1, 8)
    f.xaxis.minorticks_on()
    f.xaxis.grid(which="major", color="w", linewidth=1.0)
    f.xaxis.grid(which="minor", color="w", linewidth=0.5)
    f.set(
        box_aspect=1,
    )
    if ax is None:
        f.figure.savefig(outfile, dpi=150, bbox_inches="tight")
        f.figure.clf()
    return


def histplot(df_list, x, hue, outfile):
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    df = pd.concat(df_list, ignore_index=True)
    if hue:
        f = sns.histplot(
            df,
            x=x,
            hue=hue,
            element="step",
            stat="density",
            common_norm=False,
            palette=[
                roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
                roche_branding.RocheColours().roche_colours_rgb_dict["dark_orange"],
            ],
        )
    else:
        f = sns.histplot(
            df,
            x=x,
            element="step",
            stat="density",
            common_norm=False,
            color=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
        )
    f.set(box_aspect=1)
    f.figure.savefig(outfile, dpi=300, bbox_inches="tight")
    f.figure.clf()
    return


def compare_high_low_rmsd(df):
    # low_rmsd = df[df["rmsd"] <= 1.0]
    # high_rmsd = df[df["rmsd"] > 1.0]
    # low_rmsd["RMSD"] = "<= 1.0 Å"
    # high_rmsd["RMSD"] = "> 1.0 Å"
    # histplot([low_rmsd, high_rmsd], "atom_rf_gmean", "RMSD", "high_vs_low_rmsd_rf")
    # histplot([low_rmsd, high_rmsd], "q_gmean", "RMSD", "high_vs_low_rmsd_q")
    df = df.rename(
        columns={"atom_rf_gmean": "Atom R$_{F}$ Gmean", "q_gmean": "q$_{gmean}$"}
    )
    df = df[df["dataset"] == "native"]
    histplot([df], "Atom R$_{F}$ Gmean", None, "atom_rf_gmean_dist")
    histplot([df], "q$_{gmean}$", None, "q_gmean_dist")

    df = df[(df["MW"] > 300) & (df["MW"] < 650)]
    df = df[(df["NumRotatableBonds"] < 10)]
    df = df[df["HBD"] < 4]
    df = df[df["HBA"] < 10]
    df = df[df["NumAromaticRings"] < 5]
    df = df[df["num_pos_charges"] < 2]
    df = df[df["num_neg_charges"] < 2]
    df = df[df["bs_id"].str.len() == 8]
    df = df[df["Atom R$_{F}$ Gmean"] < 0.5]
    return


def write_rmsd_df():
    # gold_ligands = io.MoleculeReader(
    #     [
    #         "../molecular_recognition_gold_2024-12-03/gold_top_docked_ligands_pub.csdsql",
    #         "../molecular_recognition_gold_2024-12-03/gold_top_docked_ligands_roche.csdsql",
    #     ]
    # )
    # vina_ligands = io.MoleculeReader(
    #     ["vina_top_docked_ligands_pub.csdsql", "vina_top_docked_ligands_roche.csdsql"]
    # )
    # reference_ligands = io.MoleculeReader(
    #     ["reference_ligands_pub.csdsql", "reference_ligands_roche.csdsql"]
    # )
    # rmsd_dict = defaultdict(list)
    # sdf = pd.read_parquet("pub_roche_project_annotation_2024-12-03.gz")
    # sdf = sdf.rename(columns={"bs_id": "identifier"})
    # struc_qual = structure_quality.StructureQuality(sdf)

    # for vina_ligand in vina_ligands:
    #     identifier = vina_ligand.identifier
    #     if identifier in struc_qual.low_quality_bs_ids:
    #         print(identifier, " is low quality.")
    #         continue
    #     try:
    #         reference = reference_ligands.molecule(identifier)
    #         gold_ligand = gold_ligands.molecule(identifier)
    #     except:
    #         continue
    #     gold_ligand.remove_hydrogens()
    #     vina_ligand.remove_hydrogens()
    #     reference.remove_hydrogens()
    #     vina_rmsd = descriptors.MolecularDescriptors.rmsd(vina_ligand, reference)
    #     gold_rmsd = descriptors.MolecularDescriptors.rmsd(gold_ligand, reference)
    #     rmsd_dict["identifier"].append(identifier)
    #     rmsd_dict["vina_rmsd"].append(vina_rmsd)
    #     rmsd_dict["gold_rmsd"].append(gold_rmsd)
    # df = pd.DataFrame(rmsd_dict)
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    df = pd.read_csv("rmsd.csv")
    df = df.rename(columns={"vina_rmsd": "Vina RMSD [Å]", "gold_rmsd": "GOLD RMSD [Å]"})
    rmsd_plot = rmsd_correlation_plot(df, ax=ax2)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.185, 0.02, 0.6])
    fig.colorbar(rmsd_plot, cax=cbar_ax)
    cbar_ax.set_ylabel("Number of binding sites")

    vina_df = df[["identifier", "Vina RMSD [Å]"]]
    gold_df = df[["identifier", "GOLD RMSD [Å]"]]

    vina_df["Docking method"] = "Vina"
    gold_df["Docking method"] = "GOLD"

    vina_df = vina_df.rename(columns={"Vina RMSD [Å]": "RMSD [Å]"})
    gold_df = gold_df.rename(columns={"GOLD RMSD [Å]": "RMSD [Å]"})

    rmsd_df = pd.concat([vina_df, gold_df], ignore_index=True)
    dist_plot = plot_distributions(rmsd_df, ax1)
    fig.figure.savefig("figure1.tif", dpi=150, bbox_inches="tight")
    fig.figure.clf()
    df.to_csv("rmsd.csv", index=False)

    return rmsd_df


def rmsd_comparison():
    rmsd_df = write_rmsd_df()
    # df = pd.read_csv("rmsd.csv")
    # # Add docking RMSD to los_df
    # for los_file in Path(".").glob("rf_statistics/*/query_atom/los*filtered.csv"):
    #     los_df = pd.read_csv(los_file)
    #     complex_df = pd.read_csv(los_file.parent / "complex_filtered.csv")
    #     los_df = los_df.drop(
    #         columns=["vina_rmsd", "Vina RMSD [Å]", "ligand_smiles"], errors="ignore"
    #     )
    #     los_df = los_df.join(
    #         df[["identifier", "Vina RMSD [Å]"]].set_index("identifier"),
    #         on="molecule_name",
    #     )

    #     los_df = los_df.join(
    #         complex_df[["molecule_name", "ligand_smiles"]].set_index("molecule_name"),
    #         on="molecule_name",
    #     )
    #     los_df["is_public"] = False
    #     los_df.loc[los_df["molecule_name"].str.len() == 8, "is_public"] = True
    #     los_df.to_csv(los_file, index=False)

    # # write docking RMSD data also to PDB files
    # for los_file in Path("../protonate3d_2024-12-03/rf_statistics").glob(
    #     "*/query_atom/los*filtered.csv"
    # ):
    #     los_df = pd.read_csv(los_file)
    #     complex_df = pd.read_csv(los_file.parent / "complex_filtered.csv")
    #     los_df = los_df.drop(
    #         columns=["vina_rmsd", "Vina RMSD [Å]", "ligand_smiles"], errors="ignore"
    #     )
    #     los_df = los_df.join(
    #         df[["identifier", "Vina RMSD [Å]"]].set_index("identifier"),
    #         on="molecule_name",
    #     )

    #     los_df = los_df.join(
    #         complex_df[["molecule_name", "ligand_smiles"]].set_index("molecule_name"),
    #         on="molecule_name",
    #     )
    #     los_df["is_public"] = False
    #     los_df.loc[los_df["molecule_name"].str.len() == 8, "is_public"] = True
    #     los_df.to_csv(los_file, index=False)

    # return rmsd_df


def rmsd_normalized_prop_plot(
    hist2d,
    x="Vina RMSD [Å]",
    y="GOLD RMSD [Å]",
    ax=None,
    axlim: tuple = (0, 20),
    cbar_ax=None,
):
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())

    # Define the colors for the colormap: light blue and darker blue
    dark_blue = roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"]
    light_blue = roche_branding.RocheColours().roche_colours_rgb_dict[
        "extra_light_blue"
    ]

    # Create a custom colormap
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "custom_cmap", [light_blue, dark_blue, "black"], N=512
    )
    plot = sns.heatmap(
        hist2d,
        cmap=cmap,
        ax=ax,
        vmin=0.01,
        vmax=1,
        norm="log",
        cbar=True,
        cbar_ax=cbar_ax,
    )
    # plot = ax.imshow(
    #     hist2d,
    #     cmap=cmap,
    #     vmin=0.01,
    #     vmax=1,
    #     norm="log",
    # )
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if axlim:
        ax.set_xlim(*axlim)
        ax.set_ylim(*axlim)
    ax.grid(True, color="lightgrey")
    ax.set_facecolor("white")
    # ax.set_aspect("equal")

    return plot


def rmsd_correlation_plot(
    df,
    x="Vina RMSD [Å]",
    y="GOLD RMSD [Å]",
    ax=None,
    axlim: tuple = (0, 20),
    unity: bool = True,
    outname: str = "rmsd_vina_vs_gold.tif",
):

    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())

    # Define the colors for the colormap: light blue and darker blue
    dark_blue = roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"]
    light_blue = roche_branding.RocheColours().roche_colours_rgb_dict[
        "extra_light_blue"
    ]

    # Create a custom colormap
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "custom_cmap", [light_blue, dark_blue, "black"], N=512
    )
    plot = ax.hexbin(
        df[x], df[y], gridsize=20, bins="log", cmap=cmap, vmin=1, vmax=100000
    )
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    if axlim:
        ax.set_xlim(*axlim)
        ax.set_ylim(*axlim)
    if unity:
        plot_unity(df[x], df[y])
    ax.grid(True, color="lightgrey")
    ax.set_facecolor("white")
    ax.set_aspect("equal")

    return plot


def rf_correlation_plot(df):
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())

    plot = sns.scatterplot(
        df,
        x="rf_pdb",
        y="rf_vina",
        palette=[
            roche_branding.RocheColours().roche_colours_rgb_dict["extra_light_blue"],
            roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
        ],
        hue="is significant",
        edgecolor="grey",
        s=8,
    )
    plot_unity(df["rf_pdb"], df["rf_vina"])
    plot.set(
        box_aspect=1,
    )
    xdata = df["rf_pdb"]
    ydata = df["rf_pdb"]
    plot.text(
        xdata.max() - 1,
        ydata.max() - 3,
        "Underpredicted by docking",
        fontsize=5,
        rotation=45,
    )
    plot.text(
        xdata.max() - 3,
        ydata.max() - 1,
        "Overpredicted by docking",
        fontsize=5,
        rotation=45,
    )
    plot.set_xscale("log")
    plot.set_yscale("log")
    # plt.arrow(
    #     xdata.max() - 2, ydata.max() - 2, -1.5, 1.5, color="black", head_width=0.1
    # )
    # plt.arrow(
    #     xdata.max() - 2, ydata.max() - 2, 1.5, -1.5, color="black", head_width=0.1
    # )
    plot.set_ylabel("R$_{F,ligand}^{docking}$")
    plot.set_xlabel("R$_{F,ligand}^{experimental}$")
    plot.figure.savefig("rf_vina_vs_pdb.tif", dpi=300, bbox_inches="tight")
    plot.figure.clf()
    return


def rf_comparison():
    dfs = []
    for rf_dir in Path("rf_statistics").glob("*"):
        vina_file = rf_dir / "query_atom" / "rf.csv"
        pdb_file = Path("../protonate3d_2024-12-03") / rf_dir / "query_atom" / "rf.csv"
        if vina_file.is_file() and pdb_file.is_file():
            ligand_atom_type = rf_dir.name
            vina_rf_df = pd.read_csv(vina_file).drop(
                columns=["type", "alpha_min", "alpha_max"]
            )
            pdb_rf_df = pd.read_csv(pdb_file).drop(
                columns=["type", "alpha_min", "alpha_max"]
            )
            df = vina_rf_df.join(
                pdb_rf_df.set_index("atom_type"),
                on="atom_type",
                lsuffix="_vina",
                rsuffix="_pdb",
            )
            df["relative_rf"] = df["rf_vina"] / df["rf_pdb"]
            df["relative_rf_error"] = relative_rf_error(df)
            df["ligand_atom_type"] = ligand_atom_type
            dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    df = df[df["relative_rf"].isna() == False]
    df = df[df["atom_type"].isin(["other_ligand", "metal"]) == False]
    df = df[(df["expected_vina"] > 10) & (df["size_vina"] > 10)]
    df["is significant"] = False

    df.loc[
        (df["relative_rf"] > 1) & ((df["relative_rf"] - (df["relative_rf_error"])) > 1),
        "is significant",
    ] = True

    df.loc[
        (df["relative_rf"] < 1) & ((df["relative_rf"] + (df["relative_rf_error"])) < 1),
        "is significant",
    ] = True

    rf_correlation_plot(df)
    return df


def plot_relative_rf(df):
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    df = df[df["relative_rf_error"] / df["relative_rf"] < 1]
    # df = df[df["is significant"] == True]
    df = df[df["relative_rf_error"].isna() == False]
    for ligand_atom_type, gdf in df.groupby("ligand_atom_type"):
        gdf = gdf.sort_values("relative_rf", ascending=False)
        plot = sns.barplot(
            gdf,
            x="atom_type",
            y="relative_rf",
            yerr=gdf["relative_rf_error"],
            color=roche_branding.RocheColours().roche_colours_rgb_dict["purple"],
        )
        plot.set_ylabel("rR$_{F,ligand}$")
        gdf["relative_rf_error"]
        plot.tick_params(axis="x", rotation=90)
        plot_path = Path(f"rf_statistics/{ligand_atom_type}/query_atom/relative_rf.tif")
        plot.figure.savefig(plot_path, dpi=300, bbox_inches="tight")
        plot.figure.clf()
    return


def rdkit_descs_from_smiles(smiles):
    if type(smiles) == str:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            atoms = mol.GetAtoms()
            formal_charges = [a.GetFormalCharge() for a in atoms]
            num_pos_charges = sum([fc for fc in formal_charges if fc > 0])
            num_neg_charges = sum([fc for fc in formal_charges if fc < 0])
            return (
                Chem.rdMolDescriptors.CalcExactMolWt(mol),
                Chem.rdMolDescriptors.CalcNumRotatableBonds(mol),
                Chem.rdMolDescriptors.CalcNumHBD(mol),
                Chem.rdMolDescriptors.CalcNumHBA(mol),
                Chem.rdMolDescriptors.CalcNumAromaticRings(mol),
                num_pos_charges,
                num_neg_charges,
            )

    return tuple([pd.NA for a in range(7)])


def getMolDescriptors(mol, missingVal=None):
    """calculate the full list of descriptors for a molecule

    missingVal is used if the descriptor cannot be calculated
    """
    res = {}
    for nm, fn in Descriptors._descList:
        # some of the descriptor fucntions can throw errors if they fail, catch those here:
        try:
            val = fn(mol)
        except:
            # print the error message:
            import traceback

            traceback.print_exc()
            # and set the descriptor value to whatever missingVal is
            val = missingVal
        res[nm] = val
    return res


def rdkit_descs_from_smiles_list(smiles_list):
    # dfs = []
    # for smiles in smiles_list:
    #     dfs.append(rdkit_descs_from_smiles(smiles))
    # pd.concat(dfs)
    return run_multiprocessing(24, smiles_list, rdkit_descs_from_smiles)


def return_recovery(df):
    df = df.drop_duplicates("bs_id")
    return 100 * df["dataset"].value_counts()["native"] / df.shape[0]


def return_recovery_dict(_, df, randomize=False):
    if randomize:
        random_bs_ids = pd.Series(df["bs_id"].unique()).sample(frac=1, replace=True)

    recovery = defaultdict(list)
    # RF recovery
    rf_df = df.sort_values("atom_rf_gmean", ascending=False).drop_duplicates("bs_id")
    if randomize:
        rf_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            rf_df.set_index("bs_id"), on="bs_id"
        )
    print(rf_df.shape[0])
    pct_recovery = return_recovery(rf_df)
    recovery["Scoring method"].append("Atom R$_{F}$ Gmean")
    recovery["% recovered experimental poses"].append(pct_recovery)

    # RF * q_gmean recovery
    rfq_df = df.sort_values(
        ["rfq", "atom_rf_gmean"], ascending=[False, False]
    ).drop_duplicates("bs_id")
    if randomize:
        rfq_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            rfq_df.set_index("bs_id"), on="bs_id"
        )
    print(rfq_df.shape[0])
    pct_recovery = return_recovery(rfq_df)
    recovery["Scoring method"].append("R$_{F}$Q")
    recovery["% recovered experimental poses"].append(pct_recovery)

    # q_gmean recovery
    q_gmean_df = df.sort_values(
        ["q_gmean", "vina_score"], ascending=[False, True]
    ).drop_duplicates("bs_id")
    if randomize:
        q_gmean_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            q_gmean_df.set_index("bs_id"), on="bs_id"
        )
    print(q_gmean_df.shape[0])
    pct_recovery = return_recovery(q_gmean_df)
    recovery["Scoring method"].append("q$_{gmean}$")
    recovery["% recovered experimental poses"].append(pct_recovery)

    # Vina recovery
    vina_df = df.sort_values("vina_score", ascending=True).drop_duplicates("bs_id")
    if randomize:
        vina_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            vina_df.set_index("bs_id"), on="bs_id"
        )
    print(vina_df.shape[0])
    pct_recovery = return_recovery(vina_df)
    recovery["Scoring method"].append("docking")
    recovery["% recovered experimental poses"].append(pct_recovery)

    # Random recovery
    random_df = df.sample(frac=1).drop_duplicates("bs_id")
    if randomize:
        random_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            random_df.set_index("bs_id"), on="bs_id"
        )
    print(random_df.shape[0])
    pct_recovery = return_recovery(random_df)
    recovery["Scoring method"].append("random")
    recovery["% recovered experimental poses"].append(pct_recovery)
    recovery = pd.DataFrame(recovery)
    return recovery, vina_df, rf_df, rfq_df, q_gmean_df, random_df


def return_sorted_dfs(_, df, randomize=False):
    if randomize:
        random_bs_ids = pd.Series(df["bs_id"].unique()).sample(frac=1, replace=True)

    # RF recovery
    rf_df = df.sort_values("atom_rf_gmean", ascending=False)
    if randomize:
        rf_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            rf_df.set_index("bs_id"), on="bs_id"
        )
    print(rf_df.shape[0])

    # RF * q_gmean recovery
    rfq_df = df.sort_values(["rfq", "atom_rf_gmean"], ascending=[False, False])
    if randomize:
        rfq_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            rfq_df.set_index("bs_id"), on="bs_id"
        )
    print(rfq_df.shape[0])

    # q_gmean recovery
    q_gmean_df = df.sort_values(["q_gmean", "vina_score"], ascending=[False, True])
    if randomize:
        q_gmean_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            q_gmean_df.set_index("bs_id"), on="bs_id"
        )
    print(q_gmean_df.shape[0])

    # Vina recovery
    vina_df = df.sort_values("vina_score", ascending=True)
    if randomize:
        vina_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            vina_df.set_index("bs_id"), on="bs_id"
        )
    print(vina_df.shape[0])

    # Random recovery
    random_df = df.sample(frac=1)
    if randomize:
        random_df = pd.DataFrame({"bs_id": random_bs_ids}).join(
            random_df.set_index("bs_id"), on="bs_id"
        )
    print(random_df.shape[0])

    vina_df["rescoring_method"] = "docking"
    rf_df["rescoring_method"] = "Atom R$_{F}$ Gmean"
    rfq_df["rescoring_method"] = "R$_{F}$Q"
    q_gmean_df["rescoring_method"] = "q$_{gmean}$"
    random_df["rescoring_method"] = "random"

    return vina_df, rf_df, rfq_df, q_gmean_df, random_df


def return_recovery_plot(recovery, outname, ax=None):
    # generate plot
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    f = sns.barplot(
        recovery,
        x="Scoring method",
        y="% recovered experimental poses",
        color=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
        ax=ax,
        edgecolor=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
        linewidth=0,
    )
    for index, row in recovery.iterrows():
        f.errorbar(
            row["Scoring method"],
            row["% recovered experimental poses"],
            yerr=[[row["low"]], [row["high"]]],
            fmt="none",
            c="black",
            capsize=5,
        )

    f.set(ylim=(0, 100), box_aspect=1)
    f.tick_params(axis="x", labelrotation=45)
    if ax is None:
        f.figure.savefig(
            f"experimental_pose_recovery_{outname}.tif", dpi=300, bbox_inches="tight"
        )
        f.figure.clf()
    return f


def _make_recovery_plot(df, outname, confidence_level=0.05):
    df["q_gmean"] = df["q_gmean"].fillna(1)
    df["rfq"] = df["atom_rf_gmean"] * df["q_gmean"].fillna(1)
    df = df[
        [
            "bs_id",
            "dataset",
            "q_gmean",
            "atom_rf_gmean",
            "rfq",
            "vina_score",
            "rmsd",
        ]
    ]
    recovery, vina_df, rf_df, rfq_df, q_gmean_df, random_df = return_recovery_dict(
        1, df, randomize=True
    )
    bootstrap_recovery = run_multiprocessing(
        24, range(5000), return_recovery_dict, df=df, randomize=True
    )
    bootstrap_recovery = pd.concat(
        [r[0] for r in bootstrap_recovery], ignore_index=True
    )

    for scoring_method, gdf in bootstrap_recovery.groupby("Scoring method"):
        recoveries = list(gdf["% recovered experimental poses"].sort_values())
        recovery_low = recoveries[int((confidence_level / 2.0) * len(recoveries))]
        recovery_high = recoveries[int((1 - confidence_level / 2.0) * len(recoveries))]
        recovery.loc[recovery["Scoring method"] == scoring_method, "high"] = (
            recovery_high
        )
        recovery.loc[recovery["Scoring method"] == scoring_method, "low"] = recovery_low

    recovery["low"] = recovery["% recovered experimental poses"] - recovery["low"]
    recovery["high"] = recovery["high"] - recovery["% recovered experimental poses"]

    sort_dict = {
        "random": 0,
        "docking": 1,
        "q$_{gmean}$": 2,
        "Atom R$_{F}$ Gmean": 3,
        "R$_{F}$Q": 4,
    }
    recovery = recovery.iloc[
        recovery["Scoring method"].map(sort_dict).sort_values().index
    ].reset_index()

    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12, 6))
    return_recovery_plot(recovery, outname, ax0)

    vina_df["rescoring_method"] = "docking"
    rf_df["rescoring_method"] = "Atom R$_{F}$ Gmean"
    rfq_df["rescoring_method"] = "R$_{F}$Q"
    q_gmean_df["rescoring_method"] = "q$_{gmean}$"
    random_df["rescoring_method"] = "random"

    rescoring_df = pd.concat(
        [rf_df, rfq_df, q_gmean_df, vina_df, random_df], ignore_index=True
    )
    rescoring_df = rescoring_df.rename(
        columns={"rmsd": "RMSD [Å]", "rescoring_method": "Scoring method"}
    )
    plot_distributions(
        rescoring_df, ax1, f"rescoring_distribution_{outname}.tif", "Scoring method"
    )
    fig.figure.savefig("figure10.tif", dpi=300, bbox_inches="tight")
    fig.figure.clf()
    return


def plot_experimental_pose_recovery():
    # # Native poses
    # native_poses = PandasTools.LoadSDF(
    #     "pose_check_native_cat.sdf",
    #     molColName=None,
    #     sanitize=False,
    #     smilesName="ligand_smiles",
    # )
    # native_vina_scores = pd.read_csv("native_pose_rescore.csv")
    # native_poses = native_poses.join(
    #     native_vina_scores[["identifier", "vina_native_score"]].set_index("identifier"),
    #     on="bs_id",
    # )
    # native_poses["dataset"] = "native"
    # native_poses = native_poses.rename(columns={"vina_native_score": "vina_score"})

    # # Docked poses
    # docked_poses = PandasTools.LoadSDF(
    #     "pose_check_vina_cat.sdf",
    #     molColName=None,
    #     sanitize=False,
    #     smilesName="ligand_smiles",
    # )
    # docked_poses["dataset"] = "vina"

    # # Combine both datasets
    # native_poses = native_poses[native_poses["bs_id"].isin(docked_poses["bs_id"])]
    # docked_poses = docked_poses[docked_poses["bs_id"].isin(native_poses["bs_id"])]
    # combined = pd.concat([docked_poses, native_poses])
    # combined = combined.drop(columns="ID")
    # combined["vina_score"] = combined["vina_score"].astype(float)
    # combined = pd.read_parquet("native_and_docking_pose_check.gzip")

    # combined.to_parquet("native_and_docking_pose_check.gzip")

    # descs = rdkit_descs_from_smiles_list(combined["ligand_smiles"].to_list())
    # combined[
    #     [
    #         "MW",
    #         "NumRotatableBonds",
    #         "HBD",
    #         "HBA",
    #         "NumAromaticRings",
    #         "num_pos_charges",
    #         "num_neg_charges",
    #     ]
    # ] = descs
    # combined = combined.astype(
    #     {"atom_rf_gmean": float, "q_gmean": float, "rmsd": float}
    # )
    # combined.loc[combined["dataset"] == "native", "rmsd"] = 0
    combined = pd.read_parquet("native_and_docking_pose_check_qmean.gzip")
    # compare_high_low_rmsd(combined)
    # _make_recovery_plot(combined, "all")

    # Filter
    druglike_df = combined[(combined["MW"] > 300) & (combined["MW"] < 650)]
    druglike_df = druglike_df[(druglike_df["NumRotatableBonds"] < 10)]
    druglike_df = druglike_df[druglike_df["HBD"] < 4]
    druglike_df = druglike_df[druglike_df["HBA"] < 10]
    druglike_df = druglike_df[druglike_df["NumAromaticRings"] < 5]
    druglike_df = druglike_df[druglike_df["num_pos_charges"] < 2]
    druglike_df = druglike_df[druglike_df["num_neg_charges"] < 2]

    _make_recovery_plot(druglike_df, "druglike")

    return


def bin_dataframe(df, y, bins: list = []):
    "Returns dataframe with x on columns y on index and frequencies in cells."
    df["vina_bins"] = pd.cut(
        df["Vina RMSD [Å]"], [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5, 16]
    )
    df["y_bins"] = pd.cut(df[y], bins)
    hist_df = pd.crosstab(index=df["y_bins"], columns=df["vina_bins"])
    hist_df = hist_df.reindex(index=hist_df.index[::-1])
    return hist_df


def props_vs_rmsd():
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())

    df = pd.read_parquet("native_and_docking_pose_check_qmean.gzip")
    df = df.sort_values("vina_score").drop_duplicates("bs_id")

    adf = pd.read_parquet("pub_roche_project_annotation_2024-12-03.gz")
    adf["identifier"] = adf["bs_id"]
    sq = structure_quality.StructureQuality(adf)
    low_quality_bs_ids = sq.low_quality_bs_ids
    adf = adf[["bs_id", "RESOLUTION"]]

    descriptors = [
        "NumRotatableBonds",
        "NumAromaticRings",
        "Molecular weight [g/mol]",
        "Resolution [Å]",
    ]
    df = df.join(adf.set_index("bs_id"), on="bs_id")
    df = df.rename(
        columns={
            "rmsd": "Vina RMSD [Å]",
            "RESOLUTION": "Resolution [Å]",
            "MW": "Molecular weight [g/mol]",
        }
    )
    df = df[["bs_id", "Vina RMSD [Å]"] + descriptors]
    df = df[df["bs_id"].isin(low_quality_bs_ids) == False]
    df = df[df["Resolution [Å]"] > 0]
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])

    for cnt, descriptor in enumerate(descriptors):
        bins = list(range(int(df[descriptor].max())))
        if descriptor == "Resolution [Å]":
            bins = [x / 10 for x in range(5, 30, 5)]
        if descriptor == "Molecular weight [g/mol]":
            bins = list(range(50, 1000, 50))
        hist_2d = bin_dataframe(df, descriptor, bins)
        for i in hist_2d.index:
            hist_2d.loc[i, :] = hist_2d.loc[i] / hist_2d.loc[i].sum()
        ax = axes.flatten()[cnt]
        rmsd_plot = rmsd_normalized_prop_plot(
            hist_2d, "Vina RMSD [Å]", descriptor, ax, axlim=(), cbar_ax=cbar_ax
        )

    cbar_ax.set_ylabel("Number of poses normalized by descriptor bin")
    cbar_ax.set_label("Number of poses normalized by descriptor bin")
    fig.subplots_adjust(hspace=0.3)
    fig.figure.savefig("figureS1.tif", dpi=100, bbox_inches="tight")
    fig.figure.clf()

    fig, axes = plt.subplots(2, 2, figsize=(12, 12))

    for cnt, descriptor in enumerate(descriptors):
        bins = list(range(int(df[descriptor].max())))
        if descriptor == "Resolution [Å]":
            bins = [x / 10 for x in range(5, 30, 5)]
        if descriptor == "Molecular weight [g/mol]":
            bins = list(range(50, 1000, 50))
        ax = axes.flatten()[cnt]
        sns.histplot(
            df,
            x=descriptor,
            bins=bins,
            color=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
            ax=ax,
        )
    fig.figure.savefig("figureS2.tif", dpi=100, bbox_inches="tight")
    return


def plot_native_pose_rank():
    df = pd.read_parquet("native_and_docking_pose_check_qmean.gzip")
    df["q_gmean"] = df["q_gmean"].fillna(1)
    df["rfq"] = df["atom_rf_gmean"] * df["q_gmean"].fillna(1)
    df = df[
        [
            "bs_id",
            "dataset",
            "q_gmean",
            "atom_rf_gmean",
            "rfq",
            "vina_score",
            "rmsd",
        ]
    ]
    dfs = return_sorted_dfs(1, df, randomize=False)
    df = pd.concat(dfs, ignore_index=True)
    df["Experimental pose rank"] = (
        df.groupby(["bs_id", "rescoring_method"]).cumcount().astype(int)
    )
    df = df[df["dataset"] == "native"]
    fig, axes = plt.subplots(3, 2, figsize=(18, 12))
    fig.subplots_adjust(hspace=0.3)
    axes = axes.flatten()
    bins = range(0, 21)
    for cnt, (rescoring_method, gdf) in enumerate(df.groupby("rescoring_method")):
        hist_plot = sns.histplot(
            gdf,
            x="Experimental pose rank",
            bins=range(22),
            color=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
            ax=axes[cnt],
        )
        hist_plot.set_title(rescoring_method, y=1.0, pad=-14)
        hist_plot.set_xticks([b + 0.5 for b in bins])
        hist_plot.set_xticklabels([b + 1 for b in bins])

    axes[-1].remove()
    fig.figure.savefig(
        f"figureS5.tif",
        dpi=100,
        bbox_inches="tight",
    )
    fig.figure.clf()
    return


def plot_lowest_rmsd():
    df = pd.read_parquet("native_and_docking_pose_check_qmean.gzip")
    df = df[df["dataset"] == "vina"]
    best_rmsd_df = df.sort_values("rmsd", ascending=True).drop_duplicates("bs_id")
    best_rmsd_df["Dataset"] = "Docked pose with lowest RMSD"

    best_score_df = df.sort_values("vina_score", ascending=True).drop_duplicates(
        "bs_id"
    )
    best_score_df["Dataset"] = "Docked pose with best score"
    df = pd.concat([best_rmsd_df, best_score_df])
    df = df.rename(columns={"rmsd": "RMSD [Å]"})

    sdf = pd.read_parquet("pub_roche_project_annotation_2024-12-03.gz")
    sdf = sdf.rename(columns={"bs_id": "identifier"})
    struc_qual = structure_quality.StructureQuality(sdf)
    df = df[df["bs_id"].isin(struc_qual.low_quality_bs_ids) == False]
    dist_plot = plot_distributions(
        df,
        hue="Dataset",
        outfile="figureS4.tif",
    )

    return


def main():
    # rmsd_df = rmsd_comparison()
    # rf_df = rf_comparison()
    # # rf_df.to_csv("rf_vina_vs_pdb.csv", index=False)
    # rf_df = pd.read_csv("rf_vina_vs_pdb.csv")
    # plot_relative_rf(rf_df)
    # plot_experimental_pose_recovery()
    plot_native_pose_rank()
    # props_vs_rmsd()
    # plot_lowest_rmsd()
    return


if __name__ == "__main__":
    main()
