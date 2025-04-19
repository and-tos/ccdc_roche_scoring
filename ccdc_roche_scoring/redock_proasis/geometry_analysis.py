from collections import defaultdict
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.font_manager import fontManager, FontProperties
from pathlib import Path
from scipy.stats import gaussian_kde
from sklearn.preprocessing import MinMaxScaler
import numpy as np
from ccdc import io, entry, descriptors
from rf_statistics.utils import run_multiprocessing
from ccdc_roche_utils.plot_utils import roche_branding
from ccdc_roche_utils.plot_utils.seaborn_plots import plot_unity
from ccdc_roche_scoring import q_values


def estimate_maxima(kde):
    no_samples = 1000
    x = np.linspace(-180, 180, no_samples)
    X, Y = np.meshgrid(x, x)
    xy = np.vstack((X.flatten(), Y.flatten()))
    probs = kde.evaluate(xy)
    maximum_prob = probs.max()

    return maximum_prob


def write_gold_docking_csdsql():
    # generate csdsql
    rdr = io.EntryReader("top_poses.sdf")
    with io.EntryWriter("top_poses_.sdf") as w:
        for cnt, e in enumerate(rdr):
            m = e.molecule
            m.standardise_aromatic_bonds()
            new_e = entry.Entry.from_molecule(m)
            new_e.attributes = e.attributes
            new_e.identifier = e.attributes["Gold.Id.Protein"].split("|")[0]
            w.write(new_e)

    return


# def write_reference_csdsql():
#     # generate csdsql
#     reference_files = Path(
#         "rf_scoring/validation/redock_proasis/docking"
#     ).glob("*/ligand.sdf")
#     with io.EntryWriter("reference_ligands.csdsql") as w:
#         for reference_file in reference_files:
#             e = io.EntryReader(str(reference_file))[0]
#             e.identifier = reference_file.parent.name
#             m = e.molecule
#             m.standardise_aromatic_bonds()
#             new_e = entry.Entry.from_molecule(m)
#             w.write(new_e)

#     return


def get_2d_q_value(df):
    kde_input = df[df["Database"] == "CSD"][["TOR1", "TOR2"]].T
    kde = gaussian_kde(kde_input, bw_method=None, weights=None)
    qmax = estimate_maxima(kde)

    # docking
    docking_kde_input = df[df["Database"] == "docking"][["TOR1", "TOR2"]].T
    q = kde.evaluate(docking_kde_input) / qmax
    df.loc[df["Database"] == "docking", "q_value"] = q

    # reference
    docking_kde_input = df[df["Database"] == "reference"][["TOR1", "TOR2"]].T
    q = kde.evaluate(docking_kde_input) / qmax
    df.loc[df["Database"] == "reference", "q_value"] = q

    return


def csd_vs_docking_1D_plot(df, substructure_name):
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    geom1 = "TOR1"
    axlim = (-10, 190)
    ticks = range(0, 210, 30)

    df = df[df["Database"].isin(["CSD", "docking"])]

    # plot = sns.barplot(
    #     df,
    #     x="TOR1",
    #     y="q",
    #     hue="Database",
    #     palette=[
    #         roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
    #         roche_branding.RocheColours().roche_colours_rgb_dict["dark_orange"],
    #         roche_branding.RocheColours().roche_colours_rgb_dict["purple"],
    #     ],
    #     alpha=0.4,
    # )
    n_csd = df[df["Database"] == "CSD"].shape[0]
    n_docking = df[df["Database"] == "docking"].shape[0]
    df.loc[df["Database"] == "CSD", "Database"] = df["Database"].str.replace(
        "CSD", f"CSD, N={n_csd}"
    )
    df.loc[df["Database"] == "docking", "Database"] = df["Database"].str.replace(
        "docking", f"docking, N={n_docking}"
    )

    plot = sns.histplot(
        df,
        x="TOR1",
        bins=np.linspace(0, 180, 20),
        hue="Database",
        palette=[
            roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
            roche_branding.RocheColours().roche_colours_rgb_dict["dark_orange"],
            roche_branding.RocheColours().roche_colours_rgb_dict["purple"],
        ],
        alpha=0.4,
        multiple="layer",
        stat="probability",
        common_norm=False,
    )
    plot.set(xlim=axlim, xticks=ticks, box_aspect=1, xlabel="TOR1 [°]")
    plot.figure.savefig(
        f"{substructure_name}_{geom1}_hist.png", dpi=300, bbox_inches="tight"
    )
    plot.figure.clf()
    return


def csd_vs_docking_2D_plot(df, substructure_name):
    # set font
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())

    geom1, geom2 = "TOR1", "TOR2"
    axlim = (-10, 190)
    ticks = range(0, 210, 30)

    rmsd_palette = "coolwarm"
    marker_size = 15
    plot = q_values.QValuePlots(".").make_density_plot(df[df["Database"] == "CSD"])

    scatter_plot = sns.scatterplot(
        data=df[df["Database"] == "docking"],
        x=geom1,
        y=geom2,
        edgecolors="black",
        hue="RMSD",
        hue_norm=(0, 2.5),
        s=marker_size,
        palette=rmsd_palette,
        ax=plot.ax_joint,
    )
    plot.ax_joint.legend(
        labels=["docking"],
        loc="lower center",
        bbox_to_anchor=(0.95, -0.25),
        ncol=1,
    )

    plot.ax_joint.set(
        xlim=axlim,
        ylim=axlim,
        xticks=ticks,
        yticks=ticks,
        box_aspect=1,
        xlabel="TOR1 [°]",
        ylabel="TOR2 [°]",
    )
    plot.ax_joint.set_xticklabels(labels=ticks, rotation=90)

    norm = plt.Normalize(0, 2.5)
    sm = plt.cm.ScalarMappable(cmap=rmsd_palette, norm=norm)
    sm.set_array([])
    cbar_ax = plot.fig.add_axes([0.9, 0.1, 0.1, 0.8])  # x, y, width, height
    cbar = plt.colorbar(sm, ax=cbar_ax)
    cbar.set_label("RMSD [Å]")
    cbar_ax.grid("False")
    cbar_ax.axis("off")

    plot.figure.savefig(
        f"{substructure_name}_{geom1}_vs_{geom2}_kde.png", dpi=300, bbox_inches="tight"
    )
    plot.figure.clf()

    # Q-value distribution
    # plot = sns.histplot(
    #     data=df[df["Database"] != "CSD"],
    #     x="q",
    #     hue="Database",
    #     element="step",
    #     palette=[
    #         roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
    #         roche_branding.RocheColours().roche_colours_rgb_dict["dark_orange"],
    #     ],
    # )
    # plot.set(xlim=[0, 1])
    # plot.figure.savefig(
    #     f"{substructure_name}_{geom1}_vs_{geom2}_q_values.png",
    #     dpi=300,
    #     bbox_inches="tight",
    # )
    # plot.figure.clf()

    # Q-values vs RMSD
    # plot = sns.scatterplot(
    #     data=df[df["Database"] != "CSD"],
    #     x="q",
    #     y="RMSD",
    #     hue="Database",
    #     palette=[
    #         roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
    #         roche_branding.RocheColours().roche_colours_rgb_dict["dark_orange"],
    #     ],
    # )
    # plot.set(xlim=[0, 1])
    # plot.figure.savefig(
    #     f"{substructure_name}_{geom1}_vs_{geom2}_q_values_rmsd.png",
    #     dpi=300,
    #     bbox_inches="tight",
    # )
    # plot.figure.clf()
    return


def get_docking_df(searcher, substructure_name, db, sdf_reader):
    hits = searcher.search(db, max_hit_structures=1000, max_hits_per_structure=3)

    measurement_dict = defaultdict(list)
    for cnt, hit in enumerate(hits):
        for key, value in hit.measurements.items():
            measurement_dict[key].append(value)
        identifier = hit.identifier.split("|")[0]
        measurement_dict["identifier"].append(identifier)
        measurement_dict["RMSD"].append(sdf_reader.entry(identifier).attributes["RMSD"])
        if cnt == 0:
            filename = f"{substructure_name}_docking.sdf"
            with io.EntryWriter(filename) as w:
                e = entry.Entry.from_molecule(hit.molecule)
                for key, value in hit.measurements.items():
                    e.attributes[key] = value
                w.write(e)
    measurement_df = pd.DataFrame(measurement_dict)
    return measurement_df


def get_reference_df(searcher, substructure_name, db):
    hits = searcher.search(db, max_hit_structures=1000, max_hits_per_structure=3)

    measurement_dict = defaultdict(list)
    for cnt, hit in enumerate(hits):
        for key, value in hit.measurements.items():
            measurement_dict[key].append(value)
        identifier = hit.identifier
        measurement_dict["identifier"].append(identifier)

        if cnt == 0:
            filename = f"{substructure_name}_reference.sdf"
            with io.EntryWriter(filename) as w:
                e = entry.Entry.from_molecule(hit.molecule)
                for key, value in hit.measurements.items():
                    e.attributes[key] = value
                w.write(e)
    measurement_df = pd.DataFrame(measurement_dict)
    return measurement_df


def get_csd_df(searcher, substructure_name, db):
    hits = searcher.search(db, max_hit_structures=1000, max_hits_per_structure=3)

    measurement_dict = defaultdict(list)
    for cnt, hit in enumerate(hits):
        for key, value in hit.measurements.items():
            measurement_dict[key].append(value)
        measurement_dict["identifier"].append(hit.identifier)

        if cnt == 0:
            filename = f"{substructure_name}_csd.sdf"
            with io.EntryWriter(filename) as w:
                e = entry.Entry.from_molecule(hit.molecule)
                for key, value in hit.measurements.items():
                    e.attributes[key] = value
                w.write(e)
    measurement_df = pd.DataFrame(measurement_dict)
    return measurement_df


def calculate_rmsd(identifiers, docking_rdr, refernce_rdr):
    rmsds = []
    for identifier in identifiers:
        rmsd = 0
        try:
            docked_mol = docking_rdr.entry(identifier).molecule
            ref_mol = refernce_rdr.entry(identifier).molecule
            rmsd = descriptors.MolecularDescriptors.rmsd(docked_mol, ref_mol)
        except:
            pass
        rmsds.append(rmsd)
    return rmsds


def calculate_histogram_distance(hist1, hist2, distance_metric="euclidean"):
    """
    Calculate the distance between two 2D histograms.

    :param hist1: First histogram, a 2D NumPy array
    :param hist2: Second histogram, a 2D NumPy array
    :param distance_metric: The type of distance metric to use ('euclidean', 'manhattan')
    :return: The distance between the two histograms
    """
    if len(hist1) != len(hist2):
        raise ValueError("The histograms must have the same shape.")

    scaler = MinMaxScaler()
    scaler.fit_transform(hist1)
    hist1 = scaler.transform(hist1)

    scaler = MinMaxScaler()
    scaler.fit(hist2)
    hist2 = scaler.transform(hist2)

    if distance_metric == "euclidean":
        distance = np.sqrt(np.sum((hist1 - hist2) ** 2))
    elif distance_metric == "manhattan":
        distance = np.sum(np.abs(hist1 - hist2))
    else:
        raise ValueError(f"Unsupported distance metric: {distance_metric}")

    return distance


def euclidean_distance_correlation_plot(df):
    path = "Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    bins = [0, 50, 500, 1000, 200000]
    classification = {
        "(0, 50]": "<50",
        "(50, 500]": "50 to 500",
        "(500, 1000]": "501 to 1000",
        "(1000, 200000]": ">1000",
    }
    df["bins"] = pd.cut(df["docking_hits"], bins).apply(str)
    df["# docking hits"] = df["bins"].apply(lambda id: classification[id])
    plot = sns.scatterplot(
        df,
        x="csd_experimental_distance",
        y="csd_docking_distance",
        hue="dimensions",
        palette=[
            roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
            roche_branding.RocheColours().roche_colours_rgb_dict["dark_orange"],
            roche_branding.RocheColours().roche_colours_rgb_dict["purple"],
        ],
        edgecolor="black",
        linewidth=0.5,
        size="# docking hits",
        sizes={
            "<50": 5,
            "50 to 500": 50,
            "501 to 1000": 100,
            ">1000": 200,
        },
        legend="brief",
    )
    plot.set(xlabel="ε(CSD, experimental)", ylabel="ε(CSD, docking)")

    plot_unity(df["csd_experimental_distance"], df["csd_docking_distance"])
    plot.set(
        box_aspect=1,
    )

    plot.figure.savefig("euclidean_distance_vina_vs_pdb.png", dpi=600)
    plot.figure.clf()
    return


def _update_euclidean_distance_dict(
    q_dim, euclidean_distances, csd_df, proasis_df, docking_df, substructure_name
):
    if q_dim == 2:
        proasis_hist, fmax = q_values.get_2d_histogram(proasis_df)
        docking_hist, fmax = q_values.get_2d_histogram(docking_df)
        csd_hist, fmax = q_values.get_2d_histogram(csd_df)
        proasis_hist = proasis_hist[0]
        docking_hist = docking_hist[0]
        csd_hist = csd_hist[0]
    if q_dim == 1:
        proasis_hist, fmax = q_values.get_1D_histogram(proasis_df)
        docking_hist, fmax = q_values.get_1D_histogram(docking_df)
        csd_hist, fmax = q_values.get_1D_histogram(csd_df)

        proasis_hist = proasis_hist[0].reshape(-1, 1)
        docking_hist = docking_hist[0].reshape(-1, 1)
        csd_hist = csd_hist[0].reshape(-1, 1)

    csd_experimental_distance = calculate_histogram_distance(
        csd_hist, proasis_hist, distance_metric="euclidean"
    )
    csd_docking_distance = calculate_histogram_distance(
        csd_hist, docking_hist, distance_metric="euclidean"
    )
    proasis_docking_distance = calculate_histogram_distance(
        proasis_hist, docking_hist, distance_metric="euclidean"
    )
    euclidean_distances["substructure_name"].append(substructure_name)
    euclidean_distances["csd_experimental_distance"].append(csd_experimental_distance)
    euclidean_distances["csd_docking_distance"].append(csd_docking_distance)
    euclidean_distances["proasis_docking_distance"].append(proasis_docking_distance)
    euclidean_distances["csd_hits"].append(csd_df.shape[0])
    euclidean_distances["docking_hits"].append(docking_df.shape[0])
    euclidean_distances["proasis_hits"].append(proasis_df.shape[0])
    euclidean_distances["dimensions"].append(q_dim)
    return euclidean_distances


def calculate_q_value(cnt):
    q_evaluate = q_values.QValue()
    substructure_dict = q_evaluate.substructures
    substructure_name = substructure_dict["substructure_name"][cnt]
    print(substructure_name)

    # if substructure_name != "phenylalanine":
    #     return

    docking_rdr = io.EntryReader(
        [
            "../vina_top_docked_ligands_pub.csdsql",
            "../vina_top_docked_ligands_roche.csdsql",
        ]
    )
    reference_rdr = io.EntryReader(
        ["../reference_ligands_pub.csdsql", "../reference_ligands_roche.csdsql"]
    )
    euclidean_distances = defaultdict(list)
    # try:
    smarts = substructure_dict["smarts"][cnt]
    df_file = Path(f"{substructure_name}.csv")
    if df_file.is_file():
        df = pd.read_csv(df_file)
    else:
        docking_hits, searcher = q_evaluate.search_substructure_name(
            substructure_name, docking_rdr
        )
        if docking_hits:
            csd_hits, searcher = q_evaluate.search_substructure_name(substructure_name)
            csd_df = q_evaluate.return_hits_df(csd_hits)
            q_dim = len(
                set([tor.split("_")[0] for tor in docking_hits[0].measurements.keys()])
            )
            if csd_hits and len(csd_hits) > 50:
                if q_dim == 2:
                    csd_hist, fmax = q_values.get_2d_histogram(csd_df)
                elif q_dim == 1:
                    csd_hist, fmax = q_values.get_1D_histogram(csd_df)

                docking_df = q_evaluate.return_hits_df(docking_hits)
                docking_df["RMSD"] = calculate_rmsd(
                    docking_df["identifier"], docking_rdr, reference_rdr
                )

                proasis_hits, searcher = q_evaluate.search_substructure_name(
                    substructure_name, reference_rdr
                )
                proasis_df = q_evaluate.return_hits_df(proasis_hits)
                if csd_df.shape[0] > 0:
                    csd_df["Database"] = "CSD"
                    docking_df["Database"] = "docking"
                    proasis_df["Database"] = "reference"
                    df = pd.concat(
                        [csd_df, proasis_df, docking_df],
                        ignore_index=True,
                    )
                    if q_dim == 2:
                        df["q"] = q_evaluate.evaluate_2d_hist(df, csd_hist, fmax)
                    elif q_dim == 1:
                        df["q"] = q_evaluate.evaluate_1d_hist(df, csd_hist, fmax)
                    df.to_csv(df_file, index=False)

                if (
                    "CSD" in df["Database"].values
                    and "docking" in df["Database"].values
                    and "reference" in df["Database"].values
                ):
                    if q_dim == 2:
                        csd_vs_docking_2D_plot(df, substructure_name)
                    if q_dim == 1:
                        csd_vs_docking_1D_plot(df, substructure_name)
                    euclidean_distances = _update_euclidean_distance_dict(
                        q_dim,
                        euclidean_distances,
                        csd_df,
                        proasis_df,
                        docking_df,
                        substructure_name,
                    )

    if df_file.is_file():
        df = pd.read_csv(df_file)
        q_dim = len([c for c in df.columns if c.startswith("TOR")])
        if q_dim == 2:
            get_2d_q_value(df)
            csd_vs_docking_2D_plot(df, substructure_name)

            ref_df = df[df["Database"] == "reference"]
            low_q_entry = ref_df[ref_df["q_value"] == ref_df["q_value"].min()][
                "identifier"
            ].values[0]
            low_q_value = ref_df["q_value"].min()
            e = reference_rdr.entry(low_q_entry)
            with io.EntryWriter(
                f"{substructure_name}_{low_q_entry}_q_{low_q_value:.2f}.sdf"
            ) as w:
                w.write(e)
        elif q_dim == 1:
            csd_vs_docking_1D_plot(df, substructure_name)
        csd_df = df[df["Database"] == "CSD"]
        proasis_df = df[df["Database"] == "reference"]
        docking_df = df[df["Database"] == "docking"]
        euclidean_distances = _update_euclidean_distance_dict(
            q_dim,
            euclidean_distances,
            csd_df,
            proasis_df,
            docking_df,
            substructure_name,
        )
    # except:
    #     print("Error: ", substructure_name)
    #     return
    return euclidean_distances


def main():
    q_evaluate = q_values.QValue()
    substructure_dict = q_evaluate.substructures
    # for cnt in range(len(substructure_dict["substructure_name"])):
    #     output = calculate_q_value(cnt)
    euclidean_distances = run_multiprocessing(
        24, range(len(substructure_dict["substructure_name"])), calculate_q_value
    )
    euclidean_distances = pd.concat(
        [
            pd.DataFrame(euclidean_distance)
            for euclidean_distance in euclidean_distances
            if euclidean_distances == euclidean_distances
        ],
        ignore_index=True,
    )
    euclidean_distances.to_csv("euclidean_distances.csv", index=False, sep="\t")
    euclidean_distances = pd.read_csv("euclidean_distances.csv", sep="\t")
    euclidean_distance_correlation_plot(euclidean_distances)

    return


if __name__ == "__main__":
    main()
