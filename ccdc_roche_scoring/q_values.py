from collections import defaultdict
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import datetime
from pathlib import Path
import pickle
from scipy.stats.mstats import gmean
import seaborn as sns
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib.colors import LogNorm
from ccdc import io, search, diagram

from ccdc_roche_utils.plot_utils import roche_branding


def estimate_maxima(kde):
    no_samples = 1000
    x = np.linspace(0, 180, no_samples)
    X, Y = np.meshgrid(x, x)
    xy = np.vstack((X.flatten(), Y.flatten()))
    probs = kde.evaluate(xy)
    maximum_prob = probs.max()

    return maximum_prob


def get_2d_histogram(df):
    hist = np.histogramdd(
        df[["TOR1", "TOR2"]].to_numpy(), bins=20, range=[(0, 180), (0, 180)]
    )
    fmax = hist[0].max()
    return hist, fmax


def get_unique_measurements(measurement_items, hit, return_objects=False):
    unique_measurement_items = []
    unique_measurement_objects = {}
    measurement_ids = set([i[0].split("_")[0] for i in measurement_items])
    for measurement_id in measurement_ids:
        measurements = [
            (abs(i[1]), i[0])
            for i in measurement_items
            if i[0].startswith(measurement_id)
        ]
        measurement, _measurement_id = min(measurements, key=lambda t: t[0])
        unique_measurement_items.append((measurement_id, measurement))
        if return_objects:
            unique_measurement_objects[measurement_id] = hit.measurement_objects(
                _measurement_id
            )
    if return_objects:
        return unique_measurement_items, unique_measurement_objects
    else:
        return unique_measurement_items


def get_1D_histogram(df):
    hist = np.histogram(df["TOR1"], bins=20, range=(0, 180))
    fmax = hist[0].max()
    return hist, fmax


def draw_molecule(ccdc_mol, searcher):
    hits = searcher.search(
        ccdc_mol.components, max_hit_structures=None, max_hits_per_structure=None
    )
    diagram_generator = diagram.DiagramGenerator()
    diagram_generator.settings.element_coloring = False
    diagram_generator.settings.shrink_symbols = False
    diagram_generator.settings.font_family = "Roche Sans Light"
    hit = hits[0]
    if hasattr(searcher.substructures[0], "smarts"):
        highlight_atoms = list(hit.measurement_atoms("TOR1")) + list(
            hit.measurement_atoms("TOR2")
        )
    else:
        highlight_atoms = []
        for key in hit.measurements.keys():
            highlight_atoms.extend(list(hit.measurement_objects(key)))
    highlight_atoms = list(set(highlight_atoms))
    img = diagram_generator.image(hit.molecule, highlight_atoms=highlight_atoms)
    return img


def get_unique_hits(hits):
    unique_hits = []
    unique_hit_atoms = []
    for hit in hits:
        match_atoms = set([a for a in hit.match_atoms() if a.atomic_symbol != "H"])
        if match_atoms not in unique_hit_atoms:
            unique_hit_atoms.append(match_atoms)
            unique_hits.append(hit)
    return unique_hits


class MoleculeDiagram(object):
    def __init__(self, molecule, searcher) -> None:
        self.molecule = molecule
        self.image = draw_molecule(molecule, searcher)
        pass


class QValuePlots(object):
    def __init__(self, q_value_home) -> None:
        self.q_value_home = q_value_home
        self.font_home = Path("/Roche_fonts")

    def make_density_plot(self, df, mol_diagram=None):
        # set font
        font_path = self.font_home / "Roche Sans/RocheSansLight-Medium.ttf"
        fontManager.addfont(font_path)
        prop = FontProperties(fname=font_path)
        sns.set(font=prop.get_name(), style="white")

        # generate plot
        geom1, geom2 = "TOR1", "TOR2"
        axlim = (-10, 190)
        ticks = range(0, 210, 30)
        _bins = np.linspace(0, 180, 20)
        plot = sns.jointplot()
        plot.ax_joint.hist2d(
            data=df,
            x=geom1,
            y=geom2,
            bins=[
                _bins,
                _bins,
            ],
            cmap="Greys",
            norm=LogNorm(),
        )
        # plot.ax_joint.collections[0].set_cmap("Greys")
        plot.ax_marg_x.hist(data=df, x=geom1, bins=_bins, color="gray", alpha=0.7)
        plot.ax_marg_y.hist(
            data=df,
            x=geom2,
            bins=_bins,
            color="gray",
            alpha=0.7,
            orientation="horizontal",
        )
        plot.ax_joint.set(
            xlim=axlim,
            ylim=axlim,
            xticks=ticks,
            yticks=ticks,
            box_aspect=1,
        )
        plot.ax_joint.set_xticklabels(labels=ticks, rotation=90)

        cbar_ax = plot.fig.add_axes([0.25, -0.10, 0.4, 0.05])  # x, y, width, height
        plt.colorbar(
            plot.ax_joint.collections[0],
            cax=cbar_ax,
            orientation="horizontal",
            label="CSD count",
        )
        # plot.figure.savefig(f"test_hist.png", dpi=300, bbox_inches="tight")
        # plot.figure.clf()

        if mol_diagram:
            newax = plot.figure.add_axes([1.05, 0.35, 0.3, 0.3], anchor="NE", zorder=-1)
            newax.imshow(mol_diagram.image)
            newax.set(title=mol_diagram.molecule.identifier)
            newax.axis("off")
        return plot

    def plot_density(self, df, substructure_name, img):
        outdir = Path(self.q_value_home / "plots")
        plot = self.make_density_plot(df, img)
        plot.figure.savefig(
            outdir / f"{substructure_name}_hexbin.png",
            dpi=300,
            bbox_inches="tight",
        )
        plot.figure.clf()
        return

    def plot_hist(self, df, substructure_name, img):
        outdir = self.q_value_home / "plots"

        # set font
        font_path = self.font_home / "Roche Sans/RocheSansLight-Medium.ttf"
        fontManager.addfont(font_path)
        prop = FontProperties(fname=font_path)
        sns.set(font=prop.get_name())

        # generate plot
        axlim = (-10, 200)
        ticks = range(0, 190, 30)

        plot = sns.histplot(
            df["TOR1"],
            bins=20,
            binrange=(0, 180),
            color=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
        )

        plot.set(
            xlim=axlim,
            xticks=ticks,
            box_aspect=1,
        )
        plot.set_xticklabels(labels=ticks, rotation=90)

        newax = plot.figure.add_axes([1.05, 0.3, 0.3, 0.3], anchor="NE", zorder=-1)
        newax.imshow(img)
        newax.axis("off")
        plot.figure.savefig(
            outdir / f"{substructure_name}_hist.png",
            dpi=300,
            bbox_inches="tight",
        )
        plot.figure.clf()
        return


class QValue(object):
    def __init__(
        self, q_value_home=Path("/q_values/")
    ) -> None:
        self.smarts_dict = {}
        self.q_value_home = q_value_home
        self.quest_2D_dir = Path(q_value_home / "quest_queries_2D_symmetric")
        self.quest_1D_dir = Path(q_value_home / "quest_queries_1D_symmetric")
        self.kde_dict_path = Path(q_value_home / "kde_dict.pkl")
        self.internal_csd = io.EntryReader(
            str(Path(io.csd_directory()) / "csd_roche_all.csdsql")
        )
        self.csd = io.EntryReader(str(Path(io.csd_directory()) / "as545be_ASER.sqlite"))
        self.substructures = self._get_substructure_dict()

    def _get_substructure_dict(self):
        substructure_dict = defaultdict(list)
        for substructure_name in self.smarts_dict:
            substructure_dict["substructure_name"].append(substructure_name)
            substructure_dict["smarts"].append(self.smarts_dict[substructure_name])
            substructure_dict["q_dim"].append(2)

        for quest_query in self.quest_1D_dir.glob("*.con"):
            substructure_name = quest_query.stem
            substructure_dict["substructure_name"].append(substructure_name)
            substructure_dict["smarts"].append("quest")
            substructure_dict["q_dim"].append(1)

        for quest_query in self.quest_2D_dir.glob("*.con"):
            substructure_name = quest_query.stem
            substructure_dict["substructure_name"].append(substructure_name)
            substructure_dict["smarts"].append("quest")
            substructure_dict["q_dim"].append(2)
        return substructure_dict

    def search_substructure_name(self, substructure_name, db="CSD"):
        index = self.substructures["substructure_name"].index(substructure_name)
        smarts = self.substructures["smarts"][index]
        q_dim = self.substructures["q_dim"][index]
        if smarts == "quest":
            if q_dim == 1:
                quest_files = [self.quest_1D_dir / f"{substructure_name}.con"]
            if q_dim == 2:
                quest_files = [self.quest_2D_dir / f"{substructure_name}.con"]
            smarts = []
        else:
            smarts = [smarts]
            quest_files = []
        hits, searcher = self.csd_search(db, smarts=smarts, quest_files=quest_files)
        return hits, searcher

    def csd_search(self, db, smarts=[], quest_files=[]):
        searcher = search.SubstructureSearch()
        for quest_file in quest_files:
            substructure = search.ConnserSubstructure(str(quest_file))
            searcher.add_substructure(substructure)

        for smarts_query in smarts:
            substructure = search.SMARTSSubstructure(smarts_query)
            searcher.add_substructure(substructure)
            searcher.add_torsion_angle_measurement(
                "TOR1", (0, 1), (0, 3), (0, 4), (0, 5)
            )
            searcher.add_torsion_angle_measurement(
                "TOR2", (0, 0), (0, 1), (0, 3), (0, 4)
            )
        if db == "CSD":
            db = [self.internal_csd, self.csd]
            searcher.settings.only_organic = True
            searcher.settings.not_polymeric = True
            searcher.settings.no_errors = True
            searcher.settings.no_disorder = True
            searcher.settings.max_r_factor = 10
        hits = searcher.search(db, max_hit_structures=None, max_hits_per_structure=None)
        return hits, searcher

    def return_hits_df(self, hits, ccdc_mol=None):
        measurement_dict = defaultdict(list)
        for cnt, hit in enumerate(hits):
            measurement_items = hit.measurements.items()
            if ccdc_mol:
                (
                    unique_measurement_items,
                    unique_measurement_objects,
                ) = get_unique_measurements(measurement_items, hit, return_objects=True)
            else:
                unique_measurement_items = get_unique_measurements(
                    measurement_items, hit, return_objects=False
                )

            for key, value in unique_measurement_items:
                if value == 180.0:
                    value = 179.9  # include in histogram
                measurement_dict[key].append(abs(value))
                if ccdc_mol:
                    atom1_index = ccdc_mol.atom(
                        unique_measurement_objects[key][1].label
                    ).index
                    atom2_index = ccdc_mol.atom(
                        unique_measurement_objects[key][2].label
                    ).index
                    measurement_dict[f"{key}_atom1"].append(atom1_index)
                    measurement_dict[f"{key}_atom2"].append(atom2_index)
            measurement_dict["identifier"].append(hit.identifier)
        measurement_df = pd.DataFrame(measurement_dict).drop_duplicates()
        return measurement_df

    def evaluate_2d_hist(self, df, hist, fmax):
        torsions_1 = df["TOR1"]
        bin1 = np.digitize(torsions_1, bins=hist[1][0]) - 1
        torsions_2 = df["TOR2"]
        bin2 = np.digitize(torsions_2, bins=hist[1][1]) - 1
        bins = zip(bin1, bin2)
        q2 = []
        for bin in bins:
            f = hist[0][bin]
            q2.append(f / fmax)
        return q2

    def evaluate_1d_hist(self, df, hist, fmax):
        torsions = df["TOR1"]
        bin = np.digitize(torsions, bins=hist[1]) - 1
        f = hist[0][bin]
        q = f / fmax
        return q

    def get_q_values(self, ccdc_mol):
        _ccdc_mol = ccdc_mol.copy()
        _ccdc_mol.add_hydrogens()
        _ccdc_mol.standardise_aromatic_bonds()
        if not self.kde_dict_path.is_file():
            Exception("KDE dictionary does not exit.")
        with open(self.kde_dict_path, "rb") as file:
            kde_dict = pickle.load(file)

        measurement_dfs = []
        measurement_df = pd.DataFrame()
        for substructure_name in kde_dict:
            substructure_dict = kde_dict[substructure_name]
            smarts = substructure_dict["smarts"]
            kde = substructure_dict["kde"]
            fmax = substructure_dict["fmax"]
            q_dim = substructure_dict["q_dim"]
            if smarts == "quest":
                if q_dim == 2:
                    quest_file = self.quest_2D_dir / f"{substructure_name}.con"
                elif q_dim == 1:
                    quest_file = self.quest_1D_dir / f"{substructure_name}.con"
                hits, searcher = self.csd_search(_ccdc_mol, quest_files=[quest_file])
                measurement_df = self.return_hits_df(hits, ccdc_mol)
            else:
                hits, searcher = self.csd_search(_ccdc_mol, smarts=[smarts])
                measurement_df = self.return_hits_df(hits)
            measurement_df["smarts"] = smarts
            measurement_df["q_dim"] = q_dim
            measurement_df["substructure_name"] = substructure_name

            if measurement_df.shape[0] > 0:
                if q_dim == 2:
                    measurement_df["q2"] = self.evaluate_2d_hist(
                        measurement_df, kde, fmax
                    )
                if q_dim == 1:
                    measurement_df["q1"] = self.evaluate_1d_hist(
                        measurement_df, kde, fmax
                    )
                measurement_dfs.append(measurement_df)

        if measurement_dfs:
            measurement_df = pd.concat(measurement_dfs, ignore_index=True)
        return measurement_df.drop_duplicates().reset_index(drop=True)

    def get_q2_gmean(self, ccdc_mol):
        df = self.get_q_values(ccdc_mol)
        q_values = []
        for column in df.columns:
            if column.startswith("q1") or column.startswith("q2"):
                q_values.extend(df[df[column].isna() == False][column].tolist())
        if q_values:
            q2_gmean = gmean(q_values)
        else:
            q2_gmean = 1.1
        return q2_gmean

    def update_q2_kde_dict(self, kde_dict, hits, substructure_name, searcher):
        csd_df = self.return_hits_df(hits)
        smallest_hit = min(hits, key=lambda x: len(x.molecule.heavy_atoms))
        mol_diagram = MoleculeDiagram(smallest_hit.molecule, searcher)
        hist, fmax = get_2d_histogram(csd_df)
        substructure_dict = {}
        substructure_dict["smarts"] = "quest"
        substructure_dict["kde"] = hist
        substructure_dict["fmax"] = fmax
        substructure_dict["CSD_hits"] = csd_df.shape[0]
        substructure_dict["q_dim"] = 2
        kde_dict[substructure_name] = substructure_dict
        QValuePlots(self.q_value_home).plot_density(
            csd_df, substructure_name, mol_diagram
        )
        return

    def update_q1_kde_dict(self, kde_dict, hits, substructure_name, searcher):
        csd_df = self.return_hits_df(hits)
        smallest_hit = min(hits, key=lambda x: len(x.molecule.heavy_atoms))
        hist, fmax = get_1D_histogram(csd_df)
        img = draw_molecule(smallest_hit.molecule, searcher)
        substructure_dict = {}
        substructure_dict["smarts"] = "quest"
        substructure_dict["kde"] = hist
        substructure_dict["fmax"] = fmax
        substructure_dict["CSD_hits"] = csd_df.shape[0]
        substructure_dict["q_dim"] = 1
        kde_dict[substructure_name] = substructure_dict
        QValuePlots(self.q_value_home).plot_hist(csd_df, substructure_name, img)
        return

    def write_kde_dict(self, rewrite_dict=False):
        kde_dict_time = datetime.datetime.fromtimestamp(
            self.kde_dict_path.stat().st_mtime
        )
        if rewrite_dict == False and self.kde_dict_path.is_file():
            with open(self.kde_dict_path, "rb") as file:
                kde_dict = pickle.load(file)
        else:
            kde_dict = {}
        substructure_names = []
        for quest_query in self.quest_2D_dir.glob("*.con"):
            substructure_name = quest_query.stem
            substructure_names.append(substructure_name)

            mtime = datetime.datetime.fromtimestamp(quest_query.stat().st_mtime)
            if not rewrite_dict and kde_dict_time > mtime:
                continue

            hits, searcher = self.csd_search(
                [self.internal_csd, self.csd],
                quest_files=[quest_query],
            )
            print(substructure_name, "hits: ", len(hits))

            if len(hits) >= 60:
                self.update_q2_kde_dict(kde_dict, hits, substructure_name, searcher)

        for quest_query in self.quest_1D_dir.glob("*.con"):
            substructure_name = quest_query.stem
            substructure_names.append(substructure_name)

            mtime = datetime.datetime.fromtimestamp(quest_query.stat().st_mtime)
            if not rewrite_dict and kde_dict_time > mtime:
                continue

            self.symmetric = substructure_name.endswith("sym")

            hits, searcher = self.csd_search(
                [self.internal_csd, self.csd],
                quest_files=[quest_query],
            )

            print(substructure_name, "hits: ", len(hits))
            if len(hits) >= 50:
                self.update_q1_kde_dict(kde_dict, hits, substructure_name, searcher)

        # remove queries that were deleted
        del_substructure_names = [
            sn for sn in kde_dict.keys() if sn not in substructure_names
        ]
        for sn in del_substructure_names:
            del kde_dict[sn]
        # for substructure_name in self.smarts_dict.keys():
        #     smarts = self.smarts_dict[substructure_name]

        #     if substructure_name in kde_dict["substructure_name"]:
        #         continue

        #     hits, searcher = self.csd_search(
        #         [self.internal_csd, self.csd], smarts=[smarts]
        #     )
        #     if len(hits) >= 60:
        #         print(substructure_name, "hits: ", len(hits))
        #         self.update_q2_kde_dict(
        #             kde_dict, hits, substructure_name, searcher, smarts=smarts
        #         )

        with open(self.kde_dict_path, "wb") as file:
            pickle.dump(kde_dict, file)
        return


def main():
    q_values = QValue()
    q_values.write_kde_dict(rewrite_dict=False)
    ccdc_mol = io.MoleculeReader("/tmp/test.sdf")[0]
    q_gmean = q_values.get_q2_gmean(ccdc_mol)
    return


if __name__ == "__main__":
    main()
