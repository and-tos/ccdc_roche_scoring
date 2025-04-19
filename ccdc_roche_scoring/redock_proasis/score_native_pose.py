from pathlib import Path
import pandas as pd
from ccdc import io

# from vina import Vina
import multiprocessing
from functools import partial
import subprocess as sp
import seaborn as sns
from matplotlib.font_manager import fontManager, FontProperties
from ccdc_roche_utils.plot_utils import roche_branding
from ccdc_roche_utils.plot_utils.seaborn_plots import plot_unity


def smallest_enclosing_cuboid(ccdc_ligand):
    center = ccdc_ligand.centre_of_geometry()
    points = [a.coordinates for a in ccdc_ligand.atoms]
    xdist = 2 * max([abs(point[0] - center[0]) for point in points]) + 8
    ydist = 2 * max([abs(point[1] - center[1]) for point in points]) + 8
    zdist = 2 * max([abs(point[2] - center[2]) for point in points]) + 8
    return list(center), [xdist, ydist, zdist]


def run_multiprocessing(nproc, iterable, function, **kwargs):
    parallel_func = partial(function, **kwargs)
    pool = multiprocessing.Pool(nproc)
    output = pool.map(parallel_func, iterable)
    return output


def rescore(dock_dir):
    # dock_dir = Path(
    #     "/home/tosstora/scratch/LoS/protonate3d_2024-07-12/redocking/300086"
    # )
    rec_pdbqt = dock_dir / "protein.pdbqt"
    vina_results = dock_dir / "vina_results.sdf"
    native_ligand_sdf = dock_dir / "ligand_native.sdf"
    native_ligand_pdbqt = dock_dir / "ligand_native.pdbqt"
    if vina_results.is_file() and rec_pdbqt.is_file():
        ccdc_ligand = io.MoleculeReader(str(native_ligand_sdf))[0]
        v = Vina()
        try:
            v.set_receptor(str(rec_pdbqt))
        except TypeError:
            return None, None, None
        prepare_ligand_args = f"-i {native_ligand_sdf} -o {native_ligand_pdbqt}".split()
        sp.check_output(
            [
                "ccdc_roche_scoring/bin/mk_prepare_ligand.py"
            ]
            + prepare_ligand_args
        )
        center, side_length = smallest_enclosing_cuboid(ccdc_ligand)
        v.set_ligand_from_file(str(native_ligand_pdbqt))
        v.compute_vina_maps(center, side_length)
        identifier = io.MoleculeReader(str(dock_dir / "protein.mol2"))[0].identifier
        vina_results = io.EntryReader(str(vina_results))[0]
        top_score = float(vina_results.attributes["vina_score"])
        try:
            native_score = v.score()[0]
        except:
            print(dock_dir)
            print(identifier)
            return None, None, None
        return (identifier, native_score, top_score)
    else:
        return None, None, None


def score_correlation_plot(df, x="vina_native_score", y="vina_top_score"):
    path = "/Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    plot = sns.scatterplot(
        df,
        x=x,
        y=y,
        color=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
        s=8,
    )
    plot_unity(df[x], df[y])
    plot.set(box_aspect=1, xlim=(-20, 0), ylim=(-20, 0))

    plot.figure.savefig("native_pose_rescore.png", dpi=600)
    plot.figure.clf()
    return


def score_kde_plot(df, x="delta_score"):
    path = "/Roche_fonts/Roche Sans/RocheSansLight-Medium.ttf"
    fontManager.addfont(path)
    prop = FontProperties(fname=path)
    sns.set(font=prop.get_name())
    plot = sns.histplot(
        df,
        x=x,
        binrange=(-5, 5),
        color=roche_branding.RocheColours().roche_colours_rgb_dict["roche_blue"],
    )

    plot.set(box_aspect=1)
    plot.set_xlabel("Δdocking_score")
    plot.figure.savefig("native_pose_delta_score.png", dpi=600)
    plot.figure.clf()
    return


def main():
    # dock_dirs = Path(
    #     "/home/tosstora/scratch/LoS/protonate3d_2024-12-03/redocking"
    # ).glob("*")
    # output = run_multiprocessing(32, dock_dirs, rescore)
    # # for dock_dir in dock_dirs:
    # #     rescore(dock_dir)
    # rescore_dict = {
    #     "identifier": [a[0] for a in output],
    #     "vina_native_score": [a[1] for a in output],
    #     "vina_top_score": [a[2] for a in output],
    # }
    # df = pd.DataFrame(rescore_dict).dropna()
    # df.to_csv("native_pose_rescore.csv", index=False)
    df = pd.read_csv("native_pose_rescore.csv")
    df["delta_score"] = df["vina_native_score"] - df["vina_top_score"]
    score_correlation_plot(df)
    score_kde_plot(df)
    return


if __name__ == "__main__":
    main()
