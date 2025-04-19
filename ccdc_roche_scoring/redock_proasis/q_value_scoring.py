from pathlib import Path
import pandas as pd
from scipy.stats.mstats import gmean
from collections import defaultdict
from rf_statistics.utils import run_multiprocessing


def _return_q_gmean(df_file):
    df = pd.read_csv(df_file, sep="\t")

    if df.shape[0] > 0:
        q_gmean_dict = defaultdict(list)
        if "mol_cnt" in df.columns:
            for mol_cnt, gdf in df.groupby("mol_cnt"):
                _q_values = []
                for c in gdf.columns:
                    if c.startswith("q"):
                        _q_values.extend(gdf[gdf[c].isna() == False][c].values)
                q_gmean = gmean(_q_values)
                q_gmean_dict["identifier"].append(df_file.parent.name)
                q_gmean_dict["q_gmean"].append(q_gmean)
                q_gmean_dict["mol_cnt"].append(mol_cnt)
            return pd.DataFrame(q_gmean_dict)
        else:
            df_file.unlink()
    else:
        return df


def collect_results(dataset):
    df_files = Path(".").glob(f"*/q_values_{dataset}.csv")
    q_mean_dfs = []
    # for df_file in df_files:
    #     q_mean_df = _return_q_gmean(df_file)
    #     q_mean_dfs.append(pd.DataFrame(q_mean_df))
    q_mean_dfs = run_multiprocessing(24, df_files, _return_q_gmean)
    q_mean_df = pd.concat(q_mean_dfs, ignore_index=True)
    q_mean_df.to_csv(f"q_mean_{dataset}.csv", index=False)
    return q_mean_df


def main():
    rf_df = pd.read_parquet(
        "/LoS/molecular_recognition_vina_2024-12-03/native_and_docking_pose_check.gzip"
    )
    native_df = collect_results("native")
    vina_df = collect_results("vina")

    rf_df_vina = rf_df[rf_df["dataset"] == "vina"].sort_values(
        ["bs_id", "vina_score"], ascending=[True, True]
    )
    rf_df_vina["mol_cnt"] = rf_df_vina.groupby("bs_id").cumcount()
    vina_df = vina_df.sort_values(["identifier", "mol_cnt"], ascending=[True, True])
    rf_df_vina = rf_df_vina.join(
        vina_df.set_index(["identifier", "mol_cnt"]), on=["bs_id", "mol_cnt"]
    )

    rf_df_native = rf_df[rf_df["dataset"] == "native"].sort_values(
        ["bs_id"], ascending=True
    )

    native_df = native_df.sort_values("identifier", ascending=True)
    rf_df_native = rf_df_native.join(native_df.set_index("identifier"), on="bs_id")

    combined = pd.concat([rf_df_vina, rf_df_native], ignore_index=True)
    combined.to_parquet(
        "/LoS/molecular_recognition_vina_2024-12-03/native_and_docking_pose_check_qmean.gzip"
    )

    return


if __name__ == "__main__":
    main()
