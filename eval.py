import matplotlib.pyplot as plt
import pandas as pd
import pyarrow.feather as feather
import seaborn as sns
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-c", "--case", 
    help="name of the case identifier", 
    default="4bus", 
)
parser.add_argument(
    "-sf", "--suffix", 
    help="suffix for name of the figure to save", 
    default=None, 
)
args = parser.parse_args()

casename = args.case

if casename.find("4bus") != -1:
    name = "4 bus"
elif casename.find("13bus") != -1:
    name = "13 bus"
elif casename.find("123") != -1:
    name = "IEEE 123 bus"

real_act_file = f"outputs/{casename}/actual_real_PV.feather"
imag_act_file = f"outputs/{casename}/actual_imag_PV.feather"
real_est_file = f"outputs/{casename}/estimated_real_PV.feather"
imag_est_file = f"outputs/{casename}/estimated_imag_PV.feather"

df_act_p = feather.read_feather(real_act_file)
df_act_q = feather.read_feather(imag_act_file)
df_est_p = feather.read_feather(real_est_file)
df_est_q = feather.read_feather(imag_est_file)



def add_data_from_df(
        data_dict, 
        df, pq="p", type="actual",
        ):
    nodes = [n for n in df.columns.to_list() if n != "time"]
    for k in range(len(df)):
        for n in nodes:
            data_dict["nodes"].append(n)
            data_dict["value"].append(df.iloc[k][n])
            data_dict["pq"].append(pq)
            data_dict["type"].append(type)
    return data_dict


results = {
    "nodes" : [], 
    "value" : [], 
    "pq" : [],
    "type" : [],
}

results = add_data_from_df(results, df_act_p, pq="p", type="actual")
results = add_data_from_df(results, df_act_q, pq="q", type="actual")
results = add_data_from_df(results, df_est_p, pq="p", type="estimated")
results = add_data_from_df(results, df_est_q, pq="q", type="estimated")


df_results = pd.DataFrame(results)
df_pv_real = df_results.loc[df_results["pq"]=="p"]
df_pv_imag = df_results.loc[df_results["pq"]=="q"]

fig, axs = plt.subplots(2,1,figsize=(30,18))

sns.barplot(
    df_pv_real, x="nodes", y="value", hue="type",
    palette = ["royalblue", "crimson"],
    ax = axs[0]
)
axs[0].set_xlabel("Node in the network", fontsize=25)
axs[0].set_ylabel("PV System generation (p.u.)", fontsize=25)
axs[0].set_title("Real power estimation", fontsize=30)
axs[0].tick_params(labelsize=20, rotation=30)
# axs[0].set_ylim(0,0.4)
axs[0].legend(loc="upper right", fontsize=20, markerscale=2)


sns.barplot(
    df_pv_imag, x="nodes", y="value", hue="type",
    palette = ["royalblue", "crimson"],
    ax = axs[1]
)
axs[1].set_xlabel("Node in the network", fontsize=25)
axs[1].set_ylabel("PV System generation (p.u.)", fontsize=25)
axs[1].set_title("Reactive power estimation", fontsize=30)
axs[1].tick_params(labelsize=20, rotation=30)
# axs[1].set_ylim(0,0.4)
axs[1].legend(loc="upper right", fontsize=20, markerscale=2)

fig.suptitle(f"PV generation estimation for the {name} test feeder", fontsize=27)
if args.suffix == None:
    sfx = ""
else:
    sfx = f"_{args.suffix}"
fig.savefig(f"outputs/{casename}/pv_estimation{sfx}.png", bbox_inches='tight')