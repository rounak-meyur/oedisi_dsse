import sys, os
import matplotlib.pyplot as plt
import matplotlib
import pyarrow.feather as feather
from oedisi.types.data_types import Topology
from datetime import datetime
import numpy as np
import json
import pandas as pd


def get_time(x):
    return datetime.strptime(
        x, '%Y-%m-%d %H:%M:%S'
        ).time().strftime("%H:%M")

def get_vmag_from_complex(realV, imagV):
    voltage_real = feather.read_feather(realV)
    voltage_imag = feather.read_feather(imagV)
    df_voltages = np.abs(voltage_real.drop("time", axis=1) + 1j * voltage_imag.drop("time", axis=1))
    df_voltages["time"] = voltage_real["time"].apply(get_time)
    return df_voltages.set_index("time")

def get_varg_from_complex(realV, imagV):
    voltage_real = feather.read_feather(realV)
    voltage_imag = feather.read_feather(imagV)
    df_angles = np.angle(voltage_real.drop("time", axis=1) + 1j * voltage_imag.drop("time", axis=1), deg=True)
    df_angles["time"] = voltage_real["time"].apply(get_time)
    return df_angles.set_index("time")

def get_base_voltages(topology_file):
    with open(topology_file) as f:
        topology = Topology.parse_obj(json.load(f))
        base_voltage_df = pd.DataFrame(
            {
                "id": topology.base_voltage_magnitudes.ids,
                "value": topology.base_voltage_magnitudes.values,
            }
        )
        base_voltage_df.set_index("id", inplace=True)
        base_voltages = base_voltage_df["value"]
    return base_voltages

def compare_vmag(
        df_true, df_est, base_voltages, bus, 
        to_file = None, show=False, do_return=False,
        **kwargs):
    figsize = kwargs.get('figsize', (10, 10))
    constrained_layout = kwargs.get('constrained_layout', False)
    label_fontsize = kwargs.get('fontsize', 25)
    legend_fontsize = label_fontsize + 2
    ticklabel_fontsize = label_fontsize - 2
    title_fontsize = label_fontsize + 10
    suptitle_sfx = kwargs.get('suptitle_sfx',None)

    true_voltages = df_true[bus] / base_voltages[bus]
    est_voltages = df_est[bus] / base_voltages[bus]

    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=constrained_layout)
    xaxis = np.arange(df_est.shape[0])
    y_true = [true_voltages.iloc[i] for i in xaxis]
    y_est = [est_voltages.iloc[i] for i in xaxis]
    ax.plot(xaxis, y_true, color='red', marker='o', ls='dashed', lw=2.0, label="True")
    ax.plot(xaxis, y_est, color='blue', marker='o', ls='dashed', lw=2.0, label="Estimated")
    suptitle = f"Voltage magnitude comparison for bus {bus}"

    ax.set_xlabel("Time", fontsize=label_fontsize)
    ax.set_ylabel(f"Voltage Magnitudes (p.u.)", fontsize=label_fontsize)
    ax.tick_params(axis="x", labelsize=ticklabel_fontsize)
    ax.tick_params(axis="y", labelsize=ticklabel_fontsize)
    ax.legend(fontsize=legend_fontsize, markerscale=2)
    fig.suptitle(suptitle, fontsize=title_fontsize)

    if to_file:
        fig.savefig(to_file, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
    
    if do_return:
        return fig
    pass

if __name__ == "__main__":
    case = "test"
    realVfile = os.path.join(f"outputs/{case}/voltage_real.feather")
    imagVfile = os.path.join(f"outputs/{case}/voltage_imag.feather")
    Vmagfile = os.path.join(f"outputs/{case}/voltage_mag.feather")
    Vangfile = os.path.join(f"outputs/{case}/voltage_angle.feather")
    topofile = os.path.join(f"outputs/{case}/topology.json")

    base_voltages = get_base_voltages(topology_file=topofile)
    df_vmag_true = get_vmag_from_complex(realV=realVfile, imagV=imagVfile)
    # df_varg_true = get_varg_from_complex(realV=realVfile, imagV=imagVfile)

    df_vmag_est = feather.read_feather(Vmagfile)
    df_vang_est = feather.read_feather(Vangfile)

    compare_vmag(
        df_vmag_true, df_vmag_est, 
        base_voltages, 
        bus="4643_22536.1", 
        to_file=f"outputs/{case}/test_known.png"
        )
    compare_vmag(
        df_vmag_true, df_vmag_est, 
        base_voltages, 
        bus="4643_22536.3", 
        to_file=f"outputs/{case}/test_unknown.png"
        )