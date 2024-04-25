"""
Description: Post analysis of state estimation results

Author: Rounak Meyur
"""

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

def get_vmag_from_complex(realV, imagV, topo):
    base_voltages = get_base_voltages(topo)
    voltage_real = feather.read_feather(realV)
    voltage_imag = feather.read_feather(imagV)
    df_voltages = np.abs(voltage_real.drop("time", axis=1) + 1j * voltage_imag.drop("time", axis=1))
    df_voltages["time"] = voltage_real["time"].apply(get_time)
    df_voltages = df_voltages.set_index("time")
    for node_name in df_voltages.columns:
        df_voltages[node_name] = df_voltages[node_name] / base_voltages[node_name]
    return  df_voltages

def get_vmag_in_pu(Vmag, topo):
    base_voltages = get_base_voltages(topo)
    df_voltages = feather.read_feather(Vmag)
    df_voltages["time"] = df_voltages["time"].apply(get_time)
    df_voltages = df_voltages.set_index("time")
    for node_name in df_voltages.columns:
        df_voltages[node_name] = df_voltages[node_name] / base_voltages[node_name]
    return df_voltages


def get_varg_from_complex(realV, imagV, deg=True):
    voltage_real = feather.read_feather(realV)
    voltage_imag = feather.read_feather(imagV)
    df_voltages = voltage_real.drop("time", axis=1) + 1j * voltage_imag.drop("time", axis=1)
    df_angles = pd.DataFrame(np.angle(df_voltages.values, deg=deg), columns=df_voltages.columns)
    df_angles["time"] = voltage_real["time"].apply(get_time)
    return df_angles.set_index("time")

def get_varg(Varg, deg=True):
    df_angles = feather.read_feather(Varg)
    df_angles["time"] = df_angles["time"].apply(get_time)
    df_angles = df_angles.set_index("time")
    if not deg:
        df_angles = df_angles * np.pi/180.0
    return df_angles


def compare_voltages(
        df_mag_true, df_mag_est,
        df_ang_true, df_ang_est, 
        time="15:30", 
        to_file = None, show=False, do_return=False,
        **kwargs):
    figsize = kwargs.get('figsize', (20, 10))
    constrained_layout = kwargs.get('constrained_layout', False)
    label_fontsize = kwargs.get('fontsize', 25)
    legend_fontsize = label_fontsize + 2
    ticklabel_fontsize = label_fontsize - 2
    title_fontsize = label_fontsize + 10
    suptitle_sfx = kwargs.get('suptitle_sfx',None)

    true_vmag = df_mag_true.loc[time]
    est_vmag = df_mag_est.loc[time]
    true_vang = df_ang_true.loc[time]
    est_vang = df_ang_est.loc[time]

    vmag_err = true_vmag - est_vmag
    vang_err = true_vang - est_vang

    fig, axs = plt.subplots(1, 2, figsize=figsize, constrained_layout=constrained_layout)
    
    axs[0].bar(vmag_err.index.values.tolist(),vmag_err)
    axs[1].bar(vang_err.index.values.tolist(),vang_err)

    axs[0].set_xlabel("Buses", fontsize=label_fontsize)
    axs[0].set_ylabel(f"Voltage Magnitudes (p.u.)", fontsize=label_fontsize)
    axs[0].set_xticklabels([])
    axs[0].tick_params(axis="y", labelsize=ticklabel_fontsize)
    axs[0].set_title("Voltage magnitude comparison", fontsize=label_fontsize)

    axs[1].set_xlabel("Buses", fontsize=label_fontsize)
    axs[1].set_ylabel(f"Voltage Angles (degrees)", fontsize=label_fontsize)
    axs[1].set_xticklabels([])
    axs[1].tick_params(axis="y", labelsize=ticklabel_fontsize)
    axs[1].set_title("Voltage angle comparison", fontsize=label_fontsize)
    suptitle = f"Voltage magnitude and angle comparison at {time} hours"
    fig.suptitle(suptitle, fontsize=title_fontsize)

    if to_file:
        fig.savefig(to_file, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)
    
    if do_return:
        return fig
    pass

def compare_vang(
        df_true, df_est, bus, 
        to_file = None, show=False, do_return=False,
        **kwargs):
    figsize = kwargs.get('figsize', (10, 10))
    constrained_layout = kwargs.get('constrained_layout', False)
    label_fontsize = kwargs.get('fontsize', 25)
    legend_fontsize = label_fontsize + 2
    ticklabel_fontsize = label_fontsize - 2
    title_fontsize = label_fontsize + 10
    suptitle_sfx = kwargs.get('suptitle_sfx',None)

    true_voltages = df_true[bus]
    est_voltages = df_est[bus]

    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=constrained_layout)
    xaxis = np.arange(df_est.shape[0])
    y_true = [true_voltages.iloc[i] for i in xaxis]
    y_est = [est_voltages.iloc[i] for i in xaxis]
    ax.plot(xaxis, y_true, color='red', marker='o', ls='dashed', lw=2.0, label="True")
    ax.plot(xaxis, y_est, color='blue', marker='o', ls='dashed', lw=2.0, label="Estimated")
    suptitle = f"Voltage angle comparison for bus {bus}"

    ax.set_xlabel("Time", fontsize=label_fontsize)
    ax.set_ylabel(f"Voltage Angle (degrees)", fontsize=label_fontsize)
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
    case = sys.argv[1]
    realVfile = os.path.join(f"outputs/{case}/voltage_real.feather")
    imagVfile = os.path.join(f"outputs/{case}/voltage_imag.feather")
    Vmagfile = os.path.join(f"outputs/{case}/voltage_mag.feather")
    Vangfile = os.path.join(f"outputs/{case}/voltage_angle.feather")
    topofile = os.path.join(f"outputs/{case}/topology.json")

    base_voltages = get_base_voltages(topology_file=topofile)
    df_vmag_true = get_vmag_from_complex(realV=realVfile, imagV=imagVfile, topo=topofile)
    df_vang_true = get_varg_from_complex(realV=realVfile, imagV=imagVfile, deg=True)
    df_vmag_est = get_vmag_in_pu(Vmagfile, topo=topofile)
    df_vang_est = get_varg(Vangfile, deg=True)

    compare_voltages(
        df_vmag_true, df_vmag_est,
        df_vang_true, df_vang_est, 
        to_file = f"outputs/{case}/result.png"
    )