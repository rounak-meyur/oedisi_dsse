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

def compare_vmag(df_true, df_est, bus, base_voltages):
    true_voltages = df_true[bus] / base_voltages[bus]
    est_voltages = df_est[bus] / base_voltages[bus]

    return

if __name__ == "__main__":
    case = sys.argv[1]
    realVfile = os.path.join(f"./outputs/{case}/voltage_real.feather")
    imagVfile = os.path.join(f"./outputs/{case}/voltage_imag.feather")
    Vmagfile = os.path.join(f"./outputs/{case}/voltage_mag.feather")
    Vangfile = os.path.join(f"./outputs/{case}/voltage_angle.feather")
    topofile = os.path.join(f"./outputs/{case}/topology.json")

    base_voltages = get_base_voltages(topology_file=topofile)
    df_vmag_true = get_vmag_from_complex(realV=realVfile, imagV=imagVfile)
    df_varg_true = get_varg_from_complex(realV=realVfile, imagV=imagVfile)

    df_vmag_est = feather.read_feather(Vmagfile)
    df_vang_est = feather.read_feather(Vangfile)