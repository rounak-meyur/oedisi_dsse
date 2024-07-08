"""
Created on Tue March 30 11:11:21 2021

Author:  Rounak Meyur
Description: State Estimation with Linearized Distribution Flow equations
"""

import cmath
import warnings
import logging
import helics as h
import json
import numpy as np
import time
import pandas as pd
from pydantic import BaseModel
from enum import Enum
from typing import List, Optional, Union
from datetime import datetime
from oedisi.types.data_types import (
    AdmittanceMatrix,
    EquipmentNodeArray,
    MeasurementArray,
    Topology,
    Complex,
    VoltagesMagnitude,
    Injection,
    CommandList,
    PowersReal,
    PowersImaginary,
    AdmittanceSparse,
    VoltagesMagnitude,
    Command,
)
from scipy.sparse import csc_matrix, coo_matrix, diags, vstack, hstack
from scipy.sparse.linalg import svds, inv
import xarray as xr

import adapter
from area import area_info
from pv_detect_dsse import (
    get_Hmat, 
    get_pq, get_pq_forecast, get_pv,
    get_v, get_vbase, get_nodes
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)

def add_noise(array, mean=0, std=0.1):
  """Adds random noise to an array.

  Args:
    array: The array to add noise to.
    mean: The mean of the noise.
    std: The standard deviation of the noise.

  Returns:
    The array with noise added.
  """

  noise = np.random.normal(mean, std, size=array.shape)
  return array + noise


class UnitSystem(str, Enum):
    SI = "SI"
    PER_UNIT = "PER_UNIT"


class EstimatorParameters(BaseModel):
    Vmax: float = 1.05  # + 0.005  # Upper limit
    Vmin_act: float = 0.95  # + 0.005
    Vmin: float = 0.95  # + 0.005 + 0.002  # Lower limit\
    # Linearized equation is overestimating. So increase the lower limit by 0.005.
    # The problem is not solved to the optimal, so increase another 0.002.
    alpha: float = 0.5  #  learning rate for dual and primal
    epsilon: float = 1e-8  # decay of duals
    nu: float = 0.016  # A trivial penalty factor
    ratio_t_k: int = 1000
    units: UnitSystem = UnitSystem.PER_UNIT
    base_power: Optional[float] = 100.0


class Sensors(float):
    SBASE = 1e6
    ### add noises:
    Vmeas_sigma = 0.05e-5  # 5% * nominal values 1-2% (make it general based on Sbase/Vbase
    forecast_sigma = 0.005e-5  # 5-10% error (based on forecast uncertainty)
    Pmeas_sigma = 0.003e-5 # 0.1-0.5% error as the customers are billed (>= 15 min interval)
    Qmeas_sigma = 0.01e-5 # 1% of nominal KVA error probably?


def eqarray_to_xarray(eq: EquipmentNodeArray):
    return xr.DataArray(
        eq.values,
        dims=("eqnode",),
        coords={
            "equipment_ids": ("eqnode", eq.equipment_ids),
            "ids": ("eqnode", eq.ids),
        },
    )


def measurement_to_xarray(eq: MeasurementArray):
    return xr.DataArray(
        eq.values, 
        coords={"ids": eq.ids}
    )

def xarray_to_dict(data):
    """Convert xarray to dict with values and ids for JSON serialization."""
    coords = {key: list(data.coords[key].data) for key in data.coords.keys()}
    return {"values": list(data.data), **coords}

def matrix_to_numpy(admittance: List[List[Complex]]):
    "Convert list of list of our Complex type into a numpy matrix"
    return np.array([[x[0] + 1j * x[1] for x in row] for row in admittance])

def sub_to_dict(sub):
    return {k:v for (k,v) in zip(sub.ids,sub.values)}

class EstimatorFederate:
    "Estimator federate. Wraps PV generation estimation with pubs and subs"

    def __init__(self, federate_name, algorithm_parameters, input_mapping):
        "Initializes federate with name and remaps input into subscriptions"
        deltat = 0.1

        self.algorithm_parameters = algorithm_parameters

        # Create Federate Info object that describes the federate properties #
        fedinfo = h.helicsCreateFederateInfo()

        fedinfo.core_name = federate_name
        fedinfo.core_type = h.HELICS_CORE_TYPE_ZMQ
        fedinfo.core_init = "--federates=1"
        h.helicsFederateInfoSetTimeProperty(
            fedinfo, h.helics_property_time_delta, deltat
        )

        self.vfed = h.helicsCreateValueFederate(federate_name, fedinfo)
        logger.info("Value federate created")

        # Register the subscription
        self.sub_voltages_mag = self.vfed.register_subscription(
            input_mapping["voltages_magnitude"], "V"
        )
        self.sub_power_P = self.vfed.register_subscription(
            input_mapping["powers_real"], "W"
        )
        self.sub_power_Q = self.vfed.register_subscription(
            input_mapping["powers_imaginary"], "W"
        )
        self.sub_topology = self.vfed.register_subscription(
            input_mapping["topology"], ""
        )

        # Register the publications
        self.pub_pv_real = self.vfed.register_publication(
            "pv_real", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_pv_imag = self.vfed.register_publication(
            "pv_imag", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_pv_real_actual = self.vfed.register_publication(
            "pv_real_actual", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_pv_imag_actual = self.vfed.register_publication(
            "pv_imag_actual", h.HELICS_DATA_TYPE_STRING, ""
        )
        # logger.debug("algorithm_parameters")
        # logger.debug(algorithm_parameters)

    def run(self):
        "Enter execution and exchange data"
        # Enter execution mode #
        self.vfed.enter_executing_mode()
        logger.info("Entering execution mode")

        # granted_time = h.helicsFederateRequestTime(self.vfed, h.HELICS_TIME_MAXTIME)
        granted_time = h.helicsFederateRequestTime(self.vfed, 1000)

        topology = Topology.parse_obj(self.sub_topology.json)
        ids = topology.base_voltage_magnitudes.ids
        logger.info("Topology has been read")
        
        # Read the bus and branches
        SBASE = Sensors.SBASE
        [branch, bus] = adapter.extract_info(topology)
        slack = topology.slack_bus[0]
        [source_bus, phase] = slack.split('.')
        branch_info, area_bus = area_info(
                        branch, bus, source_bus)
        bus_info = adapter.extract_injection(
            area_bus, topology.injections)
        
        # with open("bus_info.json", 'w') as f:
        #     json.dump(bus_info, f)
        # with open("branch_info.json", 'w') as f:
        #     json.dump(branch_info, f)
        
        baseV = sub_to_dict(topology.base_voltage_magnitudes)
        with open("base_voltages.json", 'w') as f:
            json.dump(baseV, f)

        pq_load = get_pq_forecast(bus_info, source_bus, SBASE=SBASE)
        pq_forecast = add_noise(pq_load, mean=0, std=Sensors.forecast_sigma)

        voltages_mag = None
        power_P = None
        power_Q = None
        while granted_time < h.HELICS_TIME_MAXTIME:
            
            logger.info(f"granted_time: {granted_time}")
            
            # Keep granting time until all inputs arrive
            if not self.sub_voltages_mag.is_updated():
                granted_time = h.helicsFederateRequestTime(
                    self.vfed, h.HELICS_TIME_MAXTIME
                )
                logger.debug("Grant time for inputs to arrive")
                continue
            
            # Get input voltages data
            voltages_mag = VoltagesMagnitude.parse_obj(self.sub_voltages_mag.json)
            current_time = voltages_mag.time
            logger.info(current_time)
            bus_info = adapter.extract_voltages(bus_info, voltages_mag)
            
            # Get real and reactive power injections
            power_P = PowersReal.parse_obj(self.sub_power_P.json)
            power_Q = PowersImaginary.parse_obj(self.sub_power_Q.json)
            bus_info = adapter.extract_powers(bus_info, power_P, power_Q)

            # save the measurements as JSON files
            # vmag_info = sub_to_dict(voltages_mag)
            # pinj_info = sub_to_dict(power_P)
            # qinj_info = sub_to_dict(power_Q)
            # with open("vmag.json", 'w') as f:
            #     json.dump(vmag_info, f)
            # with open("pinj.json", 'w') as f:
            #     json.dump(pinj_info, f)
            # with open("qinj.json", 'w') as f:
            #     json.dump(qinj_info, f)
            
            H = get_Hmat(bus_info, branch_info, source_bus, SBASE=SBASE)
            pq = get_pq(bus_info, source_bus, SBASE=SBASE)
            vmag, vslack = get_v(bus_info, source_bus)
            

            # compute per unit voltage magnitudes
            vbase = get_vbase(bus_info, topology.base_voltage_magnitudes)
            vmag_pu = vmag / vbase
            nodes = [n for k,n in enumerate(get_nodes(bus_info)) if k not in vslack]
            logger.debug(f"Actual slack bus : {vmag_pu[vslack]}")
            
            ############################ estimation #############################
            ts = time.time()

            V_W = np.array([1/(Sensors.Vmeas_sigma**2)]*len(vmag_pu))
            Pl_W = np.array([1 / (Sensors.forecast_sigma ** 2)] * len(pq_load))
            Pinj_W = np.array([1 / (Sensors.Pmeas_sigma ** 2)] * int(len(pq) / 2) )
            Qinj_W = np.array([1 / (Sensors.Qmeas_sigma ** 2)] * int(len(pq) / 2) )
            

            # increase weight of substation bus measurement
            V_W[vslack] = 1e7*V_W[vslack]
            # logger.debug(V_W)
            
            # combine the weights
            W_array = np.hstack((V_W, Pl_W, Pinj_W, Qinj_W))
            W = np.diag(W_array)

            # perform estimation
            Z_meas = np.hstack((np.square(vmag_pu), pq_forecast, pq))
            G = H.T @ W @ H
            G_inv = np.linalg.inv(G)
            x_est = G_inv @ H.T @ W @ Z_meas

            te = time.time()
            #####################################################################

            #### Residual computation
            logger.debug(f"Slack bus estimates : {x_est[:len(vslack)]}")
            z_est = H @ x_est

            actual_pv = get_pv(bus_info, source_bus, SBASE=SBASE)
            Ppv_act = actual_pv[:int(len(pq)/2)]
            Qpv_act = actual_pv[int(len(pq)/2):]
            
            Ppv_est = x_est[len(vslack) + len(pq): len(vslack) + len(pq) +int(len(pq_load)/2) ]
            Qpv_est = x_est[len(vslack) + len(pq) +int(len(pq_load)/2): ]
            
            estimated_realPV_nodes = xr.DataArray(Ppv_est, coords={"ids":nodes})
            estimated_reactivePV_nodes = xr.DataArray(Qpv_est, coords={"ids":nodes})
            actual_realPV_nodes = xr.DataArray(Ppv_act, coords={"ids":nodes})
            actual_reactivePV_nodes = xr.DataArray(Qpv_act, coords={"ids":nodes})

            self.pub_pv_real.publish(
                MeasurementArray(
                    **xarray_to_dict(estimated_realPV_nodes),
                    time=current_time, 
                    units="pu"
                ).json()
            )
            self.pub_pv_imag.publish(
                MeasurementArray(
                    **xarray_to_dict(estimated_reactivePV_nodes),
                    time=current_time, 
                    units="pu"
                ).json()
            )
            self.pub_pv_real_actual.publish(
                MeasurementArray(
                    **xarray_to_dict(actual_realPV_nodes),
                    time=current_time, 
                    units="pu"
                ).json()
            )
            self.pub_pv_imag_actual.publish(
                MeasurementArray(
                    **xarray_to_dict(actual_reactivePV_nodes),
                    time=current_time, 
                    units="pu"
                ).json()
            )
            

            logger.info(f"Estimator takes {(te-ts)} seconds")

        self.destroy()

    def destroy(self):
        "Finalize and destroy the federates"
        h.helicsFederateDisconnect(self.vfed)
        logger.info("Federate disconnected")

        h.helicsFederateFree(self.vfed)
        h.helicsCloseLibrary()


if __name__ == "__main__":
    with open("static_inputs.json") as f:
        config = json.load(f)
        federate_name = config["name"]
        if "algorithm_parameters" in config:
            parameters = EstimatorParameters.parse_obj(config["algorithm_parameters"])
        else:
            parameters = EstimatorParameters.parse_obj({})

    with open("input_mapping.json") as f:
        input_mapping = json.load(f)

    sfed = EstimatorFederate(federate_name, parameters, input_mapping)
    sfed.run()