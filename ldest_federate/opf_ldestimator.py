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
from lindistflow import get_Hmat, get_pq, get_v, get_vbase, get_nodes

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)


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


def matrix_to_numpy(admittance: List[List[Complex]]):
    "Convert list of list of our Complex type into a numpy matrix"
    return np.array([[x[0] + 1j * x[1] for x in row] for row in admittance])



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
        self.pub_power_real = self.vfed.register_publication(
            "pv_real", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_power_imag = self.vfed.register_publication(
            "pv_imag", h.HELICS_DATA_TYPE_STRING, ""
        )
        logger.debug("algorithm_parameters")
        logger.debug(algorithm_parameters)

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
        [branch, bus] = adapter.extract_info(topology)
        slack = topology.slack_bus[0]
        [source_bus, phase] = slack.split('.')
        branch_info, area_bus = area_info(
                        branch, bus, source_bus)
        bus_info = adapter.extract_injection(area_bus, topology.injections)
        

        # # Get the PV system ratings from the injections field of topology data
        # ratings = eqarray_to_xarray(
        #     topology.injections.power_real
        # ) + 1j * eqarray_to_xarray(topology.injections.power_imaginary)
        # pv_ratings = ratings[ratings.equipment_ids.str.startswith("PVSystem")]

        # # Get the load ratings from the injections field of topology data
        # load_ratings = ratings[ratings.equipment_ids.str.startswith("Load")]
        # bus_to_index = {v: i for i, v in enumerate(topology.base_voltage_magnitudes.ids)}
        # p_load = np.zeros(shape=(len(ids),))
        # q_load = np.zeros(shape=(len(ids),))
        # for i,v in enumerate(load_ratings.ids.data):
        #     p_load[bus_to_index[v]] = -load_ratings.values[i].real
        #     q_load[bus_to_index[v]] = -load_ratings.values[i].imag
        

        # p_pv = np.zeros(shape=(len(ids),))
        # q_pv = np.zeros(shape=(len(ids),))
        # for i,v in enumerate(pv_ratings.ids.data):
        #     p_pv[bus_to_index[v]] = pv_ratings.values[i].real
        #     q_pv[bus_to_index[v]] = pv_ratings.values[i].imag

        # v = measurement_to_xarray(topology.base_voltage_magnitudes)

        voltages_mag = None
        power_P = None
        power_Q = None
        while granted_time < h.HELICS_TIME_MAXTIME:
            
            logger.debug("granted_time")
            logger.debug(granted_time)
            
            # Keep granting time until all inputs arrive
            if not self.sub_voltages_mag.is_updated():
                granted_time = h.helicsFederateRequestTime(
                    self.vfed, h.HELICS_TIME_MAXTIME
                )
                logger.debug("Grant time for inputs to arrive")
                continue
            
            # Get input voltages data
            voltages_mag = VoltagesMagnitude.parse_obj(self.sub_voltages_mag.json)
            assert topology.base_voltage_magnitudes.ids == voltages_mag.ids
            bus_info = adapter.extract_voltages(bus_info, voltages_mag)
            
            # Get real and reactive power injections
            power_P = PowersReal.parse_obj(self.sub_power_P.json)
            power_Q = PowersImaginary.parse_obj(self.sub_power_Q.json)
            assert topology.base_voltage_magnitudes.ids == power_P.ids
            assert topology.base_voltage_magnitudes.ids == power_Q.ids
            bus_info = adapter.extract_powers(bus_info, power_P, power_Q)
            
            
            Hv, Hpq = get_Hmat(bus_info, branch_info, source_bus)
            pq = get_pq(bus_info, source_bus, SBASE=100.0e6)
            vmag, vslack = get_v(bus_info, source_bus)
            

            # compute per unit voltage magnitudes
            vbase = get_vbase(bus_info, topology.base_voltage_magnitudes)
            vmag_pu = vmag / vbase
            # vtrue = np.delete(vmag_pu, vslack)
            
            nodes = get_nodes(bus_info)
            node_ids = voltages_mag.ids
            nodes_ord = [nodes.tolist().index(nd) for nd in node_ids]
            

            ts = time.time()
            
            # select rows corresponding to voltages and columns 
            # corresponsing to slack bus voltage and node injections
            v0 = vmag_pu[vslack]
            H_check = np.hstack((Hv,Hpq))
            z = np.hstack((np.identity(len(vslack)), np.zeros(shape=(len(vslack),Hpq.shape[1]))))
            for i in range(len(vslack)):
                H_check = np.insert(H_check, vslack[i], z[i,:], axis=0)
            x_check = np.concatenate((v0, pq))
            v_linear = H_check @ x_check

            # Get the order of nodes correct
            v_true = vmag_pu[nodes_ord]
            v_est = v_linear[nodes_ord]
            
            import matplotlib.pyplot as plt
            fig,ax = plt.subplots(1,1,figsize=(20,12))
            ax.plot(range(len(v_true)), v_true, 'b--', lw=2.0, label='true')
            ax.plot(range(len(v_est)), np.sqrt(v_est), color='crimson', ls='dashed', lw=2.0, label='estimated')
            ax.set_xticks(list(range(len(v_true))), nodes[nodes_ord], fontsize=15, rotation=30)
            ax.tick_params(axis='y', labelsize=20)
            ax.set_ylabel("Voltage (in p.u.)", fontsize=20)
            ax.set_xlabel("Nodes in network", fontsize=20)
            ax.legend(fontsize=25, markerscale=2)
            fig.suptitle("Linearized voltages obtained from LinDistFlow", fontsize=30)
            fig.savefig("check_voltage_estimate.png")

            te = time.time()
            logger.debug(f"Estimator takes {(te-ts)} seconds")

            logger.info("end time: " + str(datetime.now()))

            # There should be a HELICS way to do this? Set resolution?
            previous_time = granted_time
            while (
                granted_time <= np.floor(previous_time) + 1
            ):  # This should avoid waiting a full 15 minutes
                granted_time = h.helicsFederateRequestTime(self.vfed, 1000)

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