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

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)


class UnitSystem(str, Enum):
    SI = "SI"
    PER_UNIT = "PER_UNIT"

class Sensors(float):
    ### add noises:
    Vmeas_sigma = 0.05  # 5% * nominal values 1-2% (make it general based on Sbase/Vbase
    forecast_sigma = 0.05  # 5-10% error (based on forecast uncertainty)
    Pmeas_sigma = 0.003 # 0.1-0.5% error as the customers are billed (>= 15 min interval)
    Qmeas_sigma = 0.01 # 1% of nominal KVA error probably?

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

def xarray_to_dict(data):
    """Convert xarray to dict with values and ids for JSON serialization."""
    coords = {key: list(data.coords[key].data) for key in data.coords.keys()}
    return {"values": list(data.data), **coords}


def matrix_to_numpy(admittance: List[List[Complex]]):
    "Convert list of list of our Complex type into a numpy matrix"
    return np.array([[x[0] + 1j * x[1] for x in row] for row in admittance])


def get_y(admittance: Union[AdmittanceMatrix, AdmittanceSparse], ids: List[str]):
    """
    Converts the admittance matrix from the topology.json file to a numpy array
    """
    if type(admittance) == AdmittanceMatrix:
        assert ids == admittance.ids
        return matrix_to_numpy(admittance.admittance_matrix)
    elif type(admittance) == AdmittanceSparse:
        node_map = {name: i for (i, name) in enumerate(ids)}
        return coo_matrix(
            (
                [v[0] + 1j * v[1] for v in admittance.admittance_list],
                (
                    [node_map[r] for r in admittance.from_equipment],
                    [node_map[c] for c in admittance.to_equipment],
                ),
            )
        ).toarray()


def getLinearModel(YLL, YL0, V0):
    """
    Calculates components in equation (3)
    [Decentralized Low-Rank State Estimation for Power Distribution Systems]
    ----
    Arguments:
        YLL, YL0: A part of Y matrix. Defined above equation (4)
        V0: True voltages of slack buses.
    Returns:
        w, w_mag, My, Ky: w, |w|, N, K in (3)
    
    Reference: A. Bernstein and E. Dall'Anese, "Linear power-flow models in multiphase distribution networks," 
    2017 IEEE PES Innovative Smart Grid Technologies Conference Europe (ISGT-Europe), Turin, Italy, 2017, pp. 1-6, 
    doi: 10.1109/ISGTEurope.2017.8260205.
    """
    #
    Nnode = YLL.shape[0]
    YLLi = inv(YLL)  # quicker than spsolve(YLL_s, eye(YLL_s.shape[0])) by testing
    Ky0 = - (YLLi @ YL0)
    w = Ky0 @ V0
    w_mag = np.abs(w)
    Vk = w
    Vk_mag = np.abs(Vk)
    prody = YLLi @ inv(diags(np.conjugate(Vk).squeeze(), format="csc"))
    My = hstack([prody, -1j * prody])

    Ky = (
        inv(diags(np.conjugate(Vk_mag).squeeze(), format="csc"))
        @ (diags(np.conjugate(Vk).squeeze(), format="csc") @ My).real
    )
    Ky1, Ky2 = Ky[:, :Nnode], Ky[:, Nnode:]
    # Ky1, Ky2 are converted to dense, since they are not sparse
    return np.abs(Ky0).toarray(), Ky1.toarray(), Ky2.toarray(), w_mag.reshape(-1)


def get_H(E, R, X):
    """
    Construct the measurement to state matrix. The rows denote the measurements and the columns are the
    states of the system. 
    Measurement vector: z = [V0 Vn Pinj Qinj Pload Qload]
    State vector: x = [V0 Pinj Qinj pPV qPV]

    #                      (V0)       (Pall)          (Qall)           pPV             qPV
    #                   -------------------------------------------------------------------------
    #       (V0)        |   I   |       0       |       0       |       0       |       0       |
    #                   -------------------------------------------------------------------------
    #                   |       |               |               |               |               |
    #       (Vn)        |   E   |       R       |       X       |       0       |       0       |
    #                   |       |               |               |               |               |
    #                   -------------------------------------------------------------------------
    #                   |       |               |               |               |               |
    #     (Pall)        |   0   |       I       |       0       |       0       |       0       |
    #                   |       |               |               |               |               |
    #                   -------------------------------------------------------------------------
    #                   |       |               |               |               |               |
    #     (Qall)        |   0   |       0       |       I       |       0       |       0       |
    #                   |       |               |               |               |               |
    #                   -------------------------------------------------------------------------
    #                   |       |               |               |               |               |
    #       (Pl)        |   0   |      -I       |       0       |       I       |       0       |
    #                   |       |               |               |               |               |
    #                   -------------------------------------------------------------------------
    #                   |       |               |               |               |               |
    #       (Ql)        |   0   |       0       |      -I       |       0       |       I       |
    #                   |       |               |               |               |               |
    """
    n = E.shape[0]
    s = E.shape[1]
    m1 = R.shape[1]
    m2 = X.shape[1]

    zns = np.zeros(shape=(n,s))
    zsm1 = np.zeros(shape=(s,m1+s))
    zsm2 = np.zeros(shape=(s,m2+s))
    znm1 = np.zeros(shape=(n,m1+s))
    znm2 = np.zeros(shape=(n,m2+s))
    zm1m1 = np.zeros(shape=(m1+s,m1+s))
    zm1m2 = np.zeros(shape=(m1+s,m2+s))
    zm2m1 = np.zeros(shape=(m2+s,m1+s))
    zm2m2 = np.zeros(shape=(m2+s,m2+s))
    
    H1 = np.hstack((np.identity(s), zsm1, zsm2, zsm1, zsm2))
    H2 = np.hstack((E, zns, -R, zns, -X, znm1, znm2))
    H3 = np.hstack((zsm1.T, np.identity(m1+s), zm1m2, zm1m1, zm1m2))
    H4 = np.hstack((zsm2.T, zm2m1, np.identity(m2+s), zm2m1, zm2m2))
    H5 = np.hstack((zsm1.T, -np.identity(m1+s), zm1m2, np.identity(m1+s), zm1m2))
    H6 = np.hstack((zsm2.T, zm2m1, -np.identity(m2+s), zm2m1, np.identity(m2+s)))

    H = np.vstack((H1, H2, H3, H4, H5, H6))
    return H
        




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
        self.pub_pv_real_est = self.vfed.register_publication(
            "pv_real_estimated", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_pv_imag_est = self.vfed.register_publication(
            "pv_imag_estimated", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_pv_real_actual = self.vfed.register_publication(
            "pv_real_actual", h.HELICS_DATA_TYPE_STRING, ""
        )
        self.pub_pv_imag_actual = self.vfed.register_publication(
            "pv_imag_actual", h.HELICS_DATA_TYPE_STRING, ""
        )

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
        slack_index = None
        if not isinstance(topology.admittance, AdmittanceMatrix) and not isinstance(
            topology.admittance, AdmittanceSparse
        ):
            raise "PV generation estimator algorithm expects AdmittanceMatrix/Sparse as input"

        # Get the slack bus index
        for i in range(len(ids)):
            if ids[i] == topology.slack_bus[0]:
                slack_index = i

        if not slack_index == 0:
            raise "Slack index is not 0"
        
        slack_bus = [slack_index, slack_index + 1, slack_index + 2]
        
        # Get the admittance matrix and compute the 4 partitions
        # Divide the admittance matrix into 4 partitions as follows:
        #      |             |
        #      |  Y00   Y0L  |
        #  Y = |             |
        #      |  YL0   YLL  |
        #      |             |
        Y = get_y(topology.admittance, topology.base_voltage_magnitudes.ids)
        if self.algorithm_parameters.units == UnitSystem.PER_UNIT:
            base_power = 100
            if self.algorithm_parameters.base_power != None:
                base_power = self.algorithm_parameters.base_power
            Y = (
                np.array(topology.base_voltage_magnitudes.values).reshape(1, -1)
                * Y
                * np.array(topology.base_voltage_magnitudes.values).reshape(-1, 1)
                / (base_power * 1000)
            )
        self.YLL = csc_matrix(
            np.delete(np.delete(Y, slack_bus, axis=0), slack_bus, axis=1)
        )
        self.YL0 = csc_matrix(np.delete(Y, slack_bus, axis=0)[:, slack_bus])
        self.Y = csc_matrix(Y)
        del Y

        # Get the PV system ratings from the injections field of topology data
        ratings = eqarray_to_xarray(
            topology.injections.power_real
        ) + 1j * eqarray_to_xarray(topology.injections.power_imaginary)
        pv_ratings = ratings[ratings.equipment_ids.str.startswith("PVSystem")]

        # Get the load ratings from the injections field of topology data
        load_ratings = ratings[ratings.equipment_ids.str.startswith("Load")]
        bus_to_index = {v: i for i, v in enumerate(topology.base_voltage_magnitudes.ids)}
        
        p_load = np.zeros(shape=(len(ids),))
        q_load = np.zeros(shape=(len(ids),))
        for i,v in enumerate(load_ratings.ids.data):
            p_load[bus_to_index[v]] = -load_ratings.values[i].real
            q_load[bus_to_index[v]] = -load_ratings.values[i].imag
        pl_forecast = add_noise(p_load, mean=0, std=Sensors.forecast_sigma)
        ql_forecast = add_noise(q_load, mean=0, std=Sensors.forecast_sigma)
        

        p_pv = np.zeros(shape=(len(ids),))
        q_pv = np.zeros(shape=(len(ids),))
        for i,v in enumerate(pv_ratings.ids.data):
            p_pv[bus_to_index[v]] = pv_ratings.values[i].real
            q_pv[bus_to_index[v]] = pv_ratings.values[i].imag

        # v = measurement_to_xarray(topology.base_voltage_magnitudes)
        voltages_mag = None
        power_P = None
        power_Q = None
        while granted_time < h.HELICS_TIME_MAXTIME:
            
            logger.debug(f"granted_time : {granted_time}")
            
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
            voltages = measurement_to_xarray(voltages_mag)
            V0 = voltages.loc[topology.slack_bus].data
            
            self.V0 = (
                V0 / np.array(topology.base_voltage_magnitudes.values)[slack_bus]
            ).reshape(3, -1)

            # Get real and reactive power injections
            power_P = measurement_to_xarray(
                PowersReal.parse_obj(
                    self.sub_power_P.json
                    )
                )
            power_Q = measurement_to_xarray(
                PowersImaginary.parse_obj(
                    self.sub_power_Q.json
                    )
                )
            
            # get the vector of measurements in correct order
            base_voltages = np.array(topology.base_voltage_magnitudes.values)
            base_power = 100.0
            
            vmag = voltages.loc[ids].data / base_voltages
            vmag_sq = np.square(vmag)
            P_inj = power_P.loc[ids].data / base_power
            Q_inj = power_Q.loc[ids].data / base_power
            
            Ppv_act = p_pv / base_power
            Qpv_act = q_pv / base_power
            
            P_load = pl_forecast / base_power
            Q_load = ql_forecast / base_power
            
            Z_meas = np.concatenate((vmag_sq, P_inj, Q_inj, P_load, Q_load))
            
            ts = time.time()
            
            

            # Check the linear model
            self.E, self.G, self.H, self.w_mag = getLinearModel(self.YLL, self.YL0, self.V0)
            Hmat = get_H(self.E, self.G, self.H)
            
            
            # select rows corresponding to voltages and columns 
            # corresponsing to slack bus voltage and node injections
            # H_check = Hmat[:len(vmag)]
            # x_check = np.concatenate((self.V0.reshape(-1), P_inj, Q_inj, pPV_true, qPV_true))
            # v_est = H_check @ x_check

            # Node IDs in the network

            V_W = np.array([1/(Sensors.Vmeas_sigma**2)]*len(vmag_sq))
            Pl_W = np.array([1 / (Sensors.forecast_sigma ** 2)] * len(P_load))
            Ql_W = np.array([1 / (Sensors.forecast_sigma ** 2)] * len(Q_load))
            Pinj_W = np.array([1 / (Sensors.Pmeas_sigma ** 2)] * len(P_inj) )
            Qinj_W = np.array([1 / (Sensors.Qmeas_sigma ** 2)] * len(Q_inj) )
            W_array = np.hstack((V_W, Pl_W, Ql_W, Pinj_W, Qinj_W))
            W = np.diag(W_array)
            
            G = Hmat.T @ W @ Hmat
            G_inv = np.linalg.inv(G)

            
            x_est = G_inv @ Hmat.T @ W @ Z_meas

            Ppv_est = x_est[len(slack_bus)+len(P_inj)+len(Q_inj):len(slack_bus)+len(P_inj)+len(Q_inj)+len(P_inj)]
            Qpv_est = x_est[len(slack_bus)+len(P_inj)+len(Q_inj)+len(P_inj):]

            te = time.time()

            ##########################
            logger.debug(Ppv_est)
            logger.debug(Ppv_act)
            # logger.debug("Slack voltage estimates")
            # logger.debug(x_est[:len(slack_bus)])
            #### Residual computation
            Z_est = Hmat @ x_est
            # logger.debug("Residual computation")
            # logger.debug(Z_meas - Z_est)
            ##########################
            
            nodes = [n for n in ids]
            estimated_realPV_nodes = xr.DataArray(Ppv_est, coords={"ids":nodes})
            estimated_reactivePV_nodes = xr.DataArray(Qpv_est, coords={"ids":nodes})
            actual_realPV_nodes = xr.DataArray(Ppv_act, coords={"ids":nodes})
            actual_reactivePV_nodes = xr.DataArray(Qpv_act, coords={"ids":nodes})

            self.pub_pv_real_est.publish(
                MeasurementArray(
                    **xarray_to_dict(estimated_realPV_nodes),
                    time=current_time,
                    units="kW"
                ).json()
            )
            self.pub_pv_imag_est.publish(
                MeasurementArray(
                    **xarray_to_dict(estimated_reactivePV_nodes),
                    time=current_time,
                    units="kVAR"
                ).json()
            )
            self.pub_pv_real_actual.publish(
                MeasurementArray(
                    **xarray_to_dict(actual_realPV_nodes),
                    time=current_time,
                    units="kW"
                ).json()
            )
            self.pub_pv_imag_actual.publish(
                MeasurementArray(
                    **xarray_to_dict(actual_reactivePV_nodes),
                    time=current_time,
                    units="kVAR"
                ).json()
            )

            logger.debug(f"Estimator takes {(te-ts)} seconds")

            logger.info("end time: " + str(datetime.now()))

            # There should be a HELICS way to do this? Set resolution?
            previous_time = granted_time
            # while (
            #     granted_time <= np.floor(previous_time) + 1
            # ):  # This should avoid waiting a full 15 minutes
            #     granted_time = h.helicsFederateRequestTime(self.vfed, 1000)

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