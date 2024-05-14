from enum import Enum
import numpy as np
import math
import logging
from oedisi.types.data_types import VoltagesMagnitude

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)

"""
Author: Rounak Meyur
"""



def power_balance_relation(
        A : np.ndarray, 
        k_from : list, 
        k_to : list,
        col_offset : int,
        j : int
    ):
    for k in k_from:
        A[j, col_offset + k] = 1
    for k in k_to:
        A[j, col_offset + k] = -1

    return A


def voltage_cons_pri(
        Z : np.ndarray,
        A : np.ndarray, 
        idx, frm, to, 
        pii, qii, pij, qij, pik, qik, baseZ, 
        nbranch_ABC : int,
        brn_offset : int,
        bus_offset : int
    ):
    """
    Z: The matrix which relates voltage difference to power flows
    A: The matrix which relates voltage difference to voltages
    idx: The entry index for the branch
    """
    A[idx + brn_offset, frm + bus_offset] = 1
    A[idx + brn_offset, to + bus_offset] = -1

    # real power drop
    Z[idx + brn_offset, idx + nbranch_ABC * 0] = pii / baseZ
    Z[idx + brn_offset, idx + nbranch_ABC * 1] = pij / baseZ
    Z[idx + brn_offset, idx + nbranch_ABC * 2] = pik / baseZ
    # reactive power drop
    Z[idx + brn_offset, idx + nbranch_ABC * 3] = qii / baseZ
    Z[idx + brn_offset, idx + nbranch_ABC * 4] = qij / baseZ
    Z[idx + brn_offset, idx + nbranch_ABC * 5] = qik / baseZ
    
    return Z, A


def get_Hmat(
        bus_info : dict, 
        branch_info : dict, 
        source_bus : str
    ):

    slack_index = list(bus_info.keys()).index(source_bus)

    # System's base definition
    PRIMARY_V = 0.12
    SBASE = 100.0 # in MVA
    

    # Find the ABC phase and s1s2 phase triplex line and bus numbers
    nbranch_ABC = 0
    nbus_ABC = 0
    nbranch_s1s2 = 0
    nbus_s1s2 = 0
    secondary_model = ['TPX_LINE', 'SPLIT_PHASE']
    name = []
    for b_eq in branch_info:
        if branch_info[b_eq]['type'] in secondary_model:
            nbranch_s1s2 += 1
        else:
            nbranch_ABC += 1

    for b_eq in bus_info:
        name.append(b_eq)
        if bus_info[b_eq]['kv'] > PRIMARY_V:
            nbus_ABC += 1
        else:
            nbus_s1s2 += 1

    
    # Number of equality/inequality constraints (Injection equations (ABC) at each bus)
    #    #  Check if this is correct number or not:
    n_bus = nbus_ABC * 3 + nbus_s1s2  # Total Bus Number
    n_branch = nbranch_ABC * 3 + nbranch_s1s2  # Total Branch Number


    A1 = np.zeros(shape=(2*(nbus_ABC * 3 + nbus_s1s2), 2*(nbranch_ABC * 3 + nbranch_s1s2)))

    # # Define BFM constraints for both real and reactive power: Power flow conservaion
    # Constraint 1: Flow equation

    # sum(Sij) - sum(Sjk) == -sj

    for keyb, val_bus in bus_info.items():
        if keyb != source_bus:
            k_frm_3p = []
            k_to_3p = []
            k_frm_1p = []
            k_frm_1pa, k_frm_1pb, k_frm_1pc = [], [], []
            k_frm_1qa, k_frm_1qb, k_frm_1qc = [], [], []
            k_to_1p = []

            # Find bus idx in "from" of branch_sw_data
            ind_frm = 0
            ind_to = 0
            if val_bus['kv'] < PRIMARY_V:
                # if the bus is a part of a split phase transformer
                for key, val_br in branch_info.items():
                    if val_bus['idx'] == val_br['from']:
                        k_frm_1p.append(ind_frm - nbranch_ABC)

                    if val_bus['idx'] == val_br['to']:
                        k_to_1p.append(ind_to - nbranch_ABC)
                    ind_to += 1
                    ind_frm += 1
                # loc = (nbus_ABC * 3 + nbus_s1s2) + \
                #     (nbus_ABC * 6 + nbus_s1s2 * 2) + nbranch_ABC * 6
                # A_eq, b_eq = power_balance(A_eq, b_eq, k_frm_1p, k_to_1p, counteq, loc,
                #                             val_bus['idx'] + nbus_ABC * 3 + nbus_s1s2 + nbus_ABC * 5)
                # counteq += 1
                # A_eq, b_eq = power_balance(A_eq, b_eq, k_frm_1p, k_to_1p, counteq, loc + nbranch_s1s2,
                #                             val_bus['idx'] + nbus_ABC * 3 + nbus_s1s2 + nbus_ABC * 5 + nbus_s1s2)
                # counteq += 1
            else:
                # iterate through all the branches and find all whose 'from' or 'to' bus match the current bus 
                for key, val_br in branch_info.items():
                    if val_bus['idx'] == val_br['from']:
                        if bus_info[val_br['to_bus']]['kv'] > PRIMARY_V:
                            k_frm_3p.append(ind_frm)
                        else:
                            if key[-1] == 'a':
                                k_frm_1pa.append(
                                    nbranch_ABC * 6 + ind_frm - nbranch_ABC)
                                k_frm_1qa.append(
                                    nbranch_ABC * 3 + ind_frm - nbranch_ABC + nbranch_s1s2)
                            if key[-1] == 'b':
                                k_frm_1pb.append(
                                    nbranch_ABC * 5 + ind_frm - nbranch_ABC)
                                k_frm_1qb.append(
                                    nbranch_ABC * 2 + ind_frm - nbranch_ABC + nbranch_s1s2)
                            if key[-1] == 'c':
                                k_frm_1pc.append(
                                    nbranch_ABC * 4 + ind_frm - nbranch_ABC)
                                k_frm_1qc.append(
                                    nbranch_ABC * 1 + ind_frm - nbranch_ABC + nbranch_s1s2)

                    if val_bus['idx'] == val_br['to']:
                        if bus_info[val_br['fr_bus']]['kv'] > PRIMARY_V:
                            k_to_3p.append(ind_to)
                        else:
                            k_to_1p.append(ind_to - nbranch_ABC)
                    ind_to += 1
                    ind_frm += 1
                

                
                loc = 0
                # Finding the kfrms and ktos for the branches
                k_frm_A = k_frm_3p + k_frm_1pa
                k_frm_B = k_frm_3p + k_frm_1pb
                k_frm_C = k_frm_3p + k_frm_1pc
                k_to_A = k_to_B = k_to_C = k_to_3p

                # Real Power balance equations
                # Phase A
                A1 = power_balance_relation(
                    A1, k_frm_A, k_to_A, 
                    loc + nbranch_ABC * 0, 
                    val_bus['idx'] + nbus_ABC * 0
                    )
                # Phase B
                A1 = power_balance_relation(
                    A1, k_frm_B, k_to_B, 
                    loc + nbranch_ABC * 1, 
                    val_bus['idx'] + nbus_ABC * 1
                    )
                # Phase C
                A1 = power_balance_relation(
                    A1, k_frm_B, k_to_B, 
                    loc + nbranch_ABC * 2, 
                    val_bus['idx'] + nbus_ABC * 2
                    )

                
                # Finding the k_froms and k_tos for the branches
                k_frm_A = k_frm_3p + k_frm_1qa
                k_frm_B = k_frm_3p + k_frm_1qb
                k_frm_C = k_frm_3p + k_frm_1qc

                # Reactive Power balance equations
                # Phase A
                A1 = power_balance_relation(
                    A1, k_frm_A, k_to_A, 
                    loc + nbranch_ABC * 3, 
                    val_bus['idx'] + nbus_ABC * 3
                    )
                # Phase B
                A1 = power_balance_relation(
                    A1, k_frm_B, k_to_B, 
                    loc + nbranch_ABC * 4, 
                    val_bus['idx'] + nbus_ABC * 4
                    )
                # Phase C
                A1 = power_balance_relation(
                    A1, k_frm_B, k_to_B, 
                    loc + nbranch_ABC * 5, 
                    val_bus['idx'] + nbus_ABC * 5
                    )



    # Constraint 2: Voltage drop equation:
    # Vj = Vi - Zij Sij* - Sij Zij*


    A2 = np.zeros(shape = (n_branch, 2*n_branch))
    Av = np.zeros(shape = (n_branch, n_bus))

    # For Primary Nodes:
    idx = 0
    v_lim = []
    for k, val_br in branch_info.items():
        # compute base impedance
        basekV = bus_info[val_br['fr_bus']]['kv'] / np.sqrt(3)
        baseZ = basekV ** 2 / SBASE
        # Not writing voltage constraints for transformers
        if val_br['type'] not in secondary_model:
            z = np.asarray(val_br['zprim'])
            v_lim.append(val_br['from'])
            v_lim.append(val_br['to'])
            # Writing three phase voltage constraints
            # Phase A
            paa, qaa = -2 * z[0, 0][0], -2 * z[0, 0][1]
            pab, qab = -(- z[0, 1][0] + math.sqrt(3) * z[0, 1][1]), -(
                - z[0, 1][1] - math.sqrt(3) * z[0, 1][0])
            pac, qac = -(- z[0, 2][0] - math.sqrt(3) * z[0, 2][1]), -(
                - z[0, 2][1] + math.sqrt(3) * z[0, 2][0])
            A2, Av = voltage_cons_pri(
                A2, Av, 
                idx, val_br['from'], val_br['to'],
                paa, qaa, pab, qab, pac, qac, baseZ,
                nbranch_ABC, nbranch_ABC*0, nbus_ABC*0)

            
            # Phase B
            pbb, qbb = -2 * z[1, 1][0], -2 * z[1, 1][1]
            pba, qba = -(- z[0, 1][0] - math.sqrt(3) * z[0, 1][1]), -(
                - z[0, 1][1] + math.sqrt(3) * z[0, 1][0])
            pbc, qbc = -(- z[1, 2][0] + math.sqrt(3) * z[1, 2][1]), -(
                - z[1, 2][1] - math.sqrt(3) * z[1, 2][0])
            A2, Av = voltage_cons_pri(
                A2, Av, 
                idx, val_br['from'], val_br['to'],
                pba, qba, pbb, qbb, pbc, qbc, baseZ,
                nbranch_ABC, nbranch_ABC*1, nbus_ABC*1)
            
            
            # Phase C
            pcc, qcc = -2 * z[2, 2][0], -2 * z[2, 2][1]
            pca, qca = -(- z[0, 2][0] + math.sqrt(3) * z[0, 2][1]), -(
                - z[0, 2][1] - math.sqrt(3) * z[0, 2][0])
            pcb, qcb = -(- z[1, 2][0] - math.sqrt(3) * z[1, 2][1]), -(
                - z[0, 2][1] + math.sqrt(3) * z[1, 2][0])
            A2, Av = voltage_cons_pri(
                A2, Av, 
                idx, val_br['from'], val_br['to'],
                pca, qca, pcb, qcb, pcc, qcc, baseZ,
                nbranch_ABC, nbranch_ABC*2, nbus_ABC*2)
            
        idx += 1

    slack_node_idx = [slack_index, slack_index+nbus_ABC, slack_index+2*nbus_ABC]
    slack_node_idx_pq = slack_node_idx + [slack_index+n_bus, slack_index+n_bus+nbus_ABC, slack_index+n_bus+2*nbus_ABC]
    Av0 = Av[:,slack_node_idx]
    Avr = np.delete(Av, slack_node_idx, axis=1)
    Avr_inv = np.linalg.inv(Avr)

    H21 = - (Avr_inv @ Av0)
    H22 = (Avr_inv @ A2)

    H32 = np.delete(A1, slack_node_idx_pq, axis=0)
    H32_inv = np.linalg.inv(H32)
    return H21, H22@H32_inv


def get_pq(
        bus_info : dict, 
        source_bus : str, 
        SBASE : float = 100.0e6
    ):
    n_bus = len(bus_info) - 1
    pq = np.zeros(shape=(6*n_bus,))
    count = 0
    for keyb, val_bus in bus_info.items():
        if keyb != source_bus:
            # Real power injection at a bus
            # Phase A Real Power
            pq[count + n_bus*0] = val_bus['pv'][0][0] + val_bus['pq'][0][0]
            # Phase B Real Power
            pq[count + n_bus*1] = val_bus['pv'][1][0] + val_bus['pq'][1][0]
            # Phase C Real Power
            pq[count + n_bus*2] = val_bus['pv'][2][0] + val_bus['pq'][2][0]
            
            
            # Phase A Reactive Power
            pq[count + n_bus*3] = val_bus['pv'][0][1] + val_bus['pq'][0][1]
            # Phase B Reactive Power
            pq[count + n_bus*4] = val_bus['pv'][1][1] + val_bus['pq'][1][1]
            # Phase C Reactive Power
            pq[count + n_bus*5] = val_bus['pv'][2][1] + val_bus['pq'][2][1]

            count += 1
    return pq / SBASE

def get_v(
        bus_info : dict, 
        source_bus : str
        ):
    n_bus = len(bus_info)
    v = np.zeros(shape=(3*n_bus,))
    slack_index = []
    count = 0
    for keyb, val_bus in bus_info.items():
        if keyb == source_bus:
            slack_index = [count, count+n_bus, count+2*n_bus]
        v[count + n_bus*0] = val_bus['vmag'][0]
        v[count + n_bus*1] = val_bus['vmag'][1]
        v[count + n_bus*2] = val_bus['vmag'][2]
        count += 1
    return v, slack_index

def get_vbase(
        bus_info : dict, 
        voltages : VoltagesMagnitude
    ):
    n_bus = len(bus_info)
    v = np.zeros(shape=(3*n_bus,))
    ids = voltages.ids
    vals = voltages.values
    
    for keyb, val_bus in bus_info.items():
        idx = val_bus["idx"]
        v[idx + n_bus*0] = vals[ids.index(f"{keyb}.1")]
        v[idx + n_bus*1] = vals[ids.index(f"{keyb}.2")]
        v[idx + n_bus*2] = vals[ids.index(f"{keyb}.3")]
        
    return v

def get_nodes(bus_info:dict) -> list:
    n_bus = len(bus_info)
    nodes = np.empty(shape=(3*n_bus,), dtype='<U20')
    for keyb, val_bus in bus_info.items():
        idx = val_bus["idx"]
        nodes[idx + n_bus*0] = f"{keyb}.1"
        nodes[idx + n_bus*1] = f"{keyb}.2"
        nodes[idx + n_bus*2] = f"{keyb}.3"
    return nodes