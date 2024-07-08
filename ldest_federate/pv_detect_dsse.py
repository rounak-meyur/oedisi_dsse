import numpy as np
import math
import logging
from oedisi.types.data_types import VoltagesMagnitude

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)

"""
Author: Rounak Meyur

Adapted from code written by Shiva Poudel and Monish Mukherjee

Description: Builds a matrix to relate the linear voltage magnitude estimates of all nodes to the 
real and reactive power injections at the nodes.
"""


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

def power_balance_relation(
        A: np.ndarray,
        k_from: list,
        k_to: list,
        col_offset: int,
        j: int
):
    for k in k_from:
        A[j, col_offset + k] = 1
    for k in k_to:
        A[j, col_offset + k] = -1

    return A


def voltage_cons_pri(
        Z: np.ndarray,
        A: np.ndarray,
        idx, frm, to,
        pii, qii, pij, qij, pik, qik, baseZ,
        nbranch_ABC: int,
        brn_offset: int,
        bus_offset: int
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
        bus_info: dict,
        branch_info: dict,
        source_bus: str,
        SBASE: float = 1e6
):
    slack_index = list(bus_info.keys()).index(source_bus)

    # System's base definition
    PRIMARY_V = 0.12
    _SBASE = SBASE / 1e6

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
    

    # Constraint 2: Voltage drop equation:
    # Vj = Vi - Zij Sij* - Sij Zij*
    A2 = np.zeros(shape=(n_branch, 2 * n_branch))
    A1 = np.zeros(shape=(n_branch, n_bus))

    # For Primary Nodes:
    idx = 0
    v_lim = []
    for k, val_br in branch_info.items():
        # compute base impedance
        basekV = bus_info[val_br['to_bus']]['kv']
        baseZ = (basekV ** 2) / (_SBASE)

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
            A2, A1 = voltage_cons_pri(
                A2, A1,
                idx, val_br['from'], val_br['to'],
                paa, qaa, pab, qab, pac, qac, baseZ,
                nbranch_ABC, nbranch_ABC * 0, nbus_ABC * 0)

            # Phase B
            pbb, qbb = -2 * z[1, 1][0], -2 * z[1, 1][1]
            pba, qba = -(- z[0, 1][0] - math.sqrt(3) * z[0, 1][1]), -(
                    - z[0, 1][1] + math.sqrt(3) * z[0, 1][0])
            pbc, qbc = -(- z[1, 2][0] + math.sqrt(3) * z[1, 2][1]), -(
                    - z[1, 2][1] - math.sqrt(3) * z[1, 2][0])
            A2, A1 = voltage_cons_pri(
                A2, A1,
                idx, val_br['from'], val_br['to'],
                pba, qba, pbb, qbb, pbc, qbc, baseZ,
                nbranch_ABC, nbranch_ABC * 1, nbus_ABC * 1)

            # Phase C
            pcc, qcc = -2 * z[2, 2][0], -2 * z[2, 2][1]
            pca, qca = -(- z[0, 2][0] + math.sqrt(3) * z[0, 2][1]), -(
                    - z[0, 2][1] - math.sqrt(3) * z[0, 2][0])
            pcb, qcb = -(- z[1, 2][0] - math.sqrt(3) * z[1, 2][1]), -(
                    - z[0, 2][1] + math.sqrt(3) * z[1, 2][0])
            A2, A1 = voltage_cons_pri(
                A2, A1,
                idx, val_br['from'], val_br['to'],
                pca, qca, pcb, qcb, pcc, qcc, baseZ,
                nbranch_ABC, nbranch_ABC * 2, nbus_ABC * 2)

        idx += 1

    ########################################################################################################################
    # va,vb,vc are phase voltage vectors for non slack buses
    # v0a, v0b, v0c are phase voltages for slack bus
    # fpa, fpb, fpc are real power flow vectors for each phase
    # fqa, fqb, fqc are reactive power flow vectors for each phase
    # pa, pb, pc are real power injections at non slack buses
    # qa, qb, qc are reactive power injections at non slack buses

    # We write the vector expressions
    # vn = [va.T, vb.T, vc.T]
    # v0 = [v0a, v0b, v0c]
    # v = [v0.T vn.T]
    # fp = [fpa.T, fpb.T, fpc.T]
    # fq = [fqa.T, fqb.T, fqc.T]
    # f = [fp.T, fq.T]
    # pn = [pa.T, pb.T, pc.T]
    # qn = [qa.T, qb.T, qc.T]
    # p = [pn.T, qn.T]

    # The power balance equations are
    # p_all = A @ f
    # Remove rows of A corresponding to the slack nodes to get the square matrix H22 (for a radial system)
    # p = H22 @ f
    # f = H22_inv @ p

    # The lindistflow equations are
    # v_delta = -A2 @ f
    # where v_delta is the vector of voltage differences along the branches
    # v_delta = A1 @ v = [A1_slack A1_red] @ [v0.T vn.T] = (A1_slack @ v0) + (A1_red @ vn)
    # - A2 @ f = (A1_slack @ v0) + (A1_red @ vn)
    # vn = -(Avr_inv @ Av0) @ v0 - (Avr_inv @ A2) @ f

    # Denote the following
    # H11 = -(Avr_inv @ Av0)
    # H12 = -(Avr_inv @ A2)
    ########################################################################################################################
    slack_node_idx = [slack_index, slack_index + nbus_ABC, slack_index + 2 * nbus_ABC]
    A1_slack = A1[:, slack_node_idx]
    A1_red = np.delete(A1, slack_node_idx, axis=1)
    A1r_inv = np.linalg.inv(A1_red)

    logger.debug(A1)
    logger.debug(A1_red)
    logger.debug(A1_slack)

    # get incidence matrix as transpose of Avr
    # A_inc is nbus by nbranch
    A_inc = A1_red.T

    # get the voltage relation with power injections
    H11 = - (A1r_inv @ A1_slack)
    H12 = - (A1r_inv @ A2)
    H22 = np.kron(np.eye(2,dtype=int),A_inc)
    
    H22_inv = np.linalg.inv(H22)
    H_linear = np.hstack((H11, H12 @ H22_inv))

    # add rows to the H_linear matrix for the slack bus voltages
    # since we will be using them as measurement in our estimation
    z = np.hstack((
        np.identity(len(slack_node_idx)), 
        np.zeros(shape=(len(slack_node_idx), H22.shape[0]))
        ))
    for i in range(len(slack_node_idx)):
        H_linear = np.insert(H_linear, slack_node_idx[i], z[i, :], axis=0)
    

    ###################### Forming the big H matrix #############################
    small_o = np.zeros((A_inc.shape[1], A_inc.shape[1]))
    small_I = np.identity(A_inc.shape[1])

    # voltage of all nodes as a function of slack node voltage, injections and PV generations
    H1 = np.hstack((
        H_linear,                                           # slack node voltage, injections
        np.zeros((H_linear.shape[0], A_inc.shape[1])),      # real PV generation
        np.zeros((H_linear.shape[0], A_inc.shape[1]))       # reactive PV generation
        ))
    
    # real load forecast as a function of slack node voltage, injections and PV generations
    H2 = np.hstack((
        np.zeros((A_inc.shape[1], len(slack_node_idx))),    # slack node voltage
        -small_I, small_o,                                  # real and reactive injections
        small_I, small_o                                    # real and reactive PV generation
        ))
    # reactive load forecast as a function of slack node voltage, injections and PV generations
    H3 = np.hstack((
        np.zeros((A_inc.shape[1], len(slack_node_idx))),    # slack node voltage
        small_o, -small_I,                                  # real and reactive injections
        small_o, small_I                                    # real and reactive PV generation
        ))

    ## Pinj and Qinj measurements as functions of slack node voltage, injections and PV generations:
    H4 = np.hstack((
        np.zeros((A_inc.shape[1], len(slack_node_idx))), 
        small_I, small_o, 
        small_o, small_o 
        ))
    H5 = np.hstack((
        np.zeros((A_inc.shape[1], len(slack_node_idx))), 
        small_o, small_I, 
        small_o, small_o
        ))

    # stack them all to get the big H matrix
    H = np.vstack((H1, H2, H3, H4, H5))

    return H


def get_pq(
        bus_info: dict,
        source_bus: str,
        SBASE: float = 1e6
):
    n_bus = len(bus_info) - 1
    pq = np.zeros(shape=(6 * n_bus,))
    count = 0
    for keyb, val_bus in bus_info.items():
        if keyb != source_bus:
            # Real power injection at a bus
            # Phase A Real Power
            pq[count + n_bus * 0] = val_bus['pv'][0][0] + val_bus['pq'][0][0]
            # Phase B Real Power
            pq[count + n_bus * 1] = val_bus['pv'][1][0] + val_bus['pq'][1][0]
            # Phase C Real Power
            pq[count + n_bus * 2] = val_bus['pv'][2][0] + val_bus['pq'][2][0]

            # Phase A Reactive Power
            pq[count + n_bus * 3] = val_bus['pv'][0][1] + val_bus['pq'][0][1]
            # Phase B Reactive Power
            pq[count + n_bus * 4] = val_bus['pv'][1][1] + val_bus['pq'][1][1]
            # Phase C Reactive Power
            pq[count + n_bus * 5] = val_bus['pv'][2][1] + val_bus['pq'][2][1]

            count += 1
    return pq / (SBASE)

def get_pq_forecast(
        bus_info: dict,
        source_bus: str,
        SBASE: float = 1e6
):
    n_bus = len(bus_info) - 1
    pq = np.zeros(shape=(6 * n_bus,))
    count = 0
    for keyb, val_bus in bus_info.items():
        if keyb != source_bus:
            # Real power load at a bus
            # Phase A Real Power
            pq[count + n_bus * 0] = val_bus['pq'][0][0]
            # Phase B Real Power
            pq[count + n_bus * 1] = val_bus['pq'][1][0]
            # Phase C Real Power
            pq[count + n_bus * 2] = val_bus['pq'][2][0]

            # Reactive power load at a bus
            # Phase A Reactive Power
            pq[count + n_bus * 3] = val_bus['pq'][0][1]
            # Phase B Reactive Power
            pq[count + n_bus * 4] = val_bus['pq'][1][1]
            # Phase C Reactive Power
            pq[count + n_bus * 5] = val_bus['pq'][2][1]

            count += 1
    return pq / (SBASE)

def get_pv(
        bus_info: dict,
        source_bus: str,
        SBASE: float = 1e6
):
    n_bus = len(bus_info) - 1
    pv = np.zeros(shape=(6 * n_bus,))
    count = 0
    for keyb, val_bus in bus_info.items():
        if keyb != source_bus:
            # Real power generation at a bus
            # Phase A Real Power
            pv[count + n_bus * 0] = val_bus['pv'][0][0]
            # Phase B Real Power
            pv[count + n_bus * 1] = val_bus['pv'][1][0]
            # Phase C Real Power
            pv[count + n_bus * 2] = val_bus['pv'][2][0]

            # Reactive power generation at a bus
            # Phase A Reactive Power
            pv[count + n_bus * 3] = val_bus['pv'][0][1]
            # Phase B Reactive Power
            pv[count + n_bus * 4] = val_bus['pv'][1][1]
            # Phase C Reactive Power
            pv[count + n_bus * 5] = val_bus['pv'][2][1]

            count += 1
    return pv / (SBASE)


def get_v(
        bus_info: dict,
        source_bus: str
):
    n_bus = len(bus_info)
    v = np.zeros(shape=(3 * n_bus,))
    slack_index = []
    count = 0
    for keyb, val_bus in bus_info.items():
        if keyb == source_bus:
            slack_index = [count, count + n_bus, count + 2 * n_bus]
        v[count + n_bus * 0] = val_bus['vmag'][0]
        v[count + n_bus * 1] = val_bus['vmag'][1]
        v[count + n_bus * 2] = val_bus['vmag'][2]
        count += 1
    return v, slack_index


def get_vbase(
        bus_info: dict,
        voltages: VoltagesMagnitude
):
    n_bus = len(bus_info)
    v = np.zeros(shape=(3 * n_bus,))
    ids = voltages.ids
    vals = voltages.values

    for keyb, val_bus in bus_info.items():
        idx = val_bus["idx"]
        v[idx + n_bus * 0] = vals[ids.index(f"{keyb}.1")] if f"{keyb}.1" in ids else 1.0
        v[idx + n_bus * 1] = vals[ids.index(f"{keyb}.2")] if f"{keyb}.2" in ids else 1.0
        v[idx + n_bus * 2] = vals[ids.index(f"{keyb}.3")] if f"{keyb}.3" in ids else 1.0

    return v


def get_nodes(bus_info: dict) -> list:
    n_bus = len(bus_info)
    nodes = np.empty(shape=(3 * n_bus,), dtype='<U20')
    for keyb, val_bus in bus_info.items():
        idx = val_bus["idx"]
        nodes[idx + n_bus * 0] = f"{keyb}.1"
        nodes[idx + n_bus * 1] = f"{keyb}.2"
        nodes[idx + n_bus * 2] = f"{keyb}.3"
    return nodes


if __name__ == "__main__":

    # Inputs from the feeder federate saved in dictionary and other OEDISI data types
    bus_info = {
        'N2': {
            'idx': 0,
            'phases': ['1', '2', '3'],
            'kv': 7.199557856794634,
            'vmag': [7187.5, 7190.4, 7189.1],
            # 'vpu': [0.9870552869666603, 0.991655131792795, 0.9890323837566921],
            's_rated': 0.0,
            'pv': [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            'pq': [[-0.0, -0.0], [-0.0, -0.0], [-0.0, -0.0]],
            'eqid': 'N2'
        },

        'SOURCEBUS': {
            'idx': 1,
            'phases': ['1', '2', '3'],
            'kv': 7.199557856794634,
            'vmag': [7199.5, 7199.5, 7199.5],
            # 'vpu': [0.999971613586772, 0.9999739411968653, 0.9999725967771476],
            's_rated': 0.0,
            'pv': [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            # 'pq': [[2053709.2562983471, 1433148.6131437998], [1928820.098544116, 1308759.6483972222],
            #        [1986555.1984555677, 1390105.7640618484]],
            'pq': [[-0.0, -0.0], [-0.0, -0.0], [-0.0, -0.0]],
            'eqid': 'SOURCEBUS'
        },

        'N3': {
            'idx': 2,
            'phases': ['1', '2', '3'],
            'kv': 2.4017771198288433,
            'vmag': [ 2396.4, 2397.4, 2397.0],
            # 'vpu': [0.9357136937672939, 0.9444634479185576, 0.9392286647169615],
            's_rated': 0.0,
            'pv': [[266667, 0.0], [266667, 0.0], [266667, 0.0]],
            'pq': [[-0.0, -0.0], [-0.0, -0.0], [-0.0, -0.0]],
            'eqid': 'N3'
        },

        'N4': {
            'idx': 3,
            'phases': ['1', '2', '3'],
            'kv': 2.4017771198288433,
            'vmag': [2351.5, 2363.5, 2358.1],
            # 'vpu': [0.7984808956706196, 0.8581130396352007, 0.8247328404307817],
            's_rated': 0.0,
            'pv': [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
            'pq': [[-540000, -177489.31], [-540000, -177489.31],
                   [-540000, -177489.31]],
            'eqid': 'N4'
        }
    }

    branch_info = {
        'N2_SOURCEBUS': {
            'idx': 0,
            'type': 'LINE',
            'from': 0, 'to': 1,
            'fr_bus': 'N2', 'to_bus': 'SOURCEBUS',
            'phases': ['1', '2', '3'],
            'zprim': [
                [[0.1733114243335522, 0.4083439729266439], [0.059068554341517476, 0.1900225938743636],
                 [0.05813478498489997, 0.14580229440346207]],
                [[0.05906855434151743, 0.1900225938743635], [0.17674939336989512, 0.39702951340790393],
                 [0.05984739444692716, 0.16046725905764]],
                [[0.05813478498489982, 0.1458022944034619], [0.059847394446926966, 0.16046725905763978],
                 [0.174796719325521, 0.40342875923927407]]
            ]
        },

        # 'N3_N2': {
        #     'idx': 1,
        #     'type': 'LINE',
        #     'from': 2, 'to': 0,
        #     'fr_bus': 'N3', 'to_bus': 'N2',
        #     'phases': ['1', '2', '3'],
        #     'zprim': [
        #         [[0.0864586666666667*3, 0.5187520000000002*3], [-0.0, 0.0], [-0.0, 0.0]],
        #         [[-0.0, 0.0], [0.0864586666666667*3, 0.5187520000000002*3], [-0.0, 0.0]],
        #         [[-0.0, 0.0], [-0.0, 0.0], [0.0864586666666667*3, 0.5187520000000002*3]]
        #     ]
        # },

        'N3_N2': {
            'idx': 1,
            'type': 'LINE',
            'from': 2, 'to': 0,
            'fr_bus': 'N3', 'to_bus': 'N2',
            'phases': ['1', '2', '3'],
            'zprim': [
                [[0.0864586666666667 * 0.001, 0.5187520000000002 * 0.001], [-0.0, 0.0], [-0.0, 0.0]],
                [[-0.0, 0.0], [0.0864586666666667 * 0.001, 0.5187520000000002 * 0.001], [-0.0, 0.0]],
                [[-0.0, 0.0], [-0.0, 0.0], [0.0864586666666667 * 0.001, 0.5187520000000002 * 0.001]]
            ]
        },

        'N4_N3': {
            'idx': 2,
            'type': 'LINE',
            'from': 3, 'to': 2,
            'fr_bus': 'N4', 'to_bus': 'N3',
            'phases': ['1', '2', '3'],
            'zprim': [
                [[0.21663928041694003, 0.5104299661583045], [0.07383569292689657, 0.23752824234295394],
                 [0.07266848123112463, 0.18225286800432733]],
                [[0.07383569292689661, 0.23752824234295405], [0.22093674171236846, 0.49628689175987967],
                 [0.07480924305865859, 0.20058407382204974]],
                [[0.07266848123112463, 0.1822528680043271], [0.07480924305865855, 0.2005840738220496],
                 [0.2184958991569009, 0.5042859490490924]]
            ]
        }
    }

    source_bus = "SOURCEBUS"

    baseV = VoltagesMagnitude(
        values=[7199.557856794634, 7199.557856794634, 7199.557856794634,
                7199.557856794634, 7199.557856794634, 7199.557856794634,
                2401.7771198288433, 2401.7771198288433, 2401.7771198288433,
                2401.7771198288433, 2401.7771198288433, 2401.7771198288433],
        ids=['SOURCEBUS.1', 'SOURCEBUS.2', 'SOURCEBUS.3',
             'N2.1', 'N2.2', 'N2.3',
             'N3.1', 'N3.2', 'N3.3',
             'N4.1', 'N4.2', 'N4.3']
    )

    SBASE = 1e6

    H = get_Hmat(bus_info, branch_info, source_bus, SBASE=SBASE)

    pq = get_pq(bus_info, source_bus, SBASE=SBASE)
    pq_load = get_pq_forecast(bus_info, source_bus, SBASE=SBASE)
    vmag, vslack = get_v(bus_info, source_bus)

    # compute per unit voltage magnitudes
    vbase = get_vbase(bus_info, baseV)
    vmag_pu = vmag / vbase
    v0 = vmag_pu[vslack]

    ### add noises:
    v_sigma = 0.005  # 5% * nominal values 1-2% (make it general based on Sbase/Vbase
    pl_sigma = 0.02  # 5-10% error (based on forecast uncertainty)
    ql_sigma = pl_sigma
    p_inj_sigma = 0.001 # 0.1-0.5% error as the customers are billed (>= 15 min interval)
    q_inj_sigma = p_inj_sigma # 1% of nominal KVA error probably?

    ############ estimation:
    V_W = np.array([1/(v_sigma**2)]*len(vmag_pu))
    Pl_W = np.array([1 / (pl_sigma ** 2)] * len(pq_load))
    Pinj_W = np.array([1 / (p_inj_sigma ** 2)] * len(pq))
    W_array = np.hstack((V_W, Pl_W, Pinj_W))
    W = np.diag(W_array)
    # W = W/1e4
    # W = W @ np.linalg.inv(W)

    Z_meas = np.hstack((vmag_pu, pq_load, pq))
    G = H.T @ W @ H
    G_inv = np.linalg.inv(G)
    x_est = G_inv @ H.T @ W @ Z_meas

    v_sub_est = np.sqrt(x_est[:len(vslack)])
    p_inj_est = x_est[len(vslack): len(vslack)+int(len(pq)/2)]
    q_inj_est = x_est[len(vslack)+ int(len(pq)/2): len(vslack) + len(pq)]
    Ppv_inj_est = x_est[len(vslack) + len(pq): len(vslack) + len(pq) +int(len(pq_load)/2) ]
    Qpv_inj_est = x_est[len(vslack) + len(pq) +int(len(pq_load)/2): ]

    print([f"{p:0.4f}" for p in Ppv_inj_est])
    print([f"{q:0.4f}" for q in Qpv_inj_est])


    # ## Make it Automatic::
    # ## Checking if thats correct:
    # pv_inj = 0*pq
    # pv_inj[1:8:3] = 266667/SBASE #0.267
    # x_check_all = np.concatenate((v0, pq, pv_inj))
    # print(pv_inj)

    # measurement_all = H @ x_check_all

    

    # v_lin_all = measurement_all[:H_linear.shape[0]]
    # pl_lin_all = measurement_all[H_linear.shape[0]: H_linear.shape[0] + A_inc.shape[1]]
    # ql_lin_all = measurement_all[H_linear.shape[0] + A_inc.shape[1]: H_linear.shape[0] + (2*A_inc.shape[1])]
    # p_inj_all = measurement_all[H_linear.shape[0] + (2*A_inc.shape[1]): H_linear.shape[0] + (3*A_inc.shape[1])]
    # q_inj_all = measurement_all[H_linear.shape[0] + (3*A_inc.shape[1]): ]



    
    

    

    # # v_meas = add_noise(v_lin_all[nodes_ord], mean=0, std=v_sigma)
    # v_meas = add_noise(v_lin_all, mean=0, std=v_sigma)
    # pl_meas = add_noise(pl_lin_all, mean=0, std=pl_sigma)
    # ql_meas = add_noise(ql_lin_all, mean=0, std=ql_sigma)
    # p_inj_meas = add_noise(p_inj_all, mean=0, std=ql_sigma)
    # q_inj_meas = add_noise(q_inj_all, mean=0, std=ql_sigma)

    # ############ inv calculation:
    # Z_meas = np.hstack((v_meas, pl_meas, ql_meas, p_inj_meas, q_inj_meas))
    # H_inv = np.linalg.pinv(H)
    # x_est_noise = H_inv @ Z_meas

    # v_sub_est = np.sqrt(x_est_noise[:len(vslack)])
    # p_inj_est = x_est_noise[len(vslack): len(vslack)+(1*A_inc.shape[1])]
    # q_inj_est = x_est_noise[len(vslack)+ (1 * A_inc.shape[1]): len(vslack) + (2 * A_inc.shape[1])]
    # Ppv_inj_est = x_est_noise[len(vslack)+(2*A_inc.shape[1]): len(vslack)+(3*A_inc.shape[1]) ]
    # Qpv_inj_est = x_est_noise[len(vslack)+(3*A_inc.shape[1]): ]

    # print(Ppv_inj_est)
    # #######



    # ############ estimation:
    # V_W = np.array([1/(v_sigma**2)]*(A_inc.shape[1]+len(vslack)))
    # Pl_W = np.array([1 / (pl_sigma ** 2)] * A_inc.shape[1] * 2)
    # Pinj_W = np.array([1 / (p_inj_sigma ** 2)] * A_inc.shape[1] * 2)
    # W_array = np.hstack((V_W, Pl_W, Pinj_W))
    # W = np.diag(W_array)
    # # W = W/1e4
    # # W = W @ np.linalg.inv(W)

    # Z_meas = np.hstack((v_meas, pl_meas, ql_meas, p_inj_meas, q_inj_meas))
    # G = H.T @ W @ H
    # G_inv = np.linalg.inv(G)
    # x_est = G_inv @ H.T @ W @ Z_meas

    

    # print(Ppv_inj_est)
    ###### PV Detected:

    
    
    
    
    
    ##########################################################################################
    ########################### Plot to visualize the results ################################
    ##########################################################################################
    # x_check = np.concatenate((v0, pq))
    # v_linear = H_linear @ x_check

    # import matplotlib.pyplot as plt

    # # get the node list to sort the outputs in plotting order
    # nodes = get_nodes(bus_info)
    # node_ids = baseV.ids
    # nodes_ord = [nodes.tolist().index(nd) for nd in node_ids]

    # v_true = vmag_pu[nodes_ord]
    # v_est = np.sqrt(v_linear[nodes_ord])

    # # Check to see if the LDF performs well for the system
    # fig, ax = plt.subplots(1, 1, figsize=(20, 12))
    # ax.plot(range(len(v_true)), v_true, 'b--', lw=2.0,
    #         marker='*', markersize=20, label='true')
    # ax.plot(range(len(v_est)), v_est, color='crimson',
    #         marker='*', markersize=20,
    #         ls='dashed', lw=2.0, label='estimated')
    # ax.set_xticks(list(range(len(v_true))), nodes[nodes_ord], fontsize=15, rotation=30)
    # ax.tick_params(axis='y', labelsize=20)
    # ax.set_ylabel("Voltage (in p.u.)", fontsize=20)
    # ax.set_xlabel("Nodes in network", fontsize=20)
    # ax.legend(fontsize=25, markerscale=2)
    # fig.suptitle("Linearized voltages obtained from LinDistFlow", fontsize=30)
    # fig.savefig("check_voltage_estimatev3.png")
    ##########################################################################################
    ##########################################################################################