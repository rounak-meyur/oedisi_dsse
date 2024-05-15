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
    # compute base impedance
    basekV = bus_info[source_bus]['kv'] / np.sqrt(3)
    baseZ = basekV ** 2 / SBASE
    

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
                    A1, k_frm_C, k_to_C, 
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
                    A1, k_frm_C, k_to_C, 
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
    # p_all = A1 @ f
    # Remove rows of A1 corresponding to the slack nodes to get the square matrix H22 (for a radial system)
    # p = H22 @ f
    # f = H22_inv @ p

    # The lindistflow equations are
    # v_delta = -A2 @ f
    # where v_delta is the vector of voltage differences along the branches
    # v_delta = Av @ v = [Av0 Avr] @ [v0.T vn.T] = (Av0 @ v0) + (Avr @ vn)
    # A2 @ f = (Av0 @ v0) + (Avr @ vn)
    # vn = -(Avr_inv @ Av0) @ v0 - (Avr_inv @ A2) @ f
    
    # Denote the following
    # H11 = -(Avr_inv @ Av0)
    # H12 = -(Avr_inv @ A2)
    ########################################################################################################################
    slack_node_idx = [slack_index, slack_index+nbus_ABC, slack_index+2*nbus_ABC]
    slack_node_idx_pq = slack_node_idx + [slack_index+n_bus, slack_index+n_bus+nbus_ABC, slack_index+n_bus+2*nbus_ABC]
    Av0 = Av[:,slack_node_idx]
    Avr = np.delete(Av, slack_node_idx, axis=1)
    Avr_inv = np.linalg.inv(Avr)

    H11 = - (Avr_inv @ Av0)
    H12 = - (Avr_inv @ A2)

    H22 = np.delete(A1, slack_node_idx_pq, axis=0)
    H22_inv = np.linalg.inv(H22)
    return H11, H12@H22_inv


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



if __name__ == "__main__":

    

    # Inputs from the feeder federate saved in dictionary and other OEDISI data types
    bus_info = {
        'N2': {
            'idx': 0, 
            'phases': ['1', '2', '3'], 
            'kv': 7.199557856794634, 
            'vmag': [7106.360895653739, 7139.477719183365, 7120.597377968156], 
            's_rated': 0.0, 
            'pv': [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], 
            'pq': [[-0.0, -0.0], [-0.0, -0.0], [-0.0, -0.0]], 
            'eqid': 'N2'
            }, 

        'SOURCEBUS': {
            'idx': 1, 
            'phases': ['1', '2', '3'], 
            'kv': 7.199557856794634, 
            'vmag': [7199.353486521364, 7199.370243417042, 7199.360567670321], 
            's_rated': 0.0, 
            'pv': [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], 
            'pq': [[2053709.2562983471, 1433148.6131437998], [1928820.098544116, 1308759.6483972222], [1986555.1984555677, 1390105.7640618484]], 
            'eqid': 'SOURCEBUS'
            }, 

        'N3': {
            'idx': 2, 
            'phases': ['1', '2', '3'], 
            'kv': 2.4017771198288433, 
            'vmag': [2247.3750297296033, 2268.389577854546, 2255.819609890638], 
            's_rated': 0.0, 
            'pv': [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], 
            'pq': [[-0.0, -0.0], [-0.0, -0.0], [-0.0, -0.0]], 
            'eqid': 'N3'
            }, 

        'N4': {
            'idx': 3, 
            'phases': ['1', '2', '3'], 
            'kv': 2.4017771198288433, 
            'vmag': [1917.769330972959, 2060.9928349426377, 1980.8314278404248], 
            's_rated': 0.0, 
            'pv': [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], 
            'pq': [[-1799941.343138584, -871605.0855056401], [-1800124.9897248761, -871933.0649904874], [-1799889.7682501655, -871714.021517176]], 
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
                [[0.1733114243335522, 0.4083439729266439], [0.059068554341517476, 0.1900225938743636], [0.05813478498489997, 0.14580229440346207]], 
                [[0.05906855434151743, 0.1900225938743635], [0.17674939336989512, 0.39702951340790393], [0.05984739444692716, 0.16046725905764]], 
                [[0.05813478498489982, 0.1458022944034619], [0.059847394446926966, 0.16046725905763978], [0.174796719325521, 0.40342875923927407]]
                ]
            }, 
        
        'N3_N2': {
            'idx': 1, 
            'type': 'LINE', 
            'from': 2, 'to': 0, 
            'fr_bus': 'N3', 'to_bus': 'N2', 
            'phases': ['1', '2', '3'], 
            'zprim': [
                [[0.0864586666666667, 0.5187520000000002], [-0.0, 0.0], [-0.0, 0.0]], 
                [[-0.0, 0.0], [0.0864586666666667, 0.5187520000000002], [-0.0, 0.0]], 
                [[-0.0, 0.0], [-0.0, 0.0], [0.0864586666666667, 0.5187520000000002]]
                ]
            }, 
        
        'N4_N3': {
            'idx': 2, 
            'type': 'LINE', 
            'from': 3, 'to': 2, 
            'fr_bus': 'N4', 'to_bus': 'N3', 
            'phases': ['1', '2', '3'], 
            'zprim': [
                [[0.21663928041694003, 0.5104299661583045], [0.07383569292689657, 0.23752824234295394], [0.07266848123112463, 0.18225286800432733]], 
                [[0.07383569292689661, 0.23752824234295405], [0.22093674171236846, 0.49628689175987967], [0.07480924305865859, 0.20058407382204974]], 
                [[0.07266848123112463, 0.1822528680043271], [0.07480924305865855, 0.2005840738220496], [0.2184958991569009, 0.5042859490490924]]
                ]
            }
        }

    source_bus = "SOURCEBUS"

    baseV = VoltagesMagnitude(
        values = [7199.557856794634, 7199.557856794634, 7199.557856794634, 
                7199.557856794634, 7199.557856794634, 7199.557856794634, 
                2401.7771198288433, 2401.7771198288433, 2401.7771198288433, 
                2401.7771198288433, 2401.7771198288433, 2401.7771198288433], 
        ids= ['SOURCEBUS.1', 'SOURCEBUS.2', 'SOURCEBUS.3', 
            'N2.1', 'N2.2', 'N2.3', 
            'N3.1', 'N3.2', 'N3.3', 
            'N4.1', 'N4.2', 'N4.3']
    )


    # Get the linear matrices
    Hv, Hpq = get_Hmat(bus_info, branch_info, source_bus)
    pq = get_pq(bus_info, source_bus, SBASE=100.0e6)
    vmag, vslack = get_v(bus_info, source_bus)
    
    # compute per unit voltage magnitudes
    vbase = get_vbase(bus_info, baseV)
    vmag_pu = vmag / vbase
    
    # select rows corresponding to voltages and columns 
    # corresponsing to slack bus voltage and node injections
    v0 = vmag_pu[vslack]
    H_check = np.hstack((Hv,Hpq))
    z = np.hstack((np.identity(len(vslack)), np.zeros(shape=(len(vslack),Hpq.shape[1]))))
    for i in range(len(vslack)):
        H_check = np.insert(H_check, vslack[i], z[i,:], axis=0)
    x_check = np.concatenate((v0, pq))
    v_linear = H_check @ x_check


    ######### Plot to visualize the results #################
    import matplotlib.pyplot as plt
    # get the node list to sort the outputs in plotting order
    nodes = get_nodes(bus_info)
    node_ids = baseV.ids
    nodes_ord = [nodes.tolist().index(nd) for nd in node_ids]

    v_true = vmag_pu[nodes_ord]
    v_est = v_linear[nodes_ord]
    
    # Plot the linear estimate versus true value
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