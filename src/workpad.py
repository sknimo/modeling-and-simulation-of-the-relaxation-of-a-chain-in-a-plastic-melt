from material_functions import MaterialProperties
from chain_in_fix_matrix import FixMatrix
from chain_in_variable_matrix import VariableMatrix

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

monomer_weight = 28.05
entanglement_weight = 1250
equilibrium_time = 5.7e-09
#weight_average = 331800


ds = 5 * 10**(-3)
xgrid = np.arange(0, 1, ds)
tgrid = np.logspace(-9, 5, 100)
relaxation_modulus = 1873058.63232
dataset = pd.read_csv(
    '../data/inputs/Dataset.csv')
distribution = pd.DataFrame(
    {'M': 10**dataset.iloc[:, 0], 'Wi': dataset.iloc[:, 1]})


def diffusion_in_fix_matrix(distribution, xgrid, tgrid, ds):
    """Calculate the dilusion factor,ie the average of the survival probability is returned and used as a constraint in the variable matrix calculation
    Args:
        distribution(df): contains the molecular weight distribution of the polymer
        xgrid(array): the nodes along the chain onto which the solution of the PDE will be sought
        tgrid(array): the time steps along which the solution would be sought, ie the evolution of the diffusion over time.
        ds(float): spacing between the nodes on the chain
    returns:
        dilution_factor(array): the mean field approximation of the survival of a test chain, with the constraint that it doesn't relax faster than t^(1/2)
    Note:
        the fix matrix considers that it is only the chain that can move and its environment is fix
    """

    # unpacking the molecular weight distribution which dictates how the melt will relax when afforded the chance
    weight_fractions, weight_distributions = distribution.Wi, distribution.M
    av_probs = 0
    for chain in range(len(distribution)):
        # instantiating our class for every chain in the distribution
        chain_instance = FixMatrix(weight_distributions[chain], entanglement_weight,
                                   monomer_weight, equilibrium_time, xgrid, tgrid, ds)
        # the survival probability is the integral of the area under the curve of the solution of the PDE
        survival_prob = chain_instance.pde_solver()
        # averaging the survival probabilities
        av_probs += weight_fractions[chain] * survival_prob
        #======================Picking one chain and plot the output as a means of checking==========================#
        if chain == 50:
            check_prob = pd.DataFrame(
                {'time': tgrid, 'survival_prob': survival_prob})
            check_prob.to_excel(
                f'data/survival_probability_chain_{int(weight_distributions[chain])}.xlsx', index=False)
            early_p = chain_instance.early_fluctuation()[1:]
            late_p = chain_instance.late_fluctuation()[0][1:]
            ppf = chain_instance.contourLengthFluctuations()[0]
            ppf = ppf[:int(len(ppf) / 2)]
            check_path = pd.DataFrame(
                {'time': chain_instance.half_xgrid[1:], 'early_ppf': early_p, 'late_ppf': late_p, 'ppf': ppf})
            check_path.to_excel(
                f'data/ppf_{int(weight_distributions[chain])}.xlsx', index=False)
        #===========================================================================================================#
    # ensuring that the relaxation is not faster than allowed
    dilution_factor = chain_instance.dilusionFactor(av_probs)

    return dilution_factor


dilusion_factor = diffusion_in_fix_matrix(distribution, xgrid, tgrid, ds)


def diffusion_in_variable_matrix(distribution, xgrid, tgrid, ds, dilusion_factor, relaxation_modulus):
    """ Evaluates the diffusion of the chain when both the chain and its surrounding are in motion
    Args:
        distribution(df): contains the molecular weight distribution of the polymer
        xgrid(array): the nodes along the chain onto which the solution of the PDE will be sought
        tgrid(array): the time steps along which the solution would be sought, ie the evolution of the diffusion over time.
        ds(float): spacing between the nodes on the chain
        dilusion_factor(array): controls the contour length fluctuation of the chain
        relaxation_modulus(float): This is a material constant defined by the ratio of the stress to the strain of the material
    returns:
        chains_probability(array): the chain survival probability when both the chain and the environment are in motion
    """

    # unpacking the molecular weight distribution which dictates how the melt will relax when afforded the chance
    weight_fractions, weight_distributions = distribution.Wi, distribution.M
    av_probs = 0
    rouse_contributions = 0
    for chain in range(len(distribution)):
        chain_instance = VariableMatrix(weight_distributions[chain], entanglement_weight,
                                        monomer_weight, equilibrium_time, xgrid, tgrid, ds, dilusion_factor, relaxation_modulus)
        survival_prob = chain_instance.pde_solver()
        av_probs += weight_fractions[chain] * survival_prob
        # high frequencies contribution, not so relevant for industrial process though
        rouses = chain_instance.rouseRelaxation_Gr()
        rouse_contributions += weight_fractions[chain] * rouses
    # ensuring the relaxation is not faster than allowed
    dilution_factor = chain_instance.dilusionFactor(av_probs)
    # calculating the survival probability
    chains_probabilities = av_probs * dilution_factor
    relaxation_curve = chains_probabilities * relaxation_modulus
    total_relaxation = relaxation_curve + rouse_contributions

    # plt.plot(tgrid, relaxation_curve, 'r-', alpha=0.6)
    # plt.plot(tgrid, total_relaxation, 'b-.', alpha=0.6)
    # plt.loglog()
    # plt.xlabel('time[s]')
    # plt.ylabel('G(t)[Pa]')
    # plt.savefig('check.png')
    # plt.show()
    # plt.close()
    return total_relaxation


stress_relaxation = diffusion_in_variable_matrix(
    distribution, xgrid, tgrid, ds, dilusion_factor, relaxation_modulus)

df_relaxation = pd.DataFrame(
    {'time[s]': tgrid, 'relaxation modulus [Pa]': stress_relaxation})
df_relaxation.to_excel('../data/outputs/stress_relaxation.xlsx', index=False)
