from material_functions import MaterialProperties
from chain_in_fix_matrix import FixMatrix
from chain_in_variable_matrix import VariableMatrix

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

monomer_weight = 28.05
entanglement_weight = 1250
equilibrium_time = 5.7e-09
weight_average = 331800


ds = 5 * 10**(-3)
xgrid = np.arange(0, 1, ds)
tgrid = np.logspace(-9, 5, 512)
relaxation_modulus = 1873058.63232
distribution = pd.read_excel(
    'data/inputs/mte144257_simulation_inputs.xlsx', usecols=[0, 3])


def diffusion_in_fix_matrix(distribution, xgrid, tgrid, ds):
    """Calculate the dilusion factor,ie the average of the survival probability is returned and used as a constraint in the variable matrix calculation
    """
    weight_fractions, weight_distributions = distribution.Wi, distribution.M
    av_probs = 0
    for chain in range(len(distribution)):
        # fix matrix calculations
        chain_instance = FixMatrix(weight_distributions[chain], entanglement_weight,
                                   monomer_weight, equilibrium_time, xgrid, tgrid, ds)

        survival_prob = chain_instance.pde_solver()
        av_probs += weight_fractions[chain] * survival_prob

    dilution_factor = chain_instance.dilusionFactor(av_probs)

    return dilution_factor


dilusion_factor = diffusion_in_fix_matrix(distribution, xgrid, tgrid, ds)


def diffusion_in_variable_matrix(distribution, xgrid, tgrid, ds, dilusion_factor, relaxation_modulus):
    weight_fractions, weight_distributions = distribution.Wi, distribution.M
    av_probs = 0
    rouse_contributions = 0
    for chain in range(len(distribution)):
        chain_instance = VariableMatrix(weight_distributions[chain], entanglement_weight,
                                        monomer_weight, equilibrium_time, xgrid, tgrid, ds, dilusion_factor, relaxation_modulus)
        survival_prob = chain_instance.pde_solver()
        av_probs += weight_fractions[chain] * survival_prob
        rouses = chain_instance.rouseRelaxation_Gr()
        rouse_contributions += rouses
    dilution_factor = chain_instance.dilusionFactor(av_probs)
    chains_probabilities = av_probs * dilution_factor
    plt.plot(tgrid, av_probs, 'k')
    plt.plot(tgrid, chains_probabilities, 'r:', alpha=0.6)
    plt.semilogx()
    plt.show()
    return chains_probabilities


# sol = diffusion_in_variable_matrix(
#     distribution, xgrid, tgrid, ds, dilusion_factor, relaxation_modulus)

pd = pd.DataFrame({'t': tgrid, 'p_new': dilusion_factor})
pd.to_excel('undated.xlsx')
