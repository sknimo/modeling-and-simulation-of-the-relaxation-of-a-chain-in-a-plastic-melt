from chain_in_fix_matrix import FixMatrix
import numpy as np


class VariableMatrix(FixMatrix):
    """this class accounts for both motions of the chain and its environment
    There is a chain of inheritance here MaterialProperties and FixMatrix"""

    def __init__(self, weight_average, entanglement_weight, monomer_weight, equilibrium_time, xgrid, tgrid, ds, dilusion_factor, relaxation_modulus):
        FixMatrix.__init__(self, weight_average, entanglement_weight,
                           monomer_weight, equilibrium_time, xgrid, tgrid, ds)
        self.dilustion_factor = dilusion_factor
        self.relaxation_modulus = relaxation_modulus

    def late_fluctuation(self):
        """calculate the late fluctuations of the primitive path

        returns:
            tlate(list): time dependent primitive path fluctuation for chain i
        Notes:
            should be inserted in a loop, looping over all chains"""

        phi = self.dilustion_factor.reshape(-1, 1)
        U = ((3 * self.Z * phi * self.half_xgrid**2) /
             (2 * self.c)).astype('float64')
        A = self.rouse_relaxation_time / (self.c**2)
        return A * np.exp(U)

    def rouseRelaxation_Gr(self):
        """Calculate the high frequency rouse processes of single chain
        Args:
            self.Z(float): Molecular weight between entanglement
            N_i(float): degree of polymerization
            Gn(float): the plateau modulus
        returns:
            Gr(array): high frequency rouse
        """

        def gr_slowModes(t):
            """Matrix for the first part of Gr"""
            output1 = []
            # because loop stops one point before end value
            nen_limit = int(self.Z) + 1
            for k in range(1, nen_limit):
                expo = -(2 * k**2 * t) / self.rouse_relaxation_time
                output1.append(np.exp(expo))
            return np.sum(output1, axis=0) / 3

        def gr_fastModes(t):
            """Matrix for the second part of Gr"""
            upper_limit = int(self.degree_polymerization) + 1
            lower_limit = int(self.Z) + 1
            output2 = []

            for k in range(lower_limit, upper_limit):
                expo = -(2 * k**2 * t) / self.rouse_relaxation_time
                output2.append(np.exp(expo))
            return np.sum(output2, axis=0)

        # the created functions are only used if fraction is not zero
        slowModes = self.relaxation_modulus * gr_slowModes(self.tgrid)
        fastModes = self.relaxation_modulus * gr_fastModes(self.tgrid)
        result = (slowModes + fastModes) / self.Z
        return result
