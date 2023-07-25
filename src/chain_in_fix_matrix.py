# importing the base class
from material_functions import MaterialProperties
import numpy as np
from scipy.sparse import diags


class FixMatrix(MaterialProperties):
    """This class considers the chain in a fix matrix, ie is only the chains that moves but its surrounding does not.
    * Its base class sets all the material properties for the polymer.
    * We will employ symmetry and solve the
    """
    c = 2

    def __init__(self, weight_average, entanglement_weight, monomer_weight, equilibrium_time, xgrid, tgrid, ds):

        MaterialProperties.__init__(
            self, weight_average, entanglement_weight, monomer_weight, equilibrium_time)
        self.full_xgrid = xgrid  # grid nodes along chain.
        # Solving only half due to symmetry for fluctuation paths
        self.half_xgrid = xgrid[:int(len(xgrid) / 2)]
        self.tgrid = tgrid  # time values over which to solve the pde
        self.ns = len(self.full_xgrid)  # length of the grid along the chain
        self.nt = len(self.tgrid)
        self.ic = self.initial_conditions()
        self.ds = ds

    def early_fluctuation(self):
        """cal early primitive path fluctuations

        returns:
            list of arrays primitive path fluctuation for chain i
        Notes:
            should be inserted in a loop, looping over all chains"""

        # solving over the positive half ie from 0 to 0.5
        return (9 * np.pi**3 * self.equilibrium_time * self.Z**4 * self.half_xgrid**4) / 16

    def late_fluctuation(self):
        """calculate the late fluctuations of the primitive path
        args:
            x(array): grid nodes along chain
            Z(float): number of entanglement of chain i
            tauR(float): rouse rotational time
            phi(list): is the dilution factor
        returns:
            tlate(array): primitive path fluctuation for chain i
        Notes:
            should be inserted in a loop, looping over all chains"""

        U = ((3 * self.Z * self.half_xgrid**2) / (2 * self.c)).astype('float64')
        A = self.rouse_relaxation_time / (self.c**2)  # c**2=4
        return A * np.exp(U)

    @staticmethod
    def pathIntersections(early_path, late_path):
        """find intersection between the two paths
        returns
            interception(list): index of where curves intercept as a list
            """

        idx = np.argwhere(np.diff(np.sign(early_path - late_path))).flatten()
        return idx

    def contourLengthFluctuations(self):
        """cal the fluctuation the chain contour length
        calls the following functions: early_path(array), late_path(array),intersections(array)
        returns:
            ppfs(list(array)): primitive path fluctuations
        Notes:
            to be used as input when solving the pde"""

        early_path = self.early_fluctuation()
        late_path = self.late_fluctuation()
        intersections = self.pathIntersections(early_path, late_path)

        # if there are 2 intersection points
        if len(intersections) == 2:
            # when there are 2 intersection
            begining_section = early_path[:intersections[0]]
            # slicing to index b4 intersection
            # slicing 4rm 1 index after intersection
            ending_section = late_path[intersections[1] + 1:]
            # pulling values from the midsection using their index
            early_middle = early_path[intersections[0]:intersections[1] + 1]
            late_middle = late_path[intersections[0]:intersections[1] + 1]
            # update path for the midsection
            av_paths = np.sqrt(early_middle * late_middle)
            # actual fluctuation path
            paths = np.concatenate(
                (begining_section, av_paths, ending_section))
        # for cases where you just have a single point of intersection between the early and late path
        elif len(intersections) == 1:
            # cases where a single intersection occurs and it happens at the early section of the chain
            # condition where you have a C1 intersection but not C2
            if late_path[-1] < early_path[-1]:
                # values before intersections
                early_middle = early_path[intersections[0]:]
                late_middle = late_path[intersections[0]:]
                av_paths = np.sqrt(early_middle * late_middle)
                paths = np.concatenate(
                    (early_path[:intersections[0]], av_paths))
            else:  # if the single intersection occurs towards the end of the chain
                early_portion = early_path[:intersections[0] + 1]
                late_portion = late_path[:intersections[0] + 1]
                ending_section = late_path[intersections[0] + 1:]
                av_paths = np.sqrt(early_portion * late_portion)
                paths = np.concatenate((av_paths, ending_section))

        else:
            # when there is no intersection then use t_early
            paths = early_path
        # since symetric we calculated over half, mirrored the solution and concatenated
        ppfs = np.concatenate((paths, paths[::-1]))
        ppfs = ppfs[1:-1]  # excluding the values for the boundaries
        # list helpfull for situation where we have a time dependency
        return [ppfs]

    def initial_conditions(self):
        """initial solutions for the PDE at t=0"""
        P = np.ones(self.ns)
        P[0] = np.exp(-self.tgrid[0] / self.equilibrium_time)
        P[-1] = np.exp(-self.tgrid[0] / self.equilibrium_time)
        return P

    def settingMatrix(self, time_delta, fluctuation):
        """ 3 scenarios to consider
            * only values, mostly for cases when constraint release is not considered
            * case where the length of ppf is same as that of time
            * case where you have fewer ppfs than time data
            """

        d_coef = (time_delta) / (2 * self.reptation_time *
                                 (np.pi**2) * (self.ds**2))
        # cal the reaction term and turning it into a diagonal matrix

        # array due to fluctuation
        reactionTerm = np.diag(time_delta / (2 * fluctuation))

        # matrix is set-up during each iteration because dt changes so as to implement the time stepping in log
        a = diags([-d_coef, 1 + 2 * d_coef, -d_coef], [-1, 0, 1],
                  shape=(self.ns - 2, self.ns - 2)).toarray() + reactionTerm
        b1 = diags([d_coef, 1 - 2 * d_coef, d_coef], [-1, 0, 1],
                   shape=(self.ns - 2, self.ns - 2)).toarray() - reactionTerm
        return a, b1, d_coef

    def pde_solver(self):
        """"Solves the diffusion equation using Crank-Nicolson method
        * makes use of a log implementation
        * 2 helper functions are used: ic,contourlenghtfluctuations

        returns:
            survival_prob(array): chains survival probability at each time step
            Pout(array): solutions of the diffusion equations

        """

        survival_probability = []  # survival probability of a chain at each time step
        Pout = []  # to store the probability calculated at time step

        dt_values = np.diff(self.tgrid)
        fluctuations = self.contourLengthFluctuations()

        P = self.ic.copy()  # be storing the value of the probability calculated at the current step
        # storing the initial solutions
        Pout.append(P)
        survival_probability.append(np.trapz(y=self.ic, dx=self.ds))
        # looping over time
        for t in range(self.nt - 1):  # nt = dt+1

            # current delta t value
            dt = dt_values[t]
            if len(fluctuations) == 1:  # for the case of a fix matrix where there is no time dependency
                fluctuation = fluctuations[0]
            else:
                fluctuation = fluctuations[t]
            A, B1, diffusionTerm = self.settingMatrix(dt, fluctuation)

            # solving the diffusion equation
            Pn = P
            # multiply the B matrix by the initial values
            B = np.dot(Pn[1:-1], B1)  # excluding the boundaries

            # updating the boundaries at the next time step
            boundary_right = np.exp(-self.tgrid[t + 1] / self.equilibrium_time)
            boundary_left = np.exp(-self.tgrid[t + 1] / self.equilibrium_time)
            B[0] = B[0] + 2 * diffusionTerm * boundary_left
            B[-1] = B[-1] + 2 * diffusionTerm * boundary_right

            # solve the algebra to get values at new time point
            P[1:-1] = np.linalg.solve(A, B)

            # solution to that time step;array are mutable so we write a copy
            Pout.append(P.copy())
            # assert all(P < 0), 'negative values in probability, change ds step'

            survival_probability.append(np.trapz(y=P, dx=self.ds))

        return Pout, np.asarray(survival_probability)
