class MaterialProperties:

    def __init__(self, weight_average, entanglement_weight, monomer_weight, equilibrium_time):
        self.weight_average = weight_average  # Mw
        self.entanglement_weight = entanglement_weight  # Me
        self.monomer_weight = monomer_weight  # Mo
        self.equilibrium_time = equilibrium_time  # tau_e
        self.Z = self.N_en()
        # self.reptation_time = self.reptationTime()
        self.rouse_relaxation_time = self.rouseRelaxationTime()

    def N_en(self):
        """Cal the degree of entanglement (Z)
        Args:
            * M_w(float): molecular weight
            * M_e(float): molecular weight between entanglementequilibrium_time
        returns:equilibrium_time
            Nen(float): degree of entanglement"""
        return self.weight_average / self.entanglement_weight

    def degreePolymerization(self):
        """Cal the degree of Polymerization
        Args:
            mw(float): Molecular weight of chain
            mo(float): Monomer weight
        returns:
            Ni(float): degree of polymeriself.Zation"""
        return self.weight_average / self.monomer_weight

    def rouseRelaxationTime(self):
        """Cal the tau_r
        args:
            self.Z(float): array of floats with the number of entanglements
            t_e(float):  equilibrium time which is a material constant
        returns:
            tau_r(float): an array of relaxation time for chains"""

        return (self.Z**2) * self.equilibrium_time

    def reptationTime(self):
        """Cal the reptation time tau_d
        args:
            Z(float): array of floats with the number of entanglements
            t_e(float):  equilibrium time which is a material constant
        returns:
            tau_d(float): an array of relaxation time for chains"""
        return 3 * (self.Z**3) * self.equilibrium_time


monomer_weight = 28.05
entanglement_weight = 1250
equilibrium_time = 5.7e-09
weight_average = 331800


mat1 = MaterialProperties(
    weight_average, entanglement_weight, monomer_weight, equilibrium_time)
print(mat1.rouse_relaxation_time)
