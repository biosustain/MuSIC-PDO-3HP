from framed.bioreactor.bioreactors import IdealFedbatch
from framed.bioreactor import Organism
from global_constants import *

# kinetic parameters
ecoli_parameters ={}
ecoli_parameters['Vmax_glc_aerobic'] = -10
ecoli_parameters['Vmax_glc_anaerobic'] = -20
ecoli_parameters['Km_glc']= 0.015

ecoli_parameters['Vmax_o2'] = -15

ecoli_parameters['Vmax_glyc'] = -15
ecoli_parameters['Km_glyc']= 0.1 # mM, assumed

scere_parameters = {}
scere_parameters['Vmax_glc'] = -20
scere_parameters['Km_glc'] = -0.55

scere_parameters['Vmax_glyc'] = -4.5

scere_parameters['Km_glyc']= 0.1 # mM, assumed
scere_parameters['Vmax_o2'] = -8
                    
# Defining custom E. coli class
class Ecoli(Organism):

    def __init__(self, model, id=None, fba_objective=None, fba_constraints={}, model_deepcopy=True, titer_limit={}):
        super(Ecoli, self).__init__(model, id, fba_objective, fba_constraints, model_deepcopy)
        self.titer_limit = titer_limit
            
    def update(self):
        BR = self.environment
        
        # determine if the product titer is limiting the organism.  if it does, shut down the organism
        titer_limited = False
        for rxn, concentration in self.titer_limit.items():
            if rxn in BR.metabolites:
                rid = BR.metabolites.index(rxn)
                if BR.S[rid] > concentration:
                    titer_limited = True
                    #print 'titer limited'
                        
        # aerobic uptake rates
        if BR.oxygen_availability == AEROBIC:  
            # glucose uptake
            if r_glc in BR.metabolites and not titer_limited:
                rid = BR.metabolites.index(r_glc)
                vlb_glc = float(ecoli_parameters['Vmax_glc_aerobic'] * BR.S[rid] / (BR.S[rid] + ecoli_parameters['Km_glc']))
                self.fba_constraints[r_glc] = (vlb_glc, 0)
            else:
                self.fba_constraints[r_glc] = (0, 0)
                
            # glycerol uptake
            if r_glyc in BR.metabolites and not titer_limited:
                rid = BR.metabolites.index(r_glyc)
                vlb_glyc = float(ecoli_parameters['Vmax_glyc'] * BR.S[rid] / (BR.S[rid] + ecoli_parameters['Km_glyc']))
                self.fba_constraints[r_glyc] = (vlb_glyc, 0)
            else:
                self.fba_constraints[r_glyc] = (0, 0)
            
            # oxygen uptake
            self.fba_constraints[r_o2] = (ecoli_parameters['Vmax_o2'], None)
            
        # anaerobic uptake rates
        if BR.oxygen_availability == ANAEROBIC:  
            # glucose uptake
            if r_glc in BR.metabolites and not titer_limited:
                rid = BR.metabolites.index(r_glc)
                vlb_glc = float(ecoli_parameters['Vmax_glc_anaerobic'] * BR.S[rid] / (BR.S[rid] + ecoli_parameters['Km_glc']))
                self.fba_constraints[r_glc] = (vlb_glc, 0)
            else:
                self.fba_constraints[r_glc] = (0, 0)
                
            # glycerol uptake 
            if r_glyc in BR.metabolites and not titer_limited:
                rid = BR.metabolites.index(r_glyc)
                vlb_glyc = 0
                self.fba_constraints[r_glyc] = (vlb_glyc, 0)
            else:
                self.fba_constraints[r_glyc] = (0, 0)
                
            # oxygen uptake
            self.fba_constraints[r_o2] = (0, None)

# Defining custom S. cerevisiae class
class Scerevisiae(Organism):

    def __init__(self, model, id=None, fba_objective=None, fba_constraints={}, model_deepcopy=True, titer_limit={}):
        super(Scerevisiae, self).__init__(model, id, fba_objective, fba_constraints, model_deepcopy)
        self.titer_limit = titer_limit
            
    def update(self):
        BR = self.environment
        
        # determine if the product titer is limiting the organism.  if it does, shut down the organism
        titer_limited = False
        for rxn, concentration in self.titer_limit.items():
            if rxn in BR.metabolites:
                rid = BR.metabolites.index(rxn)
                if BR.S[rid] > concentration:
                    titer_limited = True
                    #print 'titer limited'
                        

        # glucose uptake
        if r_glc in BR.metabolites and not titer_limited:
            rid = BR.metabolites.index(r_glc)
            vlb_glc = float(scere_parameters['Vmax_glc'] * BR.S[rid] / (BR.S[rid] + scere_parameters['Km_glc']))
            self.fba_constraints[r_glc] = (vlb_glc, 0)
        else:
            self.fba_constraints[r_glc] = (0, 0)
                
        # glycerol uptake
        if r_glyc in BR.metabolites and not titer_limited:
            rid = BR.metabolites.index(r_glyc)
            vlb_glyc = float(scere_parameters['Vmax_glyc']  * BR.S[rid] / (BR.S[rid] + scere_parameters['Km_glyc']))
            self.fba_constraints[r_glyc] = (vlb_glyc, 0)
        else:
            self.fba_constraints[r_glyc] = (0, 0)
            
        # oxygen uptake
        if BR.oxygen_availability == AEROBIC:
            self.fba_constraints[r_o2] = (-8, None)
        elif BR.oxygen_availability == ANAEROBIC:  
            self.fba_constraints[r_o2] = (-0.05, None)




def convert_mol2gram(metric, A, B, chemical_mass_table):
    # for Bioprocess A->B, convert the base unit of strain performance metric from mol to gram.
    # Product Yield(s): mol/mol -> g/g (x molar_mass_B / molar_mass_A)
    # Biomass Yield: g/mol -> g/g (x 1/molar_mass_A)
    # Productivity: mol/L/hr -> g/L/hr (x molar_mass_B)
    # Titer: mol/L -> g/L (x molar_mass_B)
    
    MM_A = chemical_mass_table[A]  # molar mass of the substrate
    MM_B = chemical_mass_table[B]  # molar mass of the product
    
    converter = {'yield': MM_B / MM_A,
                 'productivity': MM_B,
                 'titer': MM_B,
                 'yield_biomass': 1/MM_A            
                 }

    for c_id, c_factor in converter.items():
        metric[c_id] = [m * c_factor for m in metric[c_id]]

        
chemical_mass_table = {
r_glc:  180.16,  # g/mol
 r_glyc: 92.09,   # g/mol
 r_pdo:  76.09,   # g/mol
 r_3hp:  90.08   # g/mol
 }
 
        