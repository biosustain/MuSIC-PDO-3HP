__author__ = 'kaizhuang'

####################
# Global Variables #
####################
from framed.bioreactor import AEROBIC, ANAEROBIC

# reaction names in metabolic models
r_glc = 'R_EX_glc_e'
r_glyc = 'R_EX_glyc_e'
r_pdo = 'R_EX_1_3_pdo_e'
r_3hp = 'R_EX_3_hp_e'
r_o2 = 'R_EX_o2_e'

conditions = [AEROBIC, ANAEROBIC]
substrates = [r_glc, r_glyc]
targets = [r_pdo, r_3hp]
