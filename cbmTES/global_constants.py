__author__ = 'kaizhuang'

####################
# Global Variables #
####################
from framed.bioreactor import AEROBIC, ANAEROBIC
from collections import OrderedDict

# reaction names in metabolic models
r_glc = 'R_EX_glc_e'
r_glyc = 'R_EX_glyc_e'
r_pdo = 'R_EX_1_3_pdo_e'
r_3hp = 'R_EX_3_hp_e'
r_o2 = 'R_EX_o2_e'
r_3HPP = 'R_EX_3HPP_e'
r_13PDO = 'R_EX_13PDO_e'

#conditions = [AEROBIC, ANAEROBIC]
#substrates = [r_glc, r_glyc]
#targets = [r_pdo, r_3hp]

fba_constraints = OrderedDict()

fba_constraints['ecoli', r_glc, AEROBIC] = {r_glc: (-10, 0), r_glyc: (0, 0), r_o2: (-15, None)}
fba_constraints['ecoli', r_glc, ANAEROBIC] = {r_glc: (-20, 0), r_glyc: (0, 0), r_o2: (-0, None)}
fba_constraints['ecoli', r_glyc, AEROBIC] = {r_glc: (0, 0), r_glyc: (-15, 0),r_o2: (-15, None)}
fba_constraints['ecoli', r_glyc, ANAEROBIC] = {r_glc: (0, 0), r_glyc: (-15, 0),r_o2: (0, None)}

# B. Sonnleitner 1986
fba_constraints['scere', r_glc, AEROBIC] = {r_glc: (-20, 0), r_glyc: (0, 0), r_o2: (-8, None)}
fba_constraints['scere', r_glc, ANAEROBIC] = {r_glc: (-20, 0), r_glyc: (0, 0), r_o2: (-0.05, None)}
# Ocha-Estopier 2010
fba_constraints['scere', r_glyc, AEROBIC] = {r_glc: (0, 0), r_glyc: (-4.5, 0),r_o2: (-8, None)}