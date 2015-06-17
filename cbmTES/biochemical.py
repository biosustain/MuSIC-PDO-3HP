import energy
from framed.core.models import Metabolite, Reaction, Compartment, ConstraintBasedModel

# Energy Intensity of Bioprocesses

Ei = {}   # energy intensities of bioprocesses
# fermentation
Ei['sterilization'] = 0.244/energy.eff_steam # MJ/L of media
Ei['agitation'] = 0.0018/energy.eff_electricity #MJ/L/hr 
Ei['aeration'] = 0.009/energy.eff_electricity #MJ/L/hr, applicable only in aerobic environments

# separation of biomass
Ei['centrifugation_bacteria'] = 0.0252/energy.eff_electricity # MJ/L
Ei['centrifugation_yeast'] = 0.0054/energy.eff_electricity # MJ/L

# purification of product
Ei['evaporation_single_stage'] = 1.2*2.44/energy.eff_steam + 0.04*3.6/energy.eff_electricity   # MJ/L
Ei['evaporation_multi_stage'] = 0.3*2.44/energy.eff_steam + 0.005*3.6/energy.eff_electricity  # MJ/L


####################################################################
# function for preprocessing the results of bioreactor simulations #
####################################################################
def get_strains_metric(metric):    
    number_of_strains = len(metric['yield'])
    strains = []
    for i in range(number_of_strains):
        if (metric['growth'][i] > 0) and (metric['yield'][i] > 0):  # if not non-growth or non-producing strains
            s={}
            for m in metric:
                s[m] = metric[m][i]
            strains.append(s)        
    return strains

#################################################
# function for adding a bioprocess to the model #
#################################################

def add_bioprocess_reaction(model, A, B, strain, fermentation, separation, purification, purification_efficiency, neutralization = [], strain_id=[], process_energy_intensity=[], compartment=[]):
    # Add a bioprocess A + E -> B + X + W to the model
        # A: Substrate
        # B: Product
        # X: Biomass
        # W: Waste

    # Process names:
        # Fermentation
            # aeroFB: aerobic fedbatch
            # anaeFB: anaerobic fedbatch
        # Separation
            # ctrf: centrifugation
        # Purification
            # evap_s: evaporation (single-stage)
            # evap_m: evaporation (multi-stage)

    # Neutralization
        # the the product is an organic acid, then this is necessary for maintaining pH during fermentation
            # 3hp_lime: lime is used to neutralize 3-hydroxyprorionic acid
            
            
    Y_B = strain['product_yield'] #kg/kg
    Y_X = strain['biomass_yield'] #kg/kg
    T = strain['product_titer']/1000.0 #kg/L
    P = strain['productivity']/1000.0 #kg/L/hr
    R_B = purification_efficiency
    if process_energy_intensity:
        Ei = process_energy_intensity
    if not strain_id:
        strain_id = strain['strain_id']
        
    # calculating energy cost yield of the process (MJ of Energy / g of A)
    if fermentation == 'aeroFB':
        Y_E_fermentation = Ei['sterilization']*Y_B/T + Ei['agitation']*Y_B/P + Ei['aeration']*Y_B/P
    elif fermentation == 'anaeFB':
        Y_E_fermentation = Ei['sterilization']*Y_B/T + Ei['agitation']*Y_B/P
        
    if separation == 'ctrf':
        if 'ecoli' in strain_id:
            Y_E_separation = Ei['centrifugation_bacteria']*Y_B/T
        if 'scere' in strain_id:
            Y_E_separation = Ei['centrifugation_yeast']*Y_B/T
            
    if purification == 'evapS':
        Y_E_purification = Ei['evaporation_single_stage']*Y_B/T
    elif purification == 'evapM':
        Y_E_purification = Ei['evaporation_multi_stage']*Y_B/T
        
        
    Y_E_overall = Y_E_fermentation + Y_E_separation + Y_E_purification

    # calculating the overall waste yield and product yield (after purfication)
    Y_W = 1 - Y_B - Y_X
    Y_B_overall = Y_B * R_B
    Y_W_overall = Y_W + (1-R_B)*Y_B
    
    
    rxn = Reaction('BPR_' +A +'_' +B +'_' + strain_id +'_'  +fermentation +'_' +separation +'_' +purification , reversible=False)
    model.add_reaction(rxn)
    model.stoichiometry[(A+'_'+compartment,rxn.id)]= -1
    model.stoichiometry[('energy'+'_'+compartment,rxn.id)]= -Y_E_overall
    model.stoichiometry[(B+'_'+compartment,rxn.id)]= Y_B_overall
    model.stoichiometry[('biomass'+'_'+compartment,rxn.id)]= Y_X
    model.stoichiometry[('waste_organic'+'_'+compartment,rxn.id)]= Y_W_overall

    # Neutralization
    # this is a hack.  this is only for 3HP neutralization with lime.
        # lime is used to neutralize 3HP during fermentation
        # sulfuric_acid is used to recover 3HP after fermentation
            # Y_B is used for calculation, not Y_B_overall.
            # This is because all organic acids need to be neutralized, not just the purified product

    if neutralization == '3hp_lime':
        model.stoichiometry['lime'+'_'+compartment, rxn.id] = - 0.4113 * Y_B
        model.stoichiometry['sulfuric_acid'+'_'+compartment, rxn.id] = - 0.5444 * Y_B
        model.stoichiometry['gypsum'+'_'+compartment, rxn.id] = 0.7557 * Y_B

