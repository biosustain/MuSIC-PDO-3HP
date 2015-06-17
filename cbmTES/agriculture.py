__author__ = 'kaizhuang'

import core
import energy

### Agricultural Sector Reaction Stoichiometries ###
stoichiometry = {}
reversibility = {}


### Soy Production ###

# Farming
    # Ref: Landis et al. 2007
stoichiometry['soy_farming'] = {}
reversibility['soy_farming'] = False
stoichiometry['soy_farming']['soy'] = 1 # kg
stoichiometry['soy_farming']['land'] = -1.0/2700 #hecta
#stoichiometry['soy_farming']['co2e'] = calculateGWP(CO2=157, CH4=0.19, N2O=0.5) #g CO2 equivalents
stoichiometry['soy_farming']['po4e'] = core.calculateEP(P=0.22, NOx=10.8)/1000 #kg PO4 equivalents
stoichiometry['soy_farming']['energy'] = -2.13 #MJ, does not include labour energy cost

# Crushing Plant
    # Ref: Sheehan et al. 1998
stoichiometry['soy_crushing']={}
reversibility['soy_crushing'] = False
stoichiometry['soy_crushing']['soy'] = -1 #kg
stoichiometry['soy_crushing']['soyoil'] = 0.18 #kg
stoichiometry['soy_crushing']['soymeal'] = 0.574 #kg
#stoichiometry['soy_crushing']['water'] = 0.16 # kg
stoichiometry['soy_crushing']['waste_cellulosic' ] = 0.074#kg
    # Ref: Hill et al. 2006
stoichiometry['soy_crushing']['energy'] = - 3.6*0.027/energy.eff_electricity - 2.44*0.260/energy.eff_steam #MJ


# Conversion Plant: Degumming and Transesterification
    # Ref: Sheehan et al. 1998
stoichiometry['soy_conversion']={}
reversibility['soy_conversion'] = False
stoichiometry['soy_conversion']['soyoil'] = -1 #kg
stoichiometry['soy_conversion']['biodiesel'] = 0.96 #kg
stoichiometry['soy_conversion']['glyc'] = 0.205 #kg
    # Ref: Hill et al. 2006
stoichiometry['soy_conversion']['energy']= - 3.6*0.024/energy.eff_electricity - 2.44*0.395/energy.eff_steam #MJ

### Corn Production ###
# Ref: Landis et al. 2007
stoichiometry['corn_farming'] = {}
reversibility['corn_farming'] = False
stoichiometry['corn_farming']['corn'] = 1 # kg
stoichiometry['corn_farming']['land'] = -1.0/8550 #hecta
#stoichiometry['corn_farming']['co2e'] = calculateGWP(CO2=175, CH4=0.32, N2O=1.54) #g CO2 equivalents
stoichiometry['corn_farming']['po4e'] = core.calculateEP(P=0.17, NOx=19.9)/1000 #kg PO4 equivalents
stoichiometry['corn_farming']['energy'] = -2.75 #MJ

# Ref: Urban et al. 2009
stoichiometry['corn_wetmill']={}
reversibility['corn_wetmill'] = False
stoichiometry['corn_wetmill']['corn'] = -1 #kg
stoichiometry['corn_wetmill']['starch'] = 1.84/3.38 #kg
stoichiometry['corn_wetmill']['waste_cellulosic' ] = 1 - stoichiometry['corn_wetmill']['starch'] # kg
stoichiometry['corn_wetmill']['energy'] = - 0.292/3.38*3.6/energy.eff_electricity - 0.088/3.38*53.6 # MJ

# Ref: Urban et al. 2009
stoichiometry['corn_hydrolysis']={}
reversibility['corn_hydrolysis'] = False
stoichiometry['corn_hydrolysis']['starch'] = -1 #kg
stoichiometry['corn_hydrolysis']['glc'] = 0.95 #kg
stoichiometry['corn_hydrolysis']['waste_organic'] = 1 - stoichiometry['corn_hydrolysis']['glc']  #kg

