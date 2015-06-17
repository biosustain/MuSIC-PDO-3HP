

## Helper Functions
def calculateGWP(CO2=0, CH4=0, N2O=0):
    """
    Calculate Global Warming Potential (CO2 equivalents)
    """
    CO2e = CO2 + 21.0 * CH4 + 310.0 * N2O
    return CO2e


def calculateEP(PO4=0, NOx=0, N=0, P=0):
    """
    Calculate Eutrophication Potential (PO4 equivalents)
    """
    PO4e = PO4 + NOx * 0.13 + N * 0.42 + P * 3.06
    return PO4e


    