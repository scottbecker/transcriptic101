"""Grow and store bacteria

From protocol: https://www.addgene.org/plasmid-protocols/create-glycerol-stock/

"""

import sys
import json
from custom_protocol import CustomProtocol as Protocol
from autoprotocol.protocol import Container, Ref
import utils
from utils import (ul, expid, init_inventory_well, touchdown,
                   dead_volume)

p = Protocol()

options = {}
for k, v in list(options.items()):
    if v is False: del options[k]
    

def measure_flourescence_and_absorption(plate,wells):
    # Absorbance at 600nm: cell growth
    # Fluorescence at 485nm/510nm: sfGFP

    p.fluorescence(plate, wells,
                   excitation="485:nanometer", emission="535:nanometer",
                   dataref=expid("fl2"), flashes=25)
    p.absorbance(plate, wells,
                 wavelength="600:nanometer",
                 dataref=expid("abs"), flashes=25)
    

    
if __name__ == '__main__':
    experiment_name = "sfgfp_puc19_alt_clone_measure"
    utils.experiment_name = experiment_name     
    
    # -----------------------------------------------------------------------
    # Inventory
    #
    
    inv = {
        "growth_plate"  : "ct1926ecqpeeb3", # inventory; Original source of bacteria
    }
    
    #highest flourescence reading well
    BEST_WELL_ID = 'B7'
    
    if "--test" in sys.argv:
        test_inv =  {
            "growth_plate"  : "ct18zub3vb6px6", # inventory; Original source of bacteria
        }
        inv.update(test_inv)
    
    #Source of bacteria
    growth_plate = p.ref('growth_plate', id=inv['growth_plate'], cont_type="96-flat", storage="cold_4") 
    
    growth_wells = growth_plate.wells(["{}{}".format(row,col) for col in range(1,3) for row in "ABCDEFGH"])
    
    measure_flourescence_and_absorption(growth_plate,growth_wells)
    
    # ---------------------------------------------------------------
    # Output protocol
    #
    jprotocol = json.dumps(p.as_dict(), indent=2)
    print(jprotocol)    