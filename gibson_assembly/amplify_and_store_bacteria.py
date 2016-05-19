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

experiment_name = "sfgfp_puc19_gibson_storage_v1"
utils.experiment_name = experiment_name 

def amplify_and_store_bacteria(source_bacteria_well):
    # Tubes and plates
    growth_plate = p.ref(expid("growth"), cont_type="96-deep",   storage="cold_80",  discard=False)
    growth_wells = growth_plate.wells(['A1','A2'])    
    
    prepare_growth_wells(growth_wells)
    
    p.distribute(source_bacteria_well,growth_wells,ul(25),mix_before=True,mix_vol=ul(150),
                 allow_carryover=True)
    
    for growth_well in growth_wells:
        p.mix(growth_well, volume=ul(450), repetitions=10)    
        
    p.cover(growth_plate)
    #grow bacteria until they are in their log phase of growth
    #https://www.qiagen.com/us/resources/technologies/plasmid-resource-center/growth%20of%20bacterial%20cultures/
    p.incubate(growth_plate, "warm_37", "{}:hour".format(15), shaking=True)
    
    p.uncover(growth_plate)
    
    #add glycerol
    p.provision(p.inv['glycerol'],growth_wells,ul(500))
    
    for growth_well in growth_wells:
        p.mix(growth_well, volume=ul(900), repetitions=10)
    

def prepare_growth_wells(growth_wells):
    #
    # To LB, add ampicillin at ~1/1000 concentration
    #
    p.provision(p.inv["LB Miller"], growth_wells, ul(474.5))
    p.provision(p.inv["Amp 100mgml"], growth_wells, ul(0.5))
    

def measure_growth_wells(growth_plate,growth_wells):
    #
    # Growth: absorbance and fluorescence over 24 hours
    # Absorbance at 600nm: cell growth
    # Fluorescence at 485nm/510nm: sfGFP
    #
    hr = 12
    for t in range(0,hr*2+1,hr):
        if t > 0:
            p.cover(growth_plate)
            p.incubate(growth_plate, "warm_37", "{}:hour".format(hr), shaking=True)
            p.uncover(growth_plate)

    
if __name__ == '__main__':
    # -----------------------------------------------------------------------
    # Inventory
    #
    inv = {
        "src_bacteria_plate"  : "ct18zdk7bp43ew", # inventory; Original source of bacteria
    }
    
    #highest flourescence reading well
    BEST_WELL_ID = 'B7'
    
    if "--test" in sys.argv:
        test_inv =  {
            "src_bacteria_plate"  : "ct18ybg4qb5kzr", # inventory; Original source of bacteria
        }
        inv.update(test_inv)
    
    #Source of bacteria
    source_bacteria_well = p.ref('src_bacteria_plate', id=inv['src_bacteria_plate'], cont_type="96-flat", storage="cold_4").well(BEST_WELL_ID)   
    
    utils.init_inventory_well(source_bacteria_well)
    
    amplify_and_store_bacteria(source_bacteria_well)
    
    # ---------------------------------------------------------------
    # Output protocol
    #
    jprotocol = json.dumps(p.as_dict(), indent=2)
    print(jprotocol)    