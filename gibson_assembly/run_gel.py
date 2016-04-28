""" 

Derived from Brian Naughton @ http://blog.booleanbiotech.com/genetic_engineering_pipeline_python.html

"""
import sys
import json
from custom_protocol import CustomProtocol as Protocol
from utils import (ul, expid, init_inventory_well, touchdown,
                   dead_volume)
import numpy


p = Protocol()

# ---------------------------------------------------
# Set up experiment
#
experiment_name = "sfgfp_pcr_v1"
template_length = 726

_options = {
    'run_gel'        : True,  # run a gel to see the plasmid size
}
options = {k for k,v in _options.items() if v is True}

# ---------------------------------------------------
# Inventory and provisioning
# https://developers.transcriptic.com/v1.0/docs/containers
#

#@ToDo: convert sfgfp from 10uM to 100uM
#@ToDo: confirm my 1ng/uL is correct

inv = {
    'sfgfp_pcr_v1_amplified':            'ct18xam4qn2ydg', # inventory: unknown concentration
    'water':                             'rs17gmh5wafm5p', # catalog; Autoclaved MilliQ H2O
}

if "--test" in sys.argv:
    test_inv = {
        'sfgfp_pcr_v1_amplified':        'ct18x9rfdw99w5', # inventory: unknown concentration
    }
    inv.update(test_inv)



#provision the pcr product
pcr_product_well = p.ref("sfgfp_pcr_v1_amplified", id=inv['sfgfp_pcr_v1_amplified'], cont_type="micro-1.5", storage="cold_20").well(0)

# Temporary tubes for use, then discarded
#diluted_product_well1 =     p.ref("diluted_1_10th_pcr_product_well", cont_type="micro-1.5", storage="cold_20").well(0)
diluted_product_well2 =     p.ref("diluted_1_20th_pcr_product_well", cont_type="micro-1.5", storage="cold_20").well(0)
#create a water tube with 36ul water and 4ul pcr product --> 40uL --> 1/10 dilution

#p.provision(inv["water"], diluted_product_well1, ul(36))
#p.transfer(pcr_product_well, diluted_product_well1, ul(4),mix_before=True,mix_after=True,mix_vol=ul(20))

p.provision(inv["water"], diluted_product_well2, ul(38))
p.transfer(pcr_product_well, diluted_product_well2, ul(2),mix_before=True,mix_after=True,mix_vol=ul(20))

# Initialize all existing inventory
all_inventory_wells = [pcr_product_well]
for well in all_inventory_wells:
    init_inventory_well(well)

# --------------------------------------------------------
# Run a gel
#

#p.gel_separate([diluted_product_well1,diluted_product_well2],
p.gel_separate([diluted_product_well2],
               ul(20), "agarose(10,1.2%)", "ladder1", "10:minute", expid("gel", experiment_name))
print(json.dumps(p.as_dict(), indent=2))
