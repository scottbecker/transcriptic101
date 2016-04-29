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

_options = {
    'run_absorbance'        : True,
}
options = {k for k,v in _options.items() if v is True}

# ---------------------------------------------------
# Inventory and provisioning
# https://developers.transcriptic.com/v1.0/docs/containers
#

inv = {
    'pcr_reagent_plate':     'ct18x95j499zx4', # inventory: A1: Q5 polymerase,
                                                           # A2: Buffer, A3: Enhancer, B1: dNTP 10mM
                                                           # C1: exosap
    'water':                 'rs17gmh5wafm5p', # catalog; Autoclaved MilliQ H2O
    'dna_to_clean':          'ct18xam4qn2ydg', # inventory: sfgfp_pcr_v1_amplified
}

if "--test" in sys.argv:
    test_inv = {
        'pcr_reagent_plate': 'ct18x92yfcbhhz', # inventory: A1: Q5 polymerase,
                                                              # A2: Buffer, A3: Enhancer, B1: dNTP 10mM
                                                              # C1: exosap
     'dna_to_clean':         'ct18x9rfdw99w5', # inventory: sfgfp_pcr_v1_amplified
    }
    inv.update(test_inv)


# Existing inventory
dna_to_clean_well = p.ref("dna_to_clean", id=inv['dna_to_clean'], cont_type="micro-1.5", storage="cold_20").well(0)
pcr_reagent_plate = p.ref("pcr_reagent_plate", id=inv['pcr_reagent_plate'], 
                          cont_type="96-pcr", storage="cold_20")

exosap_it_well = pcr_reagent_plate.wells(["C1"])[0]

# New inventory resulting from this experiment
exosap_reaction_plate = p.ref(expid("cleaned",experiment_name), cont_type="96-pcr", storage="cold_20")


if 'run_absorbance' in options:
    abs_plate = p.ref("abs_plate_%s_clean"%experiment_name, 
                      cont_type="96-flat", storage="cold_20", discard=False)

# Initialize all existing inventory
all_inventory_wells = [dna_to_clean_well]
for well in all_inventory_wells:
    init_inventory_well(well)

# -----------------------------------------------------
# ExoSAP-IT PCR product cleanyup
# http://media.affymetrix.com/support/technical/usb/brief_proto/78200B.pdf
# requires 2uL of exosap for every 5 ul of dna product


p.transfer(dna_to_clean_well, exosap_reaction_plate.wells(["A1"]), ul(15), mix_before=True, mix_vol=ul(20))
p.transfer(exosap_it_well, exosap_reaction_plate.wells(["A1"]), ul(6), mix_before=True, mix_after=True,mix_vol=ul(6))

# ---------------------------------------------------------
# Thermocycle
#


cycles = [{"cycles":  1, "steps": [{"temperature": "37:celsius", "duration": "15:minute"},
                                   {"temperature": "80:celsius", "duration": "15:minute"}]}]
p.seal(exosap_reaction_plate)
p.thermocycle(exosap_reaction_plate, cycles, volume=ul(21))
p.unseal(exosap_reaction_plate)


#---------------------------------------------------------
# Absorbance dilution series. Take 1ul out of the 25ul pcr plate wells
# Good overview here: http://bitesizebio.com/13501/dna-concentration-purity/
#
if 'run_absorbance' in options:
    abs_wells = ["A1","A2","A3"]

    p.provision(inv["water"], abs_plate.wells(abs_wells[0:2]), ul(10))
    p.provision(inv["water"], abs_plate.wells(abs_wells[2]), ul(9))

    p.transfer(exosap_reaction_plate.wells(["A1"]), abs_plate.wells(["A1"]), ul(1), mix_after=True, mix_vol=ul(5))
    p.transfer(abs_plate.wells(["A1"]), abs_plate.wells(["A2"]), ul(1), mix_after=True, mix_vol=ul(5))
    p.transfer(abs_plate.wells(["A2"]), abs_plate.wells(["A3"]), ul(1), mix_after=True, mix_vol=ul(5))

    for wavelength in [260, 280, 320]:
        p.absorbance(abs_plate, abs_plate.wells(abs_wells),
                     "{}:nanometer".format(wavelength), expid("abs_{}".format(wavelength), experiment_name),
                     flashes=25)

print(json.dumps(p.as_dict(), indent=2))
