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
    'pcr_reagent_plate':                 'ct18x95j499zx4', # inventory: A1: Q5 polymerase,
                                                           # A2: Buffer, A3: Enhancer, B1: dNTP 10mM
    'water':                             'rs17gmh5wafm5p', # catalog; Autoclaved MilliQ H2O
    'sfgfp_puc19_primer_forward_10uM':   'ct18x9wbfksyc3', # inventory; micro-1.5, cold_20, 10uM
    'sfgfp_puc19_primer_reverse_10uM':   'ct18x9wbfm4v8c', # inventory; micro-1.5, cold_20, 10uM
    'sfgfp_2nM':                         'ct18vs7cmjat2c', # inventory; sfGFP tube #1, micro-1.5, cold_20, 2.12nM, 1ng/uL
}

if "--test" in sys.argv:
    test_inv = {
        'pcr_reagent_plate':                'ct18x92yfcbhhz', # inventory: A1: Q5 polymerase,
                                                              # A2: Buffer, A3: Enhancer, B1: dNTP 10mM
        'sfgfp_puc19_primer_forward_10uM':  'ct18x626u9nvne', # inventory; micro-1.5, cold_20, 10uM
        'sfgfp_puc19_primer_reverse_10uM':  'ct18x626u9yshq', # inventory; micro-1.5, cold_20, 10uM
        'sfgfp_2nM':                        'ct18x62qg8km37', # inventory
    }
    inv.update(test_inv)


# Existing inventory
template_tube = p.ref("sfgfp_2nM", id=inv['sfgfp_2nM'], cont_type="micro-1.5", storage="cold_20").well(0)
primer_wells = [p.ref('sfgfp_puc19_primer_forward_10uM', id=inv['sfgfp_puc19_primer_forward_10uM'], 
                      cont_type="micro-1.5", storage="cold_20").well(0),
                p.ref('sfgfp_puc19_primer_reverse_10uM', id=inv['sfgfp_puc19_primer_reverse_10uM'], 
                      cont_type="micro-1.5", storage="cold_20").well(0)]

pcr_reagent_plate = p.ref("pcr_reagent_plate", id=inv['pcr_reagent_plate'], 
                          cont_type="96-pcr", storage="cold_20")

q5_poly_well = pcr_reagent_plate.wells(["A1"])[0]
q5_buffer_well = pcr_reagent_plate.wells(["A2"])[0]
dNTP_well = pcr_reagent_plate.wells(["B1"])[0]


# New inventory resulting from this experiment
sfgfp_pcroe_out_tube = p.ref(expid("amplified",experiment_name), cont_type="micro-1.5", storage="cold_20").well(0)

# Temporary tubes for use, then discarded
mastermix_well = p.ref("mastermix", cont_type="micro-1.5", storage="cold_4",  discard=True).well(0)
water_tube =     p.ref("water",     cont_type="micro-1.5", storage="ambient", discard=True).well(0)
pcr_plate =      p.ref("pcr_plate", cont_type="96-pcr",    storage="cold_4",  discard=True)
if 'run_absorbance' in options:
    abs_plate = p.ref("abs_plate", cont_type="96-flat", storage="cold_20", discard=False)

# Initialize all existing inventory
all_inventory_wells = [template_tube] + primer_wells
for well in all_inventory_wells:
    init_inventory_well(well)

# -----------------------------------------------------
# Provision water once, for general use
#
p.provision(inv["water"], water_tube, ul(500))

# -----------------------------------------------------
# Q5 PCR protocol
# www.neb.com/protocols/2013/12/13/pcr-using-q5-high-fidelity-dna-polymerase-m0491
#
# 25ul reaction (we will run it 4 times, totally 100uL of solution)
# -------------                      4rxn total
# Q5 reaction buffer      5    ul --> 20uL
# Q5 polymerase           0.25 ul --> 1uL
# 10mM dNTP               0.5  ul --> 2uL
# 10uM forward primer     1.25 ul --> 5uL
# 10uM reverse primer     1.25 ul --> 5uL
# 1ng Template (1 ng/ul) 1 ul --> 4uL (1ng/ul concentration) 
# water                   
# -------------------------------
# Sum                     9.25 ul --> 37uL
# water                   15.75 uL --> 63uL
#
#

# Mastermix tube will have 96ul of stuff, leaving space for 3x1ul aliquots of template
p.transfer(water_tube,                mastermix_well, ul(63))
p.transfer(q5_buffer_well, mastermix_well, ul(20), mix_before=True, mix_vol=ul(20), mix_after=True)
p.transfer(q5_poly_well, mastermix_well, ul(1), mix_before=True, mix_vol=ul(5), mix_after=True)
p.transfer(dNTP_well, mastermix_well, ul(2), mix_before=True, mix_vol=ul(5), mix_after=True)
p.transfer(primer_wells[0], mastermix_well, ul(5), mix_before=True, mix_vol=ul(10), mix_after=True)
p.transfer(primer_wells[1], mastermix_well, ul(5), mix_before=True, mix_vol=ul(10), mix_after=True)
p.mix(mastermix_well, volume="48:microliter", repetitions=10)

# Transfer mastermix to pcr_plate without template
p.transfer(mastermix_well, pcr_plate.wells(["A1","B1","C1"]), ul(24))
#
p.transfer(mastermix_well, pcr_plate.wells(["A2"]),           ul(24)) # acknowledged dead volume problems
p.mix(pcr_plate.wells(["A1","B1","C1","A2"]), volume=ul(12), repetitions=10)

# Finally add template
p.transfer(template_tube,  pcr_plate.wells(["A1","B1","C1"]), ul(1), mix_after=True, mix_vol=ul(5))
p.mix(pcr_plate.wells(["A1","B1","C1"]), volume=ul(12.5), repetitions=10)

# ---------------------------------------------------------
# Thermocycle with Q5 and hot start
# 61.1 annealing temperature is recommended by NEB protocol
# p.seal is enforced by transcriptic
#
extension_time = int(max(2, numpy.ceil(template_length * (11.0/1000))))
assert 0 < extension_time < 60, "extension time should be reasonable for PCR"

cycles = [{"cycles":  1, "steps": [{"temperature": "98:celsius", "duration": "30:second"}]}] + \
    touchdown(70, 61, [8, 25, extension_time], stepsize=0.5) + \
    [{"cycles": 16, "steps": [{"temperature": "98:celsius", "duration": "8:second"},
                              {"temperature": "61.1:celsius", "duration": "25:second"},
                              {"temperature": "72:celsius", "duration": "{:d}:second".format(extension_time)}]},
     {"cycles":  1, "steps": [{"temperature": "72:celsius", "duration": "2:minute"}]}]
p.seal(pcr_plate)
p.thermocycle(pcr_plate, cycles, volume=ul(25))
p.unseal(pcr_plate)
# --------------------------------------------------------
# Run a gel to hopefully see a 740bp fragment
#
if 'run_gel' in options:
    p.mix(pcr_plate.wells(["A1","B1","C1","A2"]), volume=ul(12.5), repetitions=10)
    p.transfer(water_tube, pcr_plate.wells(["D1","E1","F1","D2"]),
               [ul(18),ul(16),ul(12),ul(12)])
    p.transfer(pcr_plate.wells(["A1","B1","C1","A2"]), pcr_plate.wells(["D1","E1","F1","D2"]),
               [ul(2), ul(4), ul(8), ul(8)], mix_after=True, mix_vol=ul(10))    
    p.gel_separate(pcr_plate.wells(["D1","E1","F1","D2"]),
                   ul(20), "agarose(10,2%)", "ladder1", "10:minute", expid("gel", experiment_name))


# -------------------------------------------------------------------------
# Then consolidate to one tube. Leave at least 3ul dead volume in each tube
#
remaining_volumes = [well.volume - dead_volume['96-pcr'] for well in pcr_plate.wells(["A1","B1","C1"])]
p.consolidate(pcr_plate.wells(["A1","B1","C1"]), sfgfp_pcroe_out_tube, remaining_volumes, allow_carryover=True)

print(json.dumps(p.as_dict(), indent=2))
