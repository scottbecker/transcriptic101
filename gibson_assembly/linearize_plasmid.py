"""Protocol for cutting pUC19 with EcoRI and hindiii

Derived from Brian Naughton @ http://blog.booleanbiotech.com/genetic_engineering_pipeline_python.html

"""

import sys
import json
from custom_protocol import CustomProtocol as Protocol
from utils import (ul, expid, init_inventory_well, dead_volume,
                   assert_valid_volume)

p = Protocol()

experiment_name = "puc19_ecori_hindiii_v1"

options = {}

inv = {
    'water':    "rs17gmh5wafm5p",   # catalog; Autoclaved MilliQ H2O; ambient
    "reagent_plate":    "ct18xpszvkrxdb",   # inventory; 10 units each enzyme / ul
}

if "--test" in sys.argv:
    test_inv = {
        "reagent_plate":    "ct18xkermammek",   # inventory; 10 units each enzyme / ul
    }
    inv.update(test_inv)


reagent_plate = p.ref("re_reagent_plate", id=inv['reagent_plate'], 
                      cont_type="96-pcr", storage="cold_20")

ecori_hindiii_well = reagent_plate.wells(["A1"])[0]
cutsmart_well = reagent_plate.wells(["B1"])[0]
pUC19_well = reagent_plate.wells(["H12"])[0]

# Initialize all existing inventory
all_inventory_wells = [ecori_hindiii_well,cutsmart_well,pUC19_well]
for well in all_inventory_wells:
    init_inventory_well(well)

# Tubes and plates we use and then discard
water_tube = p.ref("water_tube", cont_type="micro-1.5", discard=True).well(0)
pcr_plate  = p.ref("pcr_plate",  cont_type="96-pcr", discard=True)

# The result of the experiment, a pUC19 cut by EcoRI and HindIII, goes in this tube for storage
puc19_cut_tube  = p.ref(expid("puc19_cut",experiment_name), cont_type="micro-1.5", storage="cold_20").well(0)

# -------------------------------------------------------------
# Provisioning and diluting.
# Diluted EcoRI can be used more than once
#
p.provision(inv["water"], water_tube, ul(500))

# -------------------------------------------------------------
# Restriction enzyme cutting pUC19
# 3 experiments (1 without re)
# 50ul total reaction volume for cutting 1ug of DNA:
# 42uL water
# 5ul CutSmart 10x
# 1ul pUC19 1ml/ml or 1ug/ul (1ug of DNA)
# 2ul EcoRI+hindiii (20 units each, >10 units per ug DNA)

p.distribute(water_tube, pcr_plate.wells(["A1","B1","A2"]), ul(42))
p.distribute(pUC19_well, pcr_plate.wells(["A1","B1","A2"]), ul(1), 
             mix_before=True, mix_after=True, mix_vol=ul(5))
p.distribute(cutsmart_well, pcr_plate.wells(["A1","B1","A2"]), ul(5), 
             mix_before=True, mix_after=True, mix_vol=ul(10))
p.distribute(water_tube, pcr_plate.wells(["A2"]), ul(2),ignore_mix_after_warning=True)
p.distribute(ecori_hindiii_well, pcr_plate.wells(["A1","B1"]), ul(2), 
             mix_before=True, mix_after=True, mix_vol=ul(5))
assert all(well.volume == ul(50) for well in pcr_plate.wells(["A1","B1","A2"]))

p.mix(pcr_plate.wells(["A1","B1","A2"]), volume=ul(25), repetitions=10)

# Incubation to induce cut, then heat inactivation of EcoRI+hindiii
p.seal(pcr_plate)
#NEB enzymes are time saver so they only need 15minutes
p.incubate(pcr_plate, "warm_37", "15:minute", shaking=False)
#inactivate RE's (see this site for thermal temps)
#https://www.neb.com/tools-and-resources/usage-guidelines/heat-inactivation
p.thermocycle(pcr_plate, [{"cycles":  1, "steps": [{"temperature": "80:celsius", "duration": "20:minute"}]}], volume=ul(50))

# --------------------------------------------------------------
# Gel electrophoresis, to ensure the cutting worked
#
p.unseal(pcr_plate)
p.mix(pcr_plate.wells(["A1","B1","A2"]), volume=ul(25), repetitions=5)
p.transfer(water_tube, pcr_plate.wells(["D1","E1","D2"]), ul(15))
p.transfer(pcr_plate.wells(["A1","B1","A2"]), pcr_plate.wells(["D1","E1","D2"]), ul(8), mix_after=True, mix_vol=ul(10))
assert all(well.volume == ul(20) + dead_volume["96-pcr"] for well in pcr_plate.wells(["D1","E1","D2"]))

p.gel_separate(pcr_plate.wells(["D1","E1","D2"]), ul(20), "agarose(10,2%)", "ladder2", "15:minute", expid("gel",experiment_name))


# ----------------------------------------------------------------------------
# Then consolidate all cut plasmid to one tube (puc19_cut_tube).
#
remaining_volumes = [well.volume - dead_volume['96-pcr'] for well in pcr_plate.wells(["A1","B1"])]
p.consolidate(pcr_plate.wells(["A1","B1"]), puc19_cut_tube, remaining_volumes, allow_carryover=True)

#verify that we haven't gone below what can actually be pulled from tubes
assert_valid_volume([water_tube, ecori_hindiii_well, cutsmart_well, pUC19_well])

# ---------------------------------------------------------------
# Test protocol
#
print(json.dumps(p.as_dict(), indent=2))