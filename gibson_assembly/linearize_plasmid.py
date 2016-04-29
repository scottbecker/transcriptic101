"""Protocol for cutting pUC19 with EcoRI and BamHI

Derived from Brian Naughton @ http://blog.booleanbiotech.com/genetic_engineering_pipeline_python.html

"""

import sys
import json
from custom_protocol import CustomProtocol as Protocol
from utils import (ul, expid, init_inventory_well, dead_volume)

p = Protocol()

experiment_name = "puc19_ecori_bamhi_v1"

options = {}

inv = {
    'water':    "rs17gmh5wafm5p",   # catalog; Autoclaved MilliQ H2O; ambient
    "pUC19":    "rs17tcqmncjfsh",   # catalog; pUC19; cold_20
    "reagent_plate":    "",   # inventory; 10 units each enzyme / ul
}

if "--test" in sys.argv:
    test_inv = {
        "reagent_plate":    "ct18xjsc7f5x6r",   # inventory; 10 units each enzyme / ul
    }
    inv.update(test_inv)

# Tubes and plates I use then discard
re_tube    = p.ref("re_tube",    cont_type="micro-1.5", storage="cold_4", discard=True).well(0)
water_tube = p.ref("water_tube", cont_type="micro-1.5", storage="cold_4", discard=True).well(0)
pcr_plate  = p.ref("pcr_plate",  cont_type="96-pcr",    storage="cold_4", discard=True)

# The result of the experiment, a pUC19 cut by EcoRI, goes in this tube for storage
puc19_cut_tube  = p.ref(expid("puc19_cut"), cont_type="micro-1.5", storage="cold_20").well(0)

# -------------------------------------------------------------
# Restriction enzyme cutting pUC19
# 3 experiments (1 without re)
# 50ul total reaction volume for cutting 1ug of DNA:
# 5ul CutSmart 10x
# 1ul pUC19 1ml/ml or 1ug/ul (1ug of DNA)
# 2ul EcoRI+BamHI (20 units each, >10 units per ug DNA)
# 
#
p.transfer(water_tube,       re_tube, ul(117))
p.provision(inv["CutSmart"], re_tube, ul(15))
p.provision(inv["pUC19"],    re_tube, ul(3))
p.mix(re_tube, volume=ul(60), repetitions=10)
assert re_tube.volume == ul(120) + dead_volume["micro-1.5"]

print("Volumes: re_tube:{} water_tube:{} EcoRI:{}".format(re_tube.volume, water_tube.volume, ecori_p10x_tube.volume))

p.distribute(re_tube,         pcr_plate.wells(["A1","B1","A2"]), ul(40))
p.distribute(water_tube,      pcr_plate.wells(["A2"]),           ul(10))
p.distribute(ecori_tube, pcr_plate.wells(["A1","B1"]),      ul(10))
assert all(well.volume == ul(50) for well in pcr_plate.wells(["A1","B1","A2"]))

p.mix(pcr_plate.wells(["A1","B1","A2"]), volume=ul(25), repetitions=10)

# Incubation to induce cut, then heat inactivation of EcoRI
p.seal(pcr_plate)
#NEB enzymes are time saver so they only need 15minutes
p.incubate(pcr_plate, "warm_37", "15:minute", shaking=False)
p.thermocycle(pcr_plate, [{"cycles":  1, "steps": [{"temperature": "65:celsius", "duration": "21:minute"}]}], volume=ul(50))

# --------------------------------------------------------------
# Gel electrophoresis, to ensure the cutting worked
#
p.unseal(pcr_plate)
p.mix(pcr_plate.wells(["A1","B1","A2"]), volume=ul(25), repetitions=5)
p.transfer(pcr_plate.wells(["A1","B1","A2"]), pcr_plate.wells(["D1","E1","D2"]), ul(8))
p.transfer(water_tube, pcr_plate.wells(["D1","E1","D2"]), ul(15), mix_after=True, mix_vol=ul(10))
assert all(well.volume == ul(20) + dead_volume["96-pcr"] for well in pcr_plate.wells(["D1","E1","D2"]))

p.gel_separate(pcr_plate.wells(["D1","E1","D2"]), ul(20), "agarose(10,2%)", "ladder2", "15:minute", expid("gel"))


# ----------------------------------------------------------------------------
# Then consolidate all cut plasmid to one tube (puc19_cut_tube).
#
remaining_volumes = [well.volume - dead_volume['96-pcr'] for well in pcr_plate.wells(["A1","B1"])]
print("Consolidated volume: {}".format(sum(remaining_volumes, ul(0))))
p.consolidate(pcr_plate.wells(["A1","B1"]), puc19_cut_tube, remaining_volumes, allow_carryover=True)

assert all(tube.volume >= dead_volume['micro-1.5'] for tube in [water_tube, re_tube, puc19_cut_tube, ecori_p10x_tube])

# ---------------------------------------------------------------
# Test protocol
#
print(json.dumps(p.as_dict(), indent=2))