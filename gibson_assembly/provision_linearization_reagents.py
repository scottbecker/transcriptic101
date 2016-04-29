import json
from autoprotocol.protocol import Protocol
import sys
from utils import ul

"""
This protoocol will provision restriction enzyme tubes

We provision in advance to reduce thaw/unthaw cycles of the main reagents

NEB recommends storing at selling concentrations (so we don't dilute)

"""


_options = {
    'use_existing_reagent_plate': False
}
options = {k for k,v in _options.items() if v is True}

plate_name = 're_reagents'
number_of_experiments = 6

inv = {
    "EcoRI":    "rs17ta8xftpdk6",   # catalog; EcoRI-HF; cold_20; 20 units/ul
    "BamHI":    "rs17ta8tz5ffby",   # catalog; BamHI-HF; cold_20; 20 units/ul
    "CutSmart": "rs17ta93g3y85t",   # catalog; CutSmart Buffer 10x; cold_20
    'reagent_plate': '' #optional: only needed after the first run
}

if "--test" in sys.argv:
    test_inv = {
        'reagent_plate': '' #optional: only needed after the first run
    }
    inv.update(test_inv)

p = Protocol()

if 'use_existing_reagent_plate' not in options:
    dead_volume = ul(3)
    reagent_plate = p.ref(plate_name, cont_type="96-pcr", storage="cold_20")
else:
    dead_volume = ul(0)
    reagent_plate = p.ref(plate_name, id=inv['reagent_plate'], 
                              cont_type="96-pcr", storage="cold_20")

ecori_bamhi_well = reagent_plate.wells(["A1"])[0]
#provisioning less than 10uL is dangerous do to pipetting mistakes (can't mix after)
re_volume = max(ul(10),ul(number_of_experiments*1)+dead_volume/2)
p.provision(inv["EcoRI"], ecori_bamhi_well, re_volume)
p.provision(inv["BamHI"], ecori_bamhi_well, re_volume)
ecori_bamhi_well.name = 'EcoRI+BamHI'
ecori_bamhi_well.properties = {'Unit Concentration':'10 units/ul'} 

#there is an extra negative control for every 2 experiments so we use 3/2*experiments
cutsmart_well = reagent_plate.wells(["B1"])[0]
p.provision(inv["CutSmart"], cutsmart_well, ul((number_of_experiments*3/2)*5)+dead_volume)
cutsmart_well.name = 'CutSmart 10x'

print(json.dumps(p.as_dict(), indent=2))