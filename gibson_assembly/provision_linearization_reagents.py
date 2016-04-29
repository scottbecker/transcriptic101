import json
from autoprotocol.protocol import Protocol
import sys
from utils import ul

"""
This protoocol will provision restriction enzyme tubes

We provision in advance to reduce thaw/unthaw cycles of the main reagents

NEB recommends storing at selling concentrations (so we don't dilute)

A1 - Mix of EcoRI and hindiii at 10 units of each re / ul
B1 - CutSmart
H12 - pUC19



To retrieve from this reagent plate:

reagent_plate = p.ref("re_reagent_plate", id=inv['reagent_plate'], 
                      cont_type="96-pcr", storage="cold_20")

ecori_hindiii_well = reagent_plate.wells(["A1"])[0]
cutsmart_well = reagent_plate.wells(["B1"])[0]
pUC19_well = reagent_plate.wells(["H12"])[0]


"""


_options = {
    'use_existing_reagent_plate': False
}
options = {k for k,v in _options.items() if v is True}

plate_name = 're_reagents'
number_of_experiments = 6

inv = {
    "EcoRI":    "rs17ta8xftpdk6",   # catalog; EcoRI-HF; cold_20; 20 units/ul
    "hindiii":    "rs18nw6kpnp44v",   # catalog; hindiii-HF; cold_20; 20 units/ul
    "CutSmart": "rs17ta93g3y85t",   # catalog; CutSmart Buffer 10x; cold_20
    "pUC19":    "rs17tcqmncjfsh",   # catalog; pUC19; cold_20
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

ecori_hindiii_well = reagent_plate.wells(["A1"])[0]
#provisioning less than 10uL is dangerous due to pipetting mistakes (can't mix after)
volume_per_experiment = ul(1)
re_volume = max(ul(10),number_of_experiments*volume_per_experiment+dead_volume/2)
p.provision(inv["EcoRI"], ecori_hindiii_well, re_volume)
p.provision(inv["hindiii"], ecori_hindiii_well, re_volume)
ecori_hindiii_well.name = 'EcoRI+hindiii'
ecori_hindiii_well.properties = {'Unit Concentration':'10 units/ul'} 

#there is an extra negative control for every 2 experiments so we use 3/2*experiments
cutsmart_well = reagent_plate.wells(["B1"])[0]
volume_per_experiment = ul(5)
p.provision(inv["CutSmart"], cutsmart_well, number_of_experiments*3/2*volume_per_experiment+dead_volume)
cutsmart_well.name = 'CutSmart 10x'

pUC19_well = reagent_plate.wells(["H12"])[0]
volume_per_experiment = ul(1)
pUC19_volume = max(ul(10),number_of_experiments*volume_per_experiment+dead_volume)
p.provision(inv["pUC19"], pUC19_well, pUC19_volume)
pUC19_well.name = 'pUC19'
pUC19_well.properties = {'Mass Concentration':'1 ug/ul'} 

print(json.dumps(p.as_dict(), indent=2))