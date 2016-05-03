import sys
import json
from custom_protocol import CustomProtocol as Protocol
from utils import (ul, expid, init_inventory_well, touchdown,
                   dead_volume)
from gibson_transform_utils import do_gibson_assembly, get_inventory

"""Full Gibson assembly and transformation protocol for sfGFP and pUC19"""


experiment_name = "sfgfp_puc19_gibson_v1"

p = Protocol()

inv = get_inventory()

puc19_cut_tube = p.ref("puc19_ecori_hindiii_puc19_cut", id=inv["puc19_ecori_hindiii_puc19_cut"],
                       cont_type="micro-1.5", storage="cold_20").well(0)
sfgfp_pcroe_amp_tube = p.ref("sfgfp_pcroe_v8_amplified", id=inv["sfgfp_pcroe_v8_amplified"],
                             cont_type="micro-1.5", storage="cold_4").well(0)
clone_plate = p.ref(expid("clone"), cont_type="96-pcr", storage="cold_4")

#
# Catalog (all to be discarded afterward)
#

water_tube       = p.ref("water",     cont_type="micro-1.5", discard=True).well(0)

#
# Initialize inventory
#

all_inventory_wells = [puc19_cut_tube, sfgfp_pcroe_amp_tube, IPTG_tube]
assert puc19_cut_tube.volume == ul(66), puc19_cut_tube.volume
assert sfgfp_pcroe_amp_tube.volume == ul(36), sfgfp_pcroe_amp_tube.volume

for well in all_inventory_wells:
    init_inventory_well(well)
    print("Inventory: {} {} {}".format(well.name, well.volume, well.properties))

    
#
# Provisioning. Water is used all over the protocol. Provision an excess since it's cheap
#
p.provision(inv["water"], water_tube, ul(500))

# ---------------------------------------------------------------
# Generate protocol
#

do_gibson_assembly(p, water_tube, clone_plate, puc19_cut_tube, sfgfp_pcroe_amp_tube)

do_sanger_seq(p)

# ---------------------------------------------------------------
# Output protocol
#
jprotocol = json.dumps(p.as_dict(), indent=2)
print(jprotocol)