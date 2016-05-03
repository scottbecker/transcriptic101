import sys
import json
from custom_protocol import CustomProtocol as Protocol
from utils import (ul, expid, init_inventory_well, touchdown,
                   dead_volume)
from gibson_transform_utils import do_transformation, get_inventory, measure_growth, do_sanger_seq

"""Transform after gibson assembly for sfGFP and pUC19"""

experiment_name = "sfgfp_puc19_transform_v1"

p = Protocol()

inv = get_inventory()

clone_plate = p.ref("sfgfp_puc19_gibson_v3_clone", id=inv["sfgfp_puc19_gibson_v3_clone"],
                    cont_type="96-pcr", storage="cold_4")
#
# Catalog (all to be discarded afterward)
#
water_tube       = p.ref("water",     cont_type="micro-1.5", storage="ambient", discard=True).well(0)
transform_plate  = p.ref("trn_plate", cont_type="96-pcr",    storage="ambient", discard=True)
transform_tube   = transform_plate.well(39) # experiment
transform_tube_L = p.ref("trn_tubeL", cont_type="micro-1.5", storage="ambient", discard=True).well(0)
transctrl_tube   = transform_plate.well(56) # control
transctrl_tube_L = p.ref("trc_tubeL", cont_type="micro-1.5", storage="ambient", discard=True).well(0)

#
# Plating according to Tali's protocol
# http://learn.transcriptic.com/blog/2015/9/9/provisioning-commercial-reagents
#
amp_6_flat = Container(None, p.container_type('6-flat'))
p.refs[expid("amp_6_flat")] = Ref(expid("amp_6_flat"),
                                  {"reserve": inv['lb-broth-100ug-ml-amp_6-flat'], "store": {"where": 'cold_4'}}, amp_6_flat)
noAB_6_flat = Container(None, p.container_type('6-flat'))
p.refs[expid("noAB_6_flat")] = Ref(expid("noAB_6_flat"),
                                   {"reserve": inv['noAB-amp_6-flat'], "store": {"where": 'cold_4'}}, noAB_6_flat)

#
# Initialize inventory
#

all_inventory_wells = [IPTG_tube, clone_plate.well(0)]

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

do_transformation(p)
measure_growth(p)

# ---------------------------------------------------------------
# Output protocol
#
jprotocol = json.dumps(p.as_dict(), indent=2)
print(jprotocol)