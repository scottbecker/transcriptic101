"""Pick colonies from plates and grow in amp media and check for fluorescence.
v2: try again with a new plate (no blue colonies)
v3: repeat with different emission and excitation wavelengths"""

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

experiment_name = "sfgfp_puc19_gibson_pick_v1"
utils.experiment_name = experiment_name 
def plate_expid(val):
    """refer to the previous plating experiment's outputs"""
    plate_exp = "sfgfp_puc19_transform_v1"
    return "{}_{}".format(plate_exp, val)

# -----------------------------------------------------------------------
# Inventory
#
inv = {
    # plates from previous experiment, must be changed every new experiment
    plate_expid("amp_6_flat")  : "ct18yn63tz9wv9", # inventory; Ampicillin plates with transformed bacteria
}

if "--test" in sys.argv:
    test_inv =  {
        plate_expid("amp_6_flat")  : "ct18y4w8dz47pa", # inventory; Ampicillin plates with transformed bacteria
    }
    inv.update(test_inv)

# Tubes and plates
lb_amp_tubes = [p.ref("lb_amp_{}".format(i+1), cont_type="micro-2.0", storage="ambient", discard=True).well(0)
                for i in range(4)]
lb_xab_tube  = p.ref("lb_xab", cont_type="micro-2.0", storage="ambient", discard=True).well(0)
growth_plate = p.ref(expid("growth"), cont_type="96-flat",   storage="cold_4",  discard=False)

# ampicillin plate
amp_6_flat = Container(None, p.container_type('6-flat'))
p.refs[plate_expid("amp_6_flat")] = Ref(plate_expid("amp_6_flat"),
                                        {"id":inv[plate_expid("amp_6_flat")], 
                                         "store": {"where": 'cold_4'}},
                                        amp_6_flat)

# Use a total of 50 wells
abs_wells = ["{}{}".format(row,col) for row in "BCDEF" for col in range(1,11)]
abs_wells_T = ["{}{}".format(row,col) for col in range(1,11) for row in "BCDEF"]
assert abs_wells[:3] == ["B1","B2","B3"] and abs_wells_T[:3] == ["B1","C1","D1"]

def prepare_growth_wells():
    #
    # To LB, add ampicillin at ~1/1000 concentration
    # Mix slowly in case of overflow
    #
    p.provision(p.inv["LB Miller"], lb_xab_tube, ul(1913))
    for lb_amp_tube in lb_amp_tubes:
        p.provision(p.inv["LB Miller"],   lb_amp_tube, ul(1911))
        p.provision(p.inv["Amp 100mgml"], lb_amp_tube, ul(2))
        p.mix(lb_amp_tube, volume=ul(800), repetitions=10)

    #
    # Add IPTG but save on X-Gal
    # http://openwetware.org/images/f/f1/Dh5a_sub.pdf
    # "If you are concerned about obtaining maximal levels of expression, add IPTG to a final concentration of 1 mM."
    # 20ul of IPTG @ 100uM in ~2000ul equals 1mM
    #
    p.provision(p.inv['IPTG'], [lb_xab_tube] + lb_amp_tubes, ul(20))

    #
    # Distribute LB among wells, row D is control (no ampicillin)
    #
    cols = range(1,11)
    row = "D" # control, no AB
    cwells = ["{}{}".format(row,col) for col in cols]
    assert set(cwells).issubset(set(abs_wells))
    p.distribute(lb_xab_tube,  growth_plate.wells(cwells), ul(191.8), allow_carryover=True)

    rows = "BCEF"
    for row, lb_amp_tube in zip(rows, lb_amp_tubes):
        cwells = ["{}{}".format(row,col) for col in cols]
        assert set(cwells).issubset(set(abs_wells))
        p.distribute(lb_amp_tube, growth_plate.wells(cwells), ul(191.8), allow_carryover=True)

    assert all(round(lb_amp_tube.volume,0) == round(lb_xab_tube.volume,0) == dead_volume['micro-2.0']
               for lb_amp_tube in lb_amp_tubes), [lb_amp_tube.volume for lb_amp_tube in lb_amp_tubes]
    return


def measure_growth_wells():
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

        p.fluorescence(growth_plate, growth_plate.wells(abs_wells).indices(),
                       excitation="485:nanometer", emission="535:nanometer",
                       dataref=expid("fl2_{}".format(t)), flashes=25)
        p.absorbance(growth_plate, growth_plate.wells(abs_wells).indices(),
                     wavelength="600:nanometer",
                     dataref=expid("abs_{}".format(t)), flashes=25)
    

# ---------------------------------------------------------------
# Protocol steps
#
prepare_growth_wells()
batch = 10
for i in range(5):
    p.autopick(amp_6_flat.well(i), growth_plate.wells(abs_wells_T[i*batch:i*batch+batch]),
               dataref=expid("autopick_{}".format(i)))
    p.image_plate(amp_6_flat, mode="top", dataref=expid("autopicked_{}".format(i)))
measure_growth_wells()

# ---------------------------------------------------------------
# Output protocol
#
jprotocol = json.dumps(p.as_dict(), indent=2)
print(jprotocol)