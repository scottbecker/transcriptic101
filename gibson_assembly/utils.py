import autoprotocol
from autoprotocol import Unit
from autoprotocol.container import Container
from autoprotocol.protocol import Protocol
from autoprotocol.protocol import Ref # "Link a ref name (string) to a Container instance."
import requests
import logging
import json

# Transcriptic authorization
org_name = 'hgbrian'
tsc_headers = {k:v for k,v in json.load(open("../auth.json")).items() if k in ["X_User_Email","X_User_Token"]}

# Transcriptic-specific dead volumes
_dead_volume = [("96-pcr",3), ("96-flat",25), ("96-flat-uv",25), ("96-deep",15),
                ("384-pcr",2), ("384-flat",5), ("384-echo",15),
                ("micro-1.5",15), ("micro-2.0",15)]
dead_volume = {k:Unit(v,"microliter") for k,v in _dead_volume}


def init_inventory_well(well, headers=tsc_headers, org_name=org_name):
    """Initialize well (set volume etc) for Transcriptic"""
    def _container_url(container_id):
        return 'https://secure.transcriptic.com/{}/samples/{}.json'.format(org_name, container_id)

    response = requests.get(_container_url(well.container.id), headers=headers)
    response.raise_for_status()

    container = response.json()
    well_data = container['aliquots'][well.index]
    well.name = "{}/{}".format(container["label"], well_data['name']) if well_data['name'] is not None else container["label"]
    well.properties = well_data['properties']
    well.volume = Unit(well_data['volume_ul'], 'microliter')

    if 'ERROR' in well.properties:
        raise ValueError("Well {} has ERROR property: {}".format(well, well.properties["ERROR"]))
    if well.volume < Unit(20, "microliter"):
        logging.warn("Low volume for well {} : {}".format(well.name, well.volume))

    return True

def touchdown(fromC, toC, durations, stepsize=2, meltC=98, extC=72):
    """Touchdown PCR protocol generator"""
    assert 0 < stepsize < toC < fromC
    def td(temp, dur): return {"temperature":"{:2g}:celsius".format(temp), "duration":"{:d}:second".format(dur)}

    return [{"cycles": 1, "steps": [td(meltC, durations[0]), td(C, durations[1]), td(extC, durations[2])]}
            for C in np.arange(fromC, toC-stepsize, -stepsize)]

def convert_ug_to_pmol(ug_dsDNA, num_nts):
    """Convert ug dsDNA to pmol"""
    return float(ug_dsDNA)/num_nts * (1e6 / 660.0)

def expid(val):
    """Generate a unique ID per experiment"""
    return "{}_{}".format(experiment_name, val)

def ul(microliters):
    """Unicode function name for creating microliter volumes"""
    return Unit(microliters,"microliter")