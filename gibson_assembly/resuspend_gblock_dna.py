import json
from autoprotocol.protocol import Protocol
from utils import ul

inv = {
    'gblock_dna':'ct18vs7cmjat2c',     # my inventory
     'te': 'rs17pwyc754v9t',           # catalog; TE
}

p = Protocol()

#existing inventory
dna_sample = p.ref("gblock_dna", id=inv['gblock_dna'], cont_type="micro-1.5", storage='cold_20')

#tubes are thawed before being placed into workcell

#Centrifuge the tube for 3-5 sec at a minimum of 3000 x g to pellet the material to the bottom of the tube.
p.spin(dna_sample, '3000:g', '5:second')

#Add 500uL TE (to reach 1 ng/uL)
p.provision(inv["te"], dna_sample.well(0), ul(500))

#Briefly vortex (here we have to use mix) and centrifuge

p.mix(dna_sample.well(0),'480:microliter',speed='480:microliter/second',repetitions=10)

p.spin(dna_sample, '3000:g', '3:second')

print(json.dumps(p.as_dict(), indent=2))