import json
from autoprotocol.protocol import Protocol

from utils import ul

"""
This protoocol will take a tube and fill it with water

"""


inv = {
     'te':'rs17pwyc754v9t'    #catalog: te
}

p = Protocol()

#create a container 

te_well = p.ref("te_tube", cont_type="micro-1.5", storage="cold_4", discard=False).well(0)

p.provision(inv["te"], te_well, ul(1500))


print(json.dumps(p.as_dict(), indent=2))