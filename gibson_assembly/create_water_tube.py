import json
from autoprotocol.protocol import Protocol

from utils import ul

"""
This protoocol will take a tube and fill it with water

"""


inv = {
     'water':'rs17gmh5wafm5p'    #autoclave milliq h2o
}

p = Protocol()

#create a container 

water_tube = p.ref("water_tube", cont_type="micro-1.5", storage="cold_4", discard=False).well(0)

p.provision(inv["water"], water_tube, ul(1500))


print(json.dumps(p.as_dict(), indent=2))