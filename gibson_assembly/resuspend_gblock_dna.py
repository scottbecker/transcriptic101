import json
from autoprotocol.protocol import Protocol


inv = {
#    'gblock_dna':'ct18vs7cmjat2c'   #do not use until we want to do this for real
     'gblock_dna':''    #

}

p = Protocol()



#ct18vs7cmjat2c

dna_sample = p.ref("bacteria", None, "micro-1.5", discard=True)
test_plate = p.ref("test_plate", None, "96-flat", storage="cold_4")

p.dispense_full_plate(test_plate, "lb-broth-noAB", "50:microliter")
w = 0
amt = 1

while amt < 20:
    p.transfer(bacterial_sample.well(0), test_plate.well(w), "%d:microliter" % amt)
    amt += 2
    w +=1

print(json.dumps(p.as_dict(), indent=2))