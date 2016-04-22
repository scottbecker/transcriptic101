import json
from autoprotocol.protocol import Protocol
from utils import ul

inv = {
    'te': 'rs17pwyc754v9t',           # catalog; TE
}

def get_synthesized_oligo_tube_protocol(tube_name, sequence):
    """

    Synthesize and prepare an oligo tube
    25nm seems to be more than enough for most downstream work (25pmol seems to be all that is used for transcriptic's pcr)
    
    """
    global inv
    
    scale = '25nm'
    
    p = Protocol()

    dna_tube = p.ref(tube_name, 
                       cont_type="micro-1.5",
                       storage="cold_20", discard=False)
    
    p.oligosynthesize([{"sequence": sequence,
                        "destination": dna_tube.well(0),
                        "scale": scale,
                        "purification": "standard"}]
                      )
    
    #spin
    p.spin(dna_tube, '2000:g', '30:second')
    
    
    #dilute to 100uM
    #safe min volume of 1.5uL is 20uL so we want a lot more than this so we don't lose too much
    #how do you go from scale and desired concentration to volume --> n = CV --> V = n/C --> V = 25nmol / 100uM = 25nmol / (100E3 nmol / 1L)
    #          = 2.5E-4 L = 250 uL
    
                    #convert to nM     #convert to uL
    te_volume = 25 / pow(100, 3) * pow(10, 6) 
    
    #add 250uL
    p.provision(inv["te"], dna_tube.well(0), ul(250))

    #spin
    p.spin(dna_tube, '2000:g', '30:second')
    
    return json.dumps(p.as_dict(), indent=2)