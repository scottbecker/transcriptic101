import json
import sys
from autoprotocol.protocol import Protocol

from utils import ul

"""
This protoocol will take create tubes with 1.5mL of reagents for PCR

"""

_options = {
    'include_Q5_set': True,
    'include_dNTP_Mixture_10mM': True,
    'include_exosap': False,
    'use_existing_pcr_reagent_plate': False
}
options = {k for k,v in _options.items() if v is True}

inv = {
    'Q5_Polymerase':                     'rs16pcce8rdytv', # catalog; Q5 High-Fidelity DNA Polymerase
    'Q5_Buffer':                         'rs16pcce8rmke3', # catalog; Q5 Reaction Buffer
    'Q5_Enhancer':                       'rs16pcce8rva4a', # catalog; Q5 Reaction Buffer    
    'dNTP_Mixture_10mM':                 'rs186wj7fvknsr', # catalog; dNTP Mixture 10mM
    'ExoSAP_IT':                         'rs18dnrskds4t6',  # catalog; ExoSAP-IT
    'pcr_reagent_plate':                  'ct18x95j499zx4' #optional: only needed after the first run
}

if "--test" in sys.argv:
    test_inv = {
        'pcr_reagent_plate':              'ct18x9bmmam5xw' #optional: only needed after the first run
    }
    inv.update(test_inv)

p = Protocol()

if 'use_existing_pcr_reagent_plate' not in options:
    #we use a PCR plate since it has less dead volume (we would be throwing out $30 if we use normal tubes)
    #pcr well has 3uL of dead volume and 5uL safe min (you can't pipette from less than this)
    pcr_reagent_plate = p.ref("pcr_reagent_plate", cont_type="96-pcr", storage="cold_20",  discard=False)
else:
    pcr_reagent_plate = p.ref("pcr_reagent_plate", id=inv['pcr_reagent_plate'], 
                              cont_type="96-pcr", storage="cold_20",  discard=False)

if 'include_Q5_set' in options:
    #We use 1uL per set of 4 experiments. To get 4 sets of experiments, we need 8uL
    q5_poly_well = pcr_reagent_plate.wells(["A1"])[0]
    p.provision(inv["Q5_Polymerase"], q5_poly_well, ul(14))
    q5_poly_well.name = 'Q5_Polymerase'

if 'include_Q5_set' in options:
    #We use 20uL per 4 experiments. 3 of dead volume + 20 uL * 4 sets of 4 experiments =  
    q5_buffer_well = pcr_reagent_plate.wells(["A2"])[0]
    p.provision(inv["Q5_Buffer"], q5_buffer_well, ul(83))
    q5_buffer_well.name = 'Q5_Buffer'

if 'include_Q5_set' in options:
    #we aren't using this but you get it for 20 cents because it comes with the Q5 kit
    q5_enhancer_well = pcr_reagent_plate.wells(["A3"])[0]
    p.provision(inv["Q5_Enhancer"], q5_enhancer_well, ul(83))
    q5_enhancer_well.name = 'Q5_Enhancer'

if 'include_dNTP_Mixture_10mM' in options:
    #We use 2uL per 4 experiments. 3 of dead volume + 2 uL * 4 sets of 4 experiments =  
    dNTP_well = pcr_reagent_plate.wells(["B1"])[0]
    p.provision(inv["dNTP_Mixture_10mM"], dNTP_well, ul(11))
    dNTP_well.name = 'dNTP_Mixture_10mM'
    dNTP_well.properties = {'Molar Concentration':'10mM'}    

if 'include_exosap' in options:
    exosap_well = pcr_reagent_plate.wells(["C1"])[0]
    p.provision(inv["ExoSAP_IT"], exosap_well, ul(11))
    exosap_well.name = 'ExoSAP_IT'    


print(json.dumps(p.as_dict(), indent=2))