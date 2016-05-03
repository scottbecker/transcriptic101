from autoprotocol.protocol import Protocol
import json
import sys

from utils import ul, init_inventory_well, dead_volume

inv = {
    'water':                                  'rs17gmh5wafm5p', # catalog: autoclave milliq h2o
    'te':                                     'rs17pwyc754v9t', # catalog; TE
    'sfgfp_puc19_primer_forward_100uM':       'ct18wxfe2dxmmd', # inventory; micro-1.5, cold_4, 100uM
    'sfgfp_puc19_primer_reverse_hindiii_100uM':  'ct18xm348hdmrb', # inventory; micro-1.5, cold_4, 100uM
}

if "--test" in sys.argv:
    test_inv = {
        'sfgfp_puc19_primer_forward_100uM':  'ct18x624ssz78h', # inventory; micro-1.5, cold_4, 100uM
        'sfgfp_puc19_primer_reverse_hindiii_100uM':  'ct18xke9k7wmkj', # inventory; micro-1.5, cold_4, 100uM
    }
    inv.update(test_inv)



def dilute_primer(source_well, destination_well, destination_volume_uL,
                  ratio, diluent):
    """
    Dilutes source_well into destination tube using diluent
    
    Args:
        source_well: the tube at a higher concentration
        destination_well: the tube that will have the dilution product
        destination_volume_uL: the volume of destination_well
        ratio:  destination_well concentration / source_well concentration
        diluent_well: a tube with enough d to be used for the diluation
    
    Returns:
        None
    
    """
    global inv
    
    #@TODO: figure out how to check if a diluent tube has enough volume in it to conduct the action 
    # (needs to account for multiple calls to this function)
    
    mix_volume = destination_volume_uL / 2
    source_volume_uL = ratio * destination_volume_uL
    diluent_volume_uL = destination_volume_uL - source_volume_uL

    p.provision(inv[diluent], destination_well, ul(diluent_volume_uL))

    p.transfer(source_well, destination_well, ul(source_volume_uL), 
               mix_before=True, 
               mix_after=True,
               mix_vol=min(ul(source_volume_uL*5),
                           #prevent drawing more from the source well than is retrievable
                           source_well.volume-source_well.container.container_type.dead_volume_ul))
    p.mix(destination_well, volume=ul(destination_volume_uL/2.0), repetitions=10)
    

#@TODO: how do you create a tube with properties

#make destination tubes
p = Protocol()
#remove ## after oligosynth is done
dilute_primer_wells = [p.ref('sfgfp_puc19_primer_forward_10uM', cont_type="micro-1.5", storage="cold_20").well(0),
                       p.ref('sfgfp_puc19_primer_reverse_hindiii_10uM', cont_type="micro-1.5", storage="cold_20").well(0)]

#set concentration property on wells
for well in dilute_primer_wells:
    well.properties = {'Molar Concentration':'10uM'}    

#remove ## after oligosynth is done
primer_wells = [p.ref('sfgfp_puc19_primer_forward_100uM', id=inv['sfgfp_puc19_primer_forward_100uM'], 
                      cont_type="micro-1.5", storage="cold_20").well(0),
                p.ref('sfgfp_puc19_primer_reverse_hindiii_100uM', id=inv['sfgfp_puc19_primer_reverse_hindiii_100uM'], 
                             cont_type="micro-1.5", storage="cold_20").well(0)]

for well in primer_wells:
    init_inventory_well(well)

for x in range(0,len(primer_wells)):
    dilute_primer(primer_wells[x], dilute_primer_wells[x], 
                  100, 0.1, 'te')


print(json.dumps(p.as_dict(), indent=2))