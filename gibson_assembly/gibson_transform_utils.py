import sys
from utils import (ul, expid, init_inventory_well, touchdown,
                   dead_volume )
def get_inventory():
    inv = {
        # my inventory
        "puc19_ecori_hindiii_puc19_cut": "ct18xzbzsz45b6", # inventory; pUC19 cut with EcoRI; cold_20
        "sfgfp_pcr_ecori_hindiii_amplified"    : "ct18xzndsv4yqs", # inventory; sfGFP amplified to 40ng/ul; cold_20
        "sfgfp_puc19_gibson_clone" : "ct18y2qnrjqz7w", # inventory; assembled sfGFP; cold_20
    }
    
    if "--test" in sys.argv:
        test_inv =  {
            # my inventory
            "puc19_ecori_hindiii_puc19_cut": "ct18xkhwbwrbj6", # inventory; pUC19 cut with EcoRI; cold_20
            "sfgfp_pcr_ecori_hindiii_amplified"    : "ct18xkegtgxzjw", # inventory; sfGFP amplified to 40ng/ul; cold_20
            "sfgfp_puc19_gibson_clone" : "ct18y29ytg9rd9", # inventory; assembled sfGFP; cold_20
        }
        inv.update(test_inv)
    
    return inv


# -----------------------------------------------------------------------------
# Cloning/assembly (see NEBuilder protocol above)
#
# "Optimized efficiency is 50–100 ng of vectors with 2 fold excess of inserts."
# pUC19 is 20ng/ul (78ul total).
# sfGFP is ~25ng/ul (48ul total)
# Therefore 2.5ul of pUC19 and 4ul sfGFP of 50ng and 100ng of vector and insert respectively
#

def do_gibson_assembly(p, water_tube, clone_plate, puc19_cut_tube, sfgfp_pcroe_amp_tube):
    #
    # Combine all the Gibson reagents in one tube and thermocycle
    #
    p.provision(p.inv["NEBuilder_Master_Mix"], clone_plate.well(0), ul(10))
    p.transfer(water_tube,           clone_plate.well(0), ul(3.5),
               mix_after=True, mix_vol=ul(6))
    p.transfer(puc19_cut_tube,       clone_plate.well(0), ul(2.5),
               mix_after=True, mix_vol=ul(6))
    p.transfer(sfgfp_pcroe_amp_tube, clone_plate.well(0), ul(4),
               mix_after=True, mix_vol=ul(10))

    p.seal(clone_plate)
    p.thermocycle(clone_plate,
                  [{"cycles":  1, "steps": [{"temperature": "50:celsius", "duration": "15:minute"}]}],
                  volume=ul(20))

    #
    # Dilute assembled plasmid 4X according to the NEB Gibson assembly protocol (20ul->80ul)
    #
    p.unseal(clone_plate)
    p.transfer(water_tube, clone_plate.well(0), ul(60), mix_after=True, mix_vol=ul(40), repetitions=5)


# --------------------------------------------------------------------------------------------------
# Transformation
# "Transform NEB 5-alpha Competent E. coli cells with 2 ul of the
#  assembled product, following the appropriate transformation protocol."
#
# Mix & Go http://www.zymoresearch.com/downloads/dl/file/id/173/t3015i.pdf
# "[After mixing] Immediately place on ice and incubate for 2-5 minutes"
# "The highest transformation efficiencies can be obtained by incubating Mix & Go cells with DNA on
#  ice for 2-5 minutes (60 minutes maximum) prior to plating."
# "It is recommended that culture plates be pre-warmed to >20°C (preferably 37°C) prior to plating."
# "Avoid exposing the cells to room temperature for more than a few seconds at a time."
#
# "If competent cells are purchased from other manufacture, dilute assembled products 4-fold
#  with H2O prior transformation. This can be achieved by mixing 5 ul of assembled products with
#  15 ul of H2O. Add 2 ul of the diluted assembled product to competent cells."
#

def _do_transformation(p, clone_plate,transform_plate, transform_tube,
                       transctrl_tube, control_pUC19):
    #
    # Combine plasmid and competent bacteria in a pcr_plate and shock
    #
    p.provision(p.inv["DH5a"], transform_tube,  ul(50))
    assert clone_plate.well(0).volume >= ul(6), clone_plate.well(0).volume
    #This was changed from the original publication of 3uL due to low colony turnout
    p.transfer(clone_plate.well(0), transform_tube, ul(5), dispense_speed="10:microliter/second",
               mix_after=True,mix_vol=ul(10))
    

    if control_pUC19:
        p.provision(p.inv["DH5a"], transctrl_tube,  ul(50))
        p.provision(p.inv["pUC19"], transctrl_tube, ul(1))

    #
    # Heatshock the bacteria to transform using a PCR machine
    #
    p.seal(transform_plate)
    p.thermocycle(transform_plate,
        [{"cycles":  1, "steps": [{"temperature":  "4:celsius", "duration": "5:minute"}]},
         {"cycles":  1, "steps": [{"temperature": "37:celsius", "duration": "30:minute"}]}],
        volume=ul(50))


def _transfer_transformed_to_plates(p, transform_plate, transform_tube,
                                    transform_tube_L, transctrl_tube,
                                    transctrl_tube_L, amp_6_flat,
                                    noAB_6_flat,control_pUC19):
    assert transform_tube.volume == ul(53), transform_tube.volume
    p.unseal(transform_plate)

    num_ab_plates = 4 # antibiotic places

    #
    # Transfer bacteria to a bigger tube for diluting
    # Then spread onto 6-flat plates
    # Generally you would spread 50-100ul of diluted bacteria
    # Put more on ampicillin plates for more opportunities to get a colony
    # I use a dilution series since it's unclear how much to plate
    #
    p.provision(p.inv["LB Miller"], transform_tube_L, ul(445))

    #
    # Add the transformed cells and mix (use new mix op in case of different pipette)
    #
    p.transfer(transform_tube, transform_tube_L, ul(50), mix_before=True)
    p.mix(transform_tube_L, volume=transform_tube_L.volume/2, repetitions=10)

    assert transform_tube.volume == dead_volume['96-pcr'] == ul(3), transform_tube.volume
    assert transform_tube_L.volume == ul(495), transform_tube_L.volume

    #
    # Spread an average of 60ul on each plate == 480ul total
    #
    for i in range(num_ab_plates):
        p.spread(transform_tube_L, amp_6_flat.well(i), ul(51+i*6))
        p.spread(transform_tube_L, noAB_6_flat.well(i), ul(51+i*6))

    assert transform_tube_L.volume == dead_volume["micro-1.5"], transform_tube_L.volume

    #
    # Controls: include 2 ordinary pUC19-transformed plates as a control
    #
    if control_pUC19:
        num_ctrl = 2
        assert num_ab_plates + num_ctrl <= 6

        p.provision(p.inv["LB Miller"], transctrl_tube_L, ul(192)+dead_volume["micro-1.5"])
        p.transfer(transctrl_tube,    transctrl_tube_L, ul(48))
        p.mix(transctrl_tube_L, volume=transctrl_tube_L.volume/2, repetitions=10)

        for i in range(num_ctrl):
            p.spread(transctrl_tube_L, amp_6_flat.well(num_ab_plates+i), ul(55+i*10))
            p.spread(transctrl_tube_L, noAB_6_flat.well(num_ab_plates+i), ul(55+i*10))

        assert transctrl_tube_L.volume >= dead_volume["micro-1.5"], transctrl_tube_L.volume

def do_transformation(p,clone_plate, transform_plate, transform_tube,
                      transform_tube_L, transctrl_tube,
                      transctrl_tube_L, amp_6_flat,
                      noAB_6_flat, control_pUC19=True):
    _do_transformation(p, clone_plate,transform_plate, transform_tube,
                       transctrl_tube, control_pUC19)
    _transfer_transformed_to_plates(p, transform_plate, transform_tube,
                                    transform_tube_L, transctrl_tube,
                                    transctrl_tube_L, amp_6_flat,
                                    noAB_6_flat,control_pUC19)


# ------------------------------------------------------
# Measure growth in plates (photograph)
#

def measure_growth(p, amp_6_flat, noAB_6_flat):
    #
    # Incubate and photograph 6-flat plates over 18 hours
    # to see blue or white colonies
    #
    for flat_name, flat in [(expid("amp_6_flat"), amp_6_flat), (expid("noAB_6_flat"), noAB_6_flat)]:
        for timepoint in [9,18]:
            p.cover(flat)
            p.incubate(flat, "warm_37", "9:hour")
            p.uncover(flat)
            p.image_plate(flat, mode="top", dataref=expid("{}_t{}".format(flat_name, timepoint)))


