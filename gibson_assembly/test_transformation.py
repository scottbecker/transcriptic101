"""Simple transformation protocol: transformation with unaltered pUC19"""

p = Protocol()

experiment_name = "debug_sfgfp_puc19_gibson_v1"

inv = {
    "water"       : "rs17gmh5wafm5p", # catalog; Autoclaved MilliQ H2O; ambient
    "DH5a"        : "rs16pbj944fnny", # catalog; Zymo DH5a; cold_80
    "LB Miller"   : "rs17bafcbmyrmh", # catalog; LB Broth Miller; cold_4
    "Amp 100mgml" : "rs17msfk8ujkca", # catalog; Ampicillin 100mg/ml; cold_20
    "pUC19"       : "rs17tcqmncjfsh", # catalog; pUC19; cold_20
}

# Catalog
transform_plate = p.ref("transform_plate", cont_type="96-pcr", storage="ambient", discard=True)
transform_tube  = transform_plate.well(0)


# ------------------------------------------------------------------------------------
# Plating transformed bacteria according to Tali's protocol (requires different code!)
# http://learn.transcriptic.com/blog/2015/9/9/provisioning-commercial-reagents
# Add 1-5ul plasmid and pre-warm culture plates to 37C before starting.
#

#
# Extra inventory for plating
#
inv["lb-broth-100ug-ml-amp_6-flat"] = "ki17sbb845ssx9" # (kit, not normal ref) from blogpost
inv["noAB-amp_6-flat"] = "ki17reefwqq3sq" # kit id
inv["LB Miller"] = "rs17bafcbmyrmh"

#
# Ampicillin and no ampicillin plates
#
amp_6_flat = Container(None, p.container_type('6-flat'))
p.refs["amp_6_flat"] = Ref('amp_6_flat',
                           {"reserve": inv['lb-broth-100ug-ml-amp_6-flat'], "store": {"where": 'cold_4'}}, amp_6_flat)
noAB_6_flat = Container(None, p.container_type('6-flat'))
p.refs["noAB_6_flat"] = Ref('noAB_6_flat',
                            {"reserve": inv['noAB-amp_6-flat'], "store": {"where": 'cold_4'}}, noAB_6_flat)

#
# Provision competent bacteria
#
p.provision(inv["DH5a"], transform_tube,  ul(50))
p.provision(inv["pUC19"], transform_tube, ul(2))

#
# Heatshock the bacteria to transform using a PCR machine
#
p.seal(transform_plate)
p.thermocycle(transform_plate,
    [{"cycles":  1, "steps": [{"temperature":  "4:celsius", "duration":  "5:minute"}]},
     {"cycles":  1, "steps": [{"temperature": "37:celsius", "duration": "30:minute"}]}],
    volume=ul(50))
p.unseal(transform_plate)

#
# Then dilute bacteria and spread onto 6-flat plates
# Put more on ampicillin plates for more opportunities to get a colony
#
p.provision(inv["LB Miller"], transform_tube, ul(355))
p.mix(transform_tube, ul(150), repetitions=5)
for i in range(6):
    p.spread(transform_tube, amp_6_flat.well(i), ul(55))
    p.spread(transform_tube, noAB_6_flat.well(i), ul(10))

assert transform_tube.volume >= ul(15), transform_tube.volume

#
# Incubate and image 6-flat plates over 18 hours
#
for flat_name, flat in [("amp_6_flat", amp_6_flat), ("noAB_6_flat", noAB_6_flat)]:
    for timepoint in [6,12,18]:
        p.cover(flat)
        p.incubate(flat, "warm_37", "6:hour")
        p.uncover(flat)
        p.image_plate(flat, mode="top", dataref=expid("{}_t{}".format(flat_name, timepoint)))

# ---------------------------------------------------------------
# Analyze protocol
#
jprotocol = json.dumps(p.as_dict(), indent=2)
print(jprotocol)