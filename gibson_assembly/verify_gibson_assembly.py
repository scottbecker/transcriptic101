"""Debugging transformation protocol: Gibson assembly followed by qPCR and a gel
v2: include v3 Gibson assembly"""

p = Protocol()
options = {}

experiment_name = "debug_sfgfp_puc19_gibson_seq_v2"

inv = {
    "water"                       : "rs17gmh5wafm5p", # catalog; Autoclaved MilliQ H2O; ambient
    "M13_F"                       : "rs17tcpqwqcaxe", # catalog; M13 Forward (-41); cold_20 (1ul = 100pmol)
    "M13_R"                       : "rs17tcph6e2qzh", # catalog; M13 Reverse (-48); cold_20 (1ul = 100pmol)
    "SensiFAST_SYBR_No-ROX"       : "rs17knkh7526ha", # catalog; SensiFAST SYBR for qPCR
    "sfgfp_puc19_gibson_v1_clone" : "ct187rzdq9kd7q", # inventory; assembled sfGFP; cold_4
    "sfgfp_puc19_gibson_v3_clone" : "ct188ejywa8jcv", # inventory; assembled sfGFP; cold_4
}

# ---------------------------------------------------------------
# First get my sfGFP pUC19 clones, assembled with Gibson assembly
#
clone_plate1 = p.ref("sfgfp_puc19_gibson_v1_clone", id=inv["sfgfp_puc19_gibson_v1_clone"],
                     cont_type="96-pcr", storage="cold_4", discard=False)
clone_plate2 = p.ref("sfgfp_puc19_gibson_v3_clone", id=inv["sfgfp_puc19_gibson_v3_clone"],
                     cont_type="96-pcr", storage="cold_4", discard=False)

water_tube = p.ref("water", cont_type="micro-1.5", storage="cold_4", discard=True).well(0)
master_tube = p.ref("master", cont_type="micro-1.5", storage="cold_4", discard=True).well(0)
primer_tube = p.ref("primer", cont_type="micro-1.5", storage="cold_4", discard=True).well(0)

pcr_plate = p.ref(expid("pcr_plate"), cont_type="96-pcr", storage="cold_4", discard=False)

init_inventory_well(clone_plate1.well("A1"))
init_inventory_well(clone_plate2.well("A1"))

seq_wells = ["B2","B4","B6", # clone_plate1
             "D2","D4","D6", # clone_plate2
             "F2","F4"] # control

# clone_plate2 was diluted 4X (20ul->80ul), according to NEB instructions
assert clone_plate1.well("A1").volume == ul(18), clone_plate1.well("A1").volume
assert clone_plate2.well("A1").volume == ul(78), clone_plate2.well("A1").volume

# --------------------------------------------------------------
# Provisioning
#
p.provision(inv["water"], water_tube, ul(500))

# primers, diluted 2X, discarded at the end
p.provision(inv["M13_F"], primer_tube, ul(13))
p.provision(inv["M13_R"], primer_tube, ul(13))
p.transfer(water_tube, primer_tube, ul(26), mix_after=True, mix_vol=ul(20), repetitions=10)

# -------------------------------------------------------------------
# PCR Master mix -- 10ul SYBR mix, plus 1ul each undiluted primer DNA (100pmol)
# Also add 15ul of dead volume
#
p.provision(inv['SensiFAST_SYBR_No-ROX'], master_tube, ul(11+len(seq_wells)*10))
p.transfer(primer_tube, master_tube, ul(4+len(seq_wells)*4))
p.mix(master_tube, volume=ul(63), repetitions=10)
assert master_tube.volume == ul(127) # 15ul dead volume

p.distribute(master_tube, pcr_plate.wells(seq_wells), ul(14), allow_carryover=True)
p.distribute(water_tube, pcr_plate.wells(seq_wells),
             [ul(ul) for ul in [5,4,2, 4,2,0, 6,6]],
             allow_carryover=True)

# Template -- starting with some small, unknown amount of DNA produced by Gibson
p.transfer(clone_plate1.well("A1"), pcr_plate.wells(seq_wells[0:3]), [ul(1),ul(2),ul(4)], one_tip=True)
p.transfer(clone_plate2.well("A1"), pcr_plate.wells(seq_wells[3:6]), [ul(2),ul(4),ul(6)], one_tip=True)

assert all(pcr_plate.well(w).volume == ul(20) for w in seq_wells)
assert clone_plate1.well("A1").volume == ul(11)
assert clone_plate2.well("A1").volume == ul(66)

# --------------------------------------------------------------
# qPCR
# standard melting curve parameters
#
p.seal(pcr_plate)
p.thermocycle(pcr_plate, [{"cycles":  1, "steps": [{"temperature": "95:celsius","duration": "2:minute"}]},
                          {"cycles": 40, "steps": [{"temperature": "95:celsius","duration": "5:second"},
                                                   {"temperature": "60:celsius","duration": "20:second"},
                                                   {"temperature": "72:celsius","duration": "15:second", "read": True}]}],
    volume=ul(20), # volume is optional
    dataref=expid("qpcr"),
    dyes={"SYBR": seq_wells}, # dye must be specified (tells transcriptic what aborbance to use?)
    melting_start="65:celsius", melting_end="95:celsius", melting_increment="0.5:celsius", melting_rate="5:second")


# --------------------------------------------------------------
# Gel -- 20ul required
# Dilute such that I have 11ul for sequencing
#
p.unseal(pcr_plate)
p.distribute(water_tube, pcr_plate.wells(seq_wells), ul(11))
p.gel_separate(pcr_plate.wells(seq_wells), ul(20), "agarose(8,0.8%)", "ladder1", "10:minute", expid("gel"))

# This appears to be a bug in Transcriptic. The actual volume should be 11ul
# but it is not updating after running a gel with 20ul.
# Primer tube should be equal to dead volume, or it's a waste
assert all(pcr_plate.well(w).volume==ul(31) for w in seq_wells)
assert primer_tube.volume == ul(16) == dead_volume['micro-1.5'] + ul(1)
assert water_tube.volume > ul(25)

# ---------------------------------------------------------------
# Test and run protocol
#
jprotocol = json.dumps(p.as_dict(), indent=2)
print(jprotocol)