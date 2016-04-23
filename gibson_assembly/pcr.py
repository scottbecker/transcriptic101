""" 

Derived from Brian Naughton @ http://blog.booleanbiotech.com/genetic_engineering_pipeline_python.html

"""

from utils import ul

p = Protocol()

# ---------------------------------------------------
# Set up experiment
#
experiment_name = "sfgfp_pcroe_v8"
template_length = 740

_options = {'dilute_primers' : False, # if working stock has not been made
            'dilute_template': False, # if working stock has not been made
            'dilute_dNTP'    : False, # if working stock has not been made
            'run_gel'        : True,  # run a gel to see the plasmid size
            'run_absorbance' : False, # check absorbance at 260/280/320
            'run_sanger'     : False} # sanger sequence the new sequence
options = {k for k,v in _options.items() if v is True}

# ---------------------------------------------------
# Inventory and provisioning
# https://developers.transcriptic.com/v1.0/docs/containers
#
# 'sfgfp2':              'ct17yx8h77tkme', # inventory; sfGFP tube #2, micro-1.5, cold_20
# 'sfgfp_puc19_primer1': 'ct17z9542mrcfv', # inventory; micro-2.0, cold_4
# 'sfgfp_puc19_primer2': 'ct17z9542m5ntb', # inventory; micro-2.0, cold_4
# 'sfgfp_idt_1ngul':     'ct184nnd3rbxfr', # inventory; micro-1.5, cold_4, (ERROR: no template)
#
inv = {
    'Q5 Polymerase':                     'rs16pcce8rdytv', # catalog; Q5 High-Fidelity DNA Polymerase
    'Q5 Buffer':                         'rs16pcce8rmke3', # catalog; Q5 Reaction Buffer
    'dNTP Mixture':                      'rs16pcb542c5rd', # catalog; dNTP Mixture (25mM?)
    'water':                             'rs17gmh5wafm5p', # catalog; Autoclaved MilliQ H2O
    'sfgfp_pcroe_v5_puc19_primer1_10uM': 'ct186cj5cqzjmr', # inventory; micro-1.5, cold_4
    'sfgfp_pcroe_v5_puc19_primer2_10uM': 'ct186cj5cq536x', # inventory; micro-1.5, cold_4
    'sfgfp1':                            'ct17yx8h759dk4', # inventory; sfGFP tube #1, micro-1.5, cold_20
}


# Existing inventory
template_tube = p.ref("sfgfp1", id=inv['sfgfp1'], cont_type="micro-1.5", storage="cold_4").well(0)
dilute_primer_tubes = [p.ref('sfgfp_pcroe_v5_puc19_primer1_10uM', id=inv['sfgfp_pcroe_v5_puc19_primer1_10uM'], cont_type="micro-1.5", storage="cold_4").well(0),
                       p.ref('sfgfp_pcroe_v5_puc19_primer2_10uM', id=inv['sfgfp_pcroe_v5_puc19_primer2_10uM'], cont_type="micro-1.5", storage="cold_4").well(0)]

# New inventory resulting from this experiment
dilute_template_tube = p.ref("sfgfp1_0.25ngul",  cont_type="micro-1.5", storage="cold_4").well(0)
dNTP_10uM_tube       = p.ref("dNTP_10uM",        cont_type="micro-1.5", storage="cold_4").well(0)
sfgfp_pcroe_out_tube = p.ref(expid("amplified"), cont_type="micro-1.5", storage="cold_4").well(0)

# Temporary tubes for use, then discarded
mastermix_tube = p.ref("mastermix", cont_type="micro-1.5", storage="cold_4",  discard=True).well(0)
water_tube =     p.ref("water",     cont_type="micro-1.5", storage="ambient", discard=True).well(0)
pcr_plate =      p.ref("pcr_plate", cont_type="96-pcr",    storage="cold_4",  discard=True)
if 'run_absorbance' in options:
    abs_plate = p.ref("abs_plate", cont_type="96-flat", storage="cold_4", discard=True)

# Initialize all existing inventory
all_inventory_wells = [template_tube] + dilute_primer_tubes
for well in all_inventory_wells:
    init_inventory_well(well)
    print(well.name, well.volume, well.properties)


# -----------------------------------------------------
# Provision water once, for general use
#
p.provision(inv["water"], water_tube, ul(500))

# -----------------------------------------------------
# Dilute primers 1/10 (100uM->10uM) and keep at 4C
#
if 'dilute_primers' in options:
    for primer_num in (0,1):
        p.transfer(water_tube, dilute_primer_tubes[primer_num], ul(90))
        p.transfer(primer_tubes[primer_num], dilute_primer_tubes[primer_num], ul(10), mix_before=True, mix_vol=ul(50))
        p.mix(dilute_primer_tubes[primer_num], volume=ul(50), repetitions=10)


# -----------------------------------------------------
# Dilute template 1/10 (10ng/ul->1ng/ul) and keep at 4C
# OR
# Dilute template 1/40 (10ng/ul->0.25ng/ul) and keep at 4C
#
if 'dilute_template' in options:
    p.transfer(water_tube, dilute_template_tube, ul(195))
    p.mix(dilute_template_tube, volume=ul(100), repetitions=10)

# Dilute dNTP to exactly 10uM
if 'dilute_DNTP' in options:
    p.transfer(water_tube,           dNTP_10uM_tube, ul(6))
    p.provision(inv["dNTP Mixture"], dNTP_10uM_tube, ul(4))


# -----------------------------------------------------
# Q5 PCR protocol
# www.neb.com/protocols/2013/12/13/pcr-using-q5-high-fidelity-dna-polymerase-m0491
#
# 25ul reaction
# -------------
# Q5 reaction buffer      5    ul
# Q5 polymerase           0.25 ul
# 10mM dNTP               0.5  ul -- 1ul = 4x12.5mM
# 10uM primer 1           1.25 ul
# 10uM primer 2           1.25 ul
# 1pg-1ng Template        1    ul -- 0.5 or 1ng/ul concentration
# -------------------------------
# Sum                     9.25 ul
#
#

# Mastermix tube will have 96ul of stuff, leaving space for 4x1ul aliquots of template
p.transfer(water_tube,             mastermix_tube, ul(64))
p.provision(inv["Q5 Buffer"],      mastermix_tube, ul(20))
p.provision(inv['Q5 Polymerase'],  mastermix_tube, ul(1))
p.transfer(dNTP_10uM_tube,         mastermix_tube, ul(1), mix_before=True, mix_vol=ul(2))
p.transfer(dilute_primer_tubes[0], mastermix_tube, ul(5), mix_before=True, mix_vol=ul(10))
p.transfer(dilute_primer_tubes[1], mastermix_tube, ul(5), mix_before=True, mix_vol=ul(10))
p.mix(mastermix_tube, volume="48:microliter", repetitions=10)

# Transfer mastermix to pcr_plate without template
p.transfer(mastermix_tube, pcr_plate.wells(["A1","B1","C1"]), ul(24))
p.transfer(mastermix_tube, pcr_plate.wells(["A2"]),           ul(24)) # acknowledged dead volume problems
p.mix(pcr_plate.wells(["A1","B1","C1","A2"]), volume=ul(12), repetitions=10)

# Finally add template
p.transfer(template_tube,  pcr_plate.wells(["A1","B1","C1"]), ul(1))
p.mix(pcr_plate.wells(["A1","B1","C1"]), volume=ul(12.5), repetitions=10)

# ---------------------------------------------------------
# Thermocycle with Q5 and hot start
# 61.1 annealing temperature is recommended by NEB protocol
# p.seal is enforced by transcriptic
#
extension_time = int(max(2, np.ceil(template_length * (11.0/1000))))
assert 0 < extension_time < 60, "extension time should be reasonable for PCR"

cycles = [{"cycles":  1, "steps": [{"temperature": "98:celsius", "duration": "30:second"}]}] + \
    touchdown(70, 61, [8, 25, extension_time], stepsize=0.5) + \
    [{"cycles": 16, "steps": [{"temperature": "98:celsius", "duration": "8:second"},
                              {"temperature": "61.1:celsius", "duration": "25:second"},
                              {"temperature": "72:celsius", "duration": "{:d}:second".format(extension_time)}]},
     {"cycles":  1, "steps": [{"temperature": "72:celsius", "duration": "2:minute"}]}]
p.seal(pcr_plate)
p.thermocycle(pcr_plate, cycles, volume=ul(25))

# --------------------------------------------------------
# Run a gel to hopefully see a 740bp fragment
#
if 'run_gel' in options:
    p.unseal(pcr_plate)
    p.mix(pcr_plate.wells(["A1","B1","C1","A2"]), volume=ul(12.5), repetitions=10)
    p.transfer(pcr_plate.wells(["A1","B1","C1","A2"]), pcr_plate.wells(["D1","E1","F1","D2"]),
               [ul(2), ul(4), ul(8), ul(8)])
    p.transfer(water_tube, pcr_plate.wells(["D1","E1","F1","D2"]),
               [ul(18),ul(16),ul(12),ul(12)], mix_after=True, mix_vol=ul(10))
    p.gel_separate(pcr_plate.wells(["D1","E1","F1","D2"]),
                   ul(20), "agarose(10,2%)", "ladder1", "10:minute", expid("gel"))

#---------------------------------------------------------
# Absorbance dilution series. Take 1ul out of the 25ul pcr plate wells
#
if 'run_absorbance' in options:
    p.unseal(pcr_plate)
    abs_wells = ["A1","B1","C1","A2","B2","C2","A3","B3","C3"]

    p.transfer(water_tube, abs_plate.wells(abs_wells[0:6]), ul(10))
    p.transfer(water_tube, abs_plate.wells(abs_wells[6:9]), ul(9))

    p.transfer(pcr_plate.wells(["A1","B1","C1"]), abs_plate.wells(["A1","B1","C1"]), ul(1), mix_after=True, mix_vol=ul(5))
    p.transfer(abs_plate.wells(["A1","B1","C1"]), abs_plate.wells(["A2","B2","C2"]), ul(1), mix_after=True, mix_vol=ul(5))
    p.transfer(abs_plate.wells(["A2","B2","C2"]), abs_plate.wells(["A3","B3","C3"]), ul(1), mix_after=True, mix_vol=ul(5))
    for wavelength in [260, 280, 320]:
        p.absorbance(abs_plate, abs_plate.wells(abs_wells),
                     "{}:nanometer".format(wavelength), exp_id("abs_{}".format(wavelength)), flashes=25)

# -----------------------------------------------------------------------------
# Sanger sequencing: https://developers.transcriptic.com/docs/sanger-sequencing
# "Each reaction should have a total volume of 15 ul and we recommend the following composition of DNA and primer:
#  PCR product (40 ng), primer (1 ul of a 10 µM stock)"
#
#  By comparing to the gel ladder concentration (175ng/lane), it looks like 5ul of PCR product has approximately 30ng of DNA
#
if 'run_sanger' in options:
    p.unseal(pcr_plate)
    seq_wells = ["G1","G2"]
    for primer_num, seq_well in [(0, seq_wells[0]),(1, seq_wells[1])]:
        p.transfer(dilute_primer_tubes[primer_num], pcr_plate.wells([seq_well]),
                   ul(1), mix_before=True, mix_vol=ul(50))
        p.transfer(pcr_plate.wells(["A1"]), pcr_plate.wells([seq_well]),
                   ul(5), mix_before=True, mix_vol=ul(10))
        p.transfer(water_tube, pcr_plate.wells([seq_well]), ul(9))

    p.mix(pcr_plate.wells(seq_wells), volume=ul(7.5), repetitions=10)
    p.sangerseq(pcr_plate, pcr_plate.wells(seq_wells[0]).indices(), expid("seq1"))
    p.sangerseq(pcr_plate, pcr_plate.wells(seq_wells[1]).indices(), expid("seq2"))

# -------------------------------------------------------------------------
# Then consolidate to one tube. Leave at least 3ul dead volume in each tube
#
remaining_volumes = [well.volume - dead_volume['96-pcr'] for well in pcr_plate.wells(["A1","B1","C1"])]
print("Consolidated volume", sum(remaining_volumes, ul(0)))
p.consolidate(pcr_plate.wells(["A1","B1","C1"]), sfgfp_pcroe_out_tube, remaining_volumes, allow_carryover=True)

uprint("\nProtocol 1. Amplify the insert (oligos previously synthesized)")
jprotocol = json.dumps(p.as_dict(), indent=2)
!echo '{jprotocol}' | transcriptic analyze
open("protocol_{}.json".format(experiment_name),'w').write(jprotocol)