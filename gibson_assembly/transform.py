"""Full Gibson assembly and transformation protocol for sfGFP and pUC19
v1: Spread IPTG and X-gal onto plates, then spread cells
v2: Mix IPTG, X-gal and cells; spread the mixture
v3: exclude X-gal so I can do colony picking better
v4: repeat v3 to try other excitation/emission wavelengths"""

p = Protocol()


do_transformation()
measure_growth()
if 'sanger' in options: do_sanger_seq()


# ---------------------------------------------------------------
# Output protocol
#
jprotocol = json.dumps(p.as_dict(), indent=2)
!echo '{jprotocol}' | transcriptic analyze

#print("\nProtocol {}\n\n{}".format(experiment_name, jprotocol))
open("protocol_{}.json".format(experiment_name),'w').write(jprotocol)