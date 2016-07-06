from autoprotocol.unit import Unit

def sample_protocol(protocol, params):
    dest_plate = params["destination_plate"]
    wells_to_measure = []

    for location in params["dye_locations"]:
        protocol.transfer(params["dye"],
                          dest_plate.well(location["well_index"]),
                          location["volume"])
        if location["volume"] != Unit(100, "microliter"):
            protocol.transfer(params["water"],
                              dest_plate.well(location["well_index"]),
                              Unit(100,"microliter") - location["volume"],
                              mix_after = True)
        wells_to_measure.append(location["well_index"])

    protocol.absorbance(dest_plate, wells_to_measure, "475:nanometer", "test")


if __name__ == '__main__':
    from autoprotocol.harness import run
    run(sample_protocol, "SampleProtocol2")
