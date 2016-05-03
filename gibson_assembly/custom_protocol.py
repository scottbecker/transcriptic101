from autoprotocol import Unit
from autoprotocol.protocol import Protocol
from utils import ul

class CustomProtocol(Protocol):
   
    def __init__(self):
        #catalog inventory
        self.inv = {
            "water"       : "rs17gmh5wafm5p", # catalog; Autoclaved MilliQ H2O; ambient
            "DH5a"        : "rs16pbj944fnny", # catalog; Zymo DH5a; cold_80
            "Gibson Mix"  : "rs16pfatkggmk5", # catalog; Gibson Mix (2X); cold_20
            "LB Miller"   : "rs17bafcbmyrmh", # catalog; LB Broth Miller; cold_4
            "Amp 100mgml" : "rs17msfk8ujkca", # catalog; Ampicillin 100mg/ml; cold_20
            "pUC19"       : "rs17tcqmncjfsh", # catalog; pUC19; cold_20            
            "M13_F"                       : "rs17tcpqwqcaxe", # catalog; M13 Forward (-41); cold_20 (1ul = 100pmol)
            "M13_R"                       : "rs17tcph6e2qzh", # catalog; M13 Reverse (-48); cold_20 (1ul = 100pmol)
            "SensiFAST_SYBR_No-ROX"       : "rs17knkh7526ha", # catalog; SensiFAST SYBR for qPCR      
            # kits (must be used differently)
            "lb-broth-100ug-ml-amp_6-flat" : "ki17sbb845ssx9", # catalog; ampicillin plates
            "noAB-amp_6-flat" : "ki17reefwqq3sq" # catalog; no antibiotic plates            
        }
        
        super(CustomProtocol,self).__init__()
       
    def ref(self, name, id=None, cont_type=None, storage=None, discard=None):
        """
        
        Add a Ref object to the dictionary of Refs associated with this protocol
        and return a Container with the id, container type and storage or
        discard conditions specified.
        
        This custom version allows discard=True to overide storage
    
        Example Usage:
    
        .. code-block:: python
    
            p = Protocol()
    
            # ref a new container (no id specified)
            sample_ref_1 = p.ref("sample_plate_1",
                                 cont_type="96-pcr",
                                 discard=True)
    
            # ref an existing container with a known id
            sample_ref_2 = p.ref("sample_plate_2",
                                 id="ct1cxae33lkj",
                                 cont_type="96-pcr",
                                 storage="ambient")
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            {
              "refs": {
                "sample_plate_1": {
                  "new": "96-pcr",
                  "discard": true
                },
                "sample_plate_2": {
                  "id": "ct1cxae33lkj",
                  "store": {
                    "where": "ambient"
                  }
                }
              },
              "instructions": []
            }
    
        Parameters
        ----------
        name : str
            name of the container/ref being created.
        id : str
            id of the container being created, from your organization's
            inventory on http://secure.transcriptic.com.  Strings representing
            ids begin with "ct".
        cont_type : str, ContainerType
            container type of the Container object that will be generated.
        storage : {"ambient", "cold_20", "cold_4", "warm_37"}, optional
            temperature the container being referenced should be stored at
            after a run is completed.  Either a storage condition must be
            specified or discard must be set to True.
        discard : bool, optional
            if no storage condition is specified and discard is set to True,
            the container being referenced will be discarded after a run.
    
        Returns
        -------
        container : Container
            Container object generated from the id and container type
             provided.
    
        Raises
        ------
        RuntimeError
            If a container previously referenced in this protocol (existant
            in refs section) has the same name as the one specified.
        RuntimeError
            If no container type is specified.
        RuntimeError
            If no valid storage or discard condition is specified.
    
        """    
        
        if discard:
            storage=None
        
        return super(CustomProtocol,self).ref(name, id=id, cont_type=cont_type, 
                                              storage=storage, discard=discard)
    
    
    def transfer(self, source, dest, volume, one_source=False, one_tip=False, 
                aspirate_speed=None, dispense_speed=None, 
                aspirate_source=None, dispense_target=None, 
                pre_buffer=None, disposal_vol=None, 
                transit_vol=None, blowout_buffer=None, 
                tip_type=None, new_group=False, **mix_kwargs):
        """
        Transfer liquid from one specific well to another.  A new pipette tip
        is used between each transfer step unless the "one_tip" parameter
        is set to True.

        Example Usage:

        .. code-block:: python

            p = Protocol()
            sample_plate = p.ref("sample_plate",
                                 ct32kj234l21g,
                                 "96-flat",
                                 storage="warm_37")


            # a basic one-to-one transfer:
            p.transfer(sample_plate.well("B3"),
                       sample_plate.well("C3"),
                       "20:microliter")

            # using a basic transfer in a loop:
            for i in xrange(1, 12):
              p.transfer(sample_plate.well(i-1),
                         sample_plate.well(i),
                         "10:microliter")

            # transfer liquid from each well in the first column of a 96-well
            # plate to each well of the second column using a new tip and
            # a different volume each time:
            volumes = ["5:microliter", "10:microliter", "15:microliter",
                       "20:microliter", "25:microliter", "30:microliter",
                       "35:microliter", "40:microliter"]

            p.transfer(sample_plate.wells_from(0,8,columnwise=True),
                       sample_plate.wells_from(1,8,columnwise=True),
                       volumes)

            # transfer liquid from wells A1 and A2 (which both contain the same
            # source) into each of the following 10 wells:
            p.transfer(sample_plate.wells_from("A1", 2),
                       sample_plate.wells_from("A3", 10),
                       "10:microliter",
                       one_source=True)

            # transfer liquid from wells containing the same source to multiple
            # other wells without discarding the tip in between:
            p.transfer(sample_plate.wells_from("A1", 2),
                       sample_plate.wells_from("A3", 10),
                       "10:microliter",
                       one_source=True,
                       one_tip=True)


        Parameters
        ----------
        source : Well, WellGroup
            Well or wells to transfer liquid from.  If multiple source wells
            are supplied and one_source is set to True, liquid will be
            transfered from each source well specified as long as it contains
            sufficient volume. Otherwise, the number of source wells specified
            must match the number of destination wells specified and liquid
            will be transfered from each source well to its corresponding
            destination well.
        dest : Well, WellGroup
            Well or WellGroup to which to transfer liquid.  The number of
            destination wells must match the number of source wells specified
            unless one_source is set to True.
        volume : str, Unit, list
            The volume(s) of liquid to be transferred from source wells to
            destination wells.  Volume can be specified as a single string or
            Unit, or can be given as a list of volumes.  The length of a list
            of volumes must match the number of destination wells given unless
            the same volume is to be transferred to each destination well.
        one_source : bool, optional
            Specify whether liquid is to be transferred to destination wells
            from a group of wells all containing the same substance.
        one_tip : bool, optional
            Specify whether all transfer steps will use the same tip or not.
        mix_after : bool, optional
            Specify whether to mix the liquid in the destination well after
            liquid is transferred.
        mix_before : bool, optional
            Specify whether to mix the liquid in the source well before
            liquid is transferred.
        mix_vol : str, Unit, optional
            Volume to aspirate and dispense in order to mix liquid in a wells
            before and/or after each transfer step.
        repetitions : int, optional
            Number of times to aspirate and dispense in order to mix
            liquid in well before and/or after each transfer step.
        flowrate : str, Unit, optional
            Speed at which to mix liquid in well before and/or after each
            transfer step.
        aspirate speed : str, Unit, optional
            Speed at which to aspirate liquid from source well.  May not be
            specified if aspirate_source is also specified. By default this is
            the maximum aspiration speed, with the start speed being half of
            the speed specified.
        dispense_speed : str, Unit, optional
            Speed at which to dispense liquid into the destination well.  May
            not be specified if dispense_target is also specified.
        aspirate_source : fn, optional
            Can't be specified if aspirate_speed is also specified.
        dispense_target : fn, optional
            Same but opposite of  aspirate_source.
        pre_buffer : str, Unit, optional
            Volume of air aspirated before aspirating liquid.
        disposal_vol : str, Unit, optional
            Volume of extra liquid to aspirate that will be dispensed into
            trash afterwards.
        transit_vol : str, Unit, optional
            Volume of air aspirated after aspirating liquid to reduce presence
            of bubbles at pipette tip.
        blowout_buffer : bool, optional
            If true the operation will dispense the pre_buffer along with the
            dispense volume. Cannot be true if disposal_vol is specified.
        tip_type : str, optional
            Type of tip to be used for the transfer operation.
        new_group : bool, optional

        Raises
        ------
        RuntimeError
            If more than one volume is specified as a list but the list length
            does not match the number of destination wells given.
        RuntimeError
            If transferring from WellGroup to WellGroup that have different
            number of wells and one_source is not True.

        """
        
        if type(volume) == list:
            min_volume = min(volume)
        else:
            min_volume = volume        
        
        if min_volume < ul(10) and not mix_kwargs.get('mix_after') \
           and not mix_kwargs.get('ignore_mix_after_warning'):
            raise Exception('mix_after required for <10uL of solution to ensure complete transfer. \n'
                            'Ensure you have are pipetting into something big enough and set this')
            
        if 'ignore_mix_after_warning' in mix_kwargs:
            del mix_kwargs['ignore_mix_after_warning']
            
        super().transfer(source, dest, volume, one_source=one_source, one_tip=one_tip, 
              aspirate_speed=aspirate_speed, dispense_speed=dispense_speed, 
              aspirate_source=aspirate_source, dispense_target=dispense_target, 
              pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
              transit_vol=transit_vol, blowout_buffer=blowout_buffer, 
              tip_type=tip_type, new_group=new_group,**mix_kwargs)
        
    def distribute(self, source, dest, volume, allow_carryover=False,
                   mix_before=False, mix_after=False, mix_vol=None, repetitions=10,
                   flowrate="100:microliter/second", aspirate_speed=None,
                   aspirate_source=None, distribute_target=None,
                   pre_buffer=None, disposal_vol=None, transit_vol=None,
                   blowout_buffer=None, tip_type=None, new_group=False,
                   ignore_mix_after_warning=False):
        """
        Distribute liquid from source well(s) to destination wells(s).
    
        For volumes under 10uL, multiple transfer operations will be used instead
    
    
        Example Usage:
    
        .. code-block:: python
    
            p = Protocol()
            sample_plate = p.ref("sample_plate",
                                 None,
                                 "96-flat",
                                 storage="warm_37")
            sample_source = p.ref("sample_source",
                                  "ct32kj234l21g",
                                  "micro-1.5",
                                  storage="cold_20")
    
            p.distribute(sample_source.well(0),
                         sample_plate.wells_from(0,8,columnwise=True),
                         "200:microliter",
                         mix_before=True,
                         mix_vol="500:microliter",
                         repetitions=20)
    
        Autoprotocol Output:
    
        .. code-block:: json
    
            "instructions": [
              {
                "groups": [
                  {
                    "distribute": {
                      "to": [
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/0"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/12"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/24"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/36"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/48"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/60"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/72"
                        },
                        {
                          "volume": "150.0:microliter",
                          "well": "sample_plate/84"
                        }
                      ],
                      "from": "sample_source/0",
                      "mix_before": {
                        "volume": "500:microliter",
                        "repetitions": 20,
                        "speed": "100:microliter/second"
                      }
                    }
                  }
                ],
                "op": "pipette"
              }
            ]
    
        Parameters
        ----------
        source : Well, WellGroup
            Well or wells to distribute liquid from.  If passed as a WellGroup
            with set_volume() called on it, liquid will be automatically be
            drawn from the wells specified using the fill_wells function.
        dest : Well, WellGroup
            Well or wells to distribute liquid to.
        volume : str, Unit, list
            Volume of liquid to be distributed to each destination well.  If a
            single string or unit is passed to represent the volume, that
            volume will be distributed to each destination well.  If a list of
            volumes is provided, that volume will be distributed to the
            corresponding well in the WellGroup provided. The length of the
            volumes list must therefore match the number of wells in the
            destination WellGroup if destination wells are recieving different
            volumes.
        allow_carryover : bool, optional
            specify whether the same pipette tip can be used to aspirate more
            liquid from source wells after the previous volume aspirated has
            been depleted.
        mix_before : bool, optional
            Specify whether to mix the liquid in the destination well before
            liquid is transferred.
        mix_after : bool, optional
            Specify whether to mix the liquid in the destination well after
            liquid is transferred.
        mix_vol : str, Unit, optional
            Volume to aspirate and dispense in order to mix liquid in a wells
            before liquid is distributed.
        repetitions : int, optional
            Number of times to aspirate and dispense in order to mix
            liquid in a well before liquid is distributed.
        flowrate : str, Unit, optional
            Speed at which to mix liquid in well before liquid is distributed.
        aspirate speed : str, Unit, optional
            Speed at which to aspirate liquid from source well.  May not be
            specified if aspirate_source is also specified. By default this is
            the maximum aspiration speed, with the start speed being half of
            the speed specified.
        aspirate_source : fn, optional
            Can't be specified if aspirate_speed is also specified.
        distribute_target : fn, optional
            A function that contains additional parameters for distributing to
            target wells including depth, dispense_speed, and calibrated
            volume.
            If this parameter is specified, the same parameters will be
            applied to every destination well.
            Can't be specified if dispense_speed is also specified.
        pre_buffer : str, Unit, optional
            Volume of air aspirated before aspirating liquid.
        disposal_vol : str, Unit, optional
            Volume of extra liquid to aspirate that will be dispensed into
            trash afterwards.
        transit_vol : str, Unit, optional
            Volume of air aspirated after aspirating liquid to reduce presence
            of bubbles at pipette tip.
        blowout_buffer : bool, optional
            If true the operation will dispense the pre_buffer along with the
            dispense volume.
            Cannot be true if disposal_vol is specified.
    
        Raises
        ------
        RuntimeError
            If no mix volume is specified for the mix_before instruction.
        ValueError
            If source and destination well(s) is/are not expressed as either
            Wells or WellGroups.
    
        """    
        

        if volume < ul(10):
            for dest_well in dest:
                self.transfer(source,dest_well,volume,mix_before=mix_before,
                              mix_after=mix_after,mix_vol=mix_vol,
                              aspirate_speed=aspirate_speed, 
                              aspirate_source=aspirate_source, 
                              pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
                              transit_vol=transit_vol, blowout_buffer=blowout_buffer,
                              tip_type=tip_type, new_group=new_group,
                              ignore_mix_after_warning=ignore_mix_after_warning)                               
        else:
            super().distribute(source, dest, volume, allow_carryover=False,
                               mix_before=False, mix_vol=None, repetitions=10,
                               flowrate="100:microliter/second", aspirate_speed=None,
                               aspirate_source=None, distribute_target=None,
                               pre_buffer=None, disposal_vol=None, transit_vol=None,
                               blowout_buffer=None, tip_type=None, new_group=False)