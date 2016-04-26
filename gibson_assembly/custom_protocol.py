from autoprotocol import Unit
from autoprotocol.protocol import Protocol

class CustomProtocol(Protocol):
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
        
        if not dispense_speed and not aspirate_speed:
            if min_volume < Unit.fromstring('100:microliter'):
                dispense_speed = str(min_volume)+'/second'
                aspirate_speed = str(min_volume)+'/second'
                
        if 'mix_before' in mix_kwargs:
            mix_volume_b = mix_kwargs.get("mix_vol_b") or mix_kwargs.get("mix_vol") or volume/2
            if mix_kwargs.get('mix_vol') < Unit.fromstring('100:microliter'):
                mix_vol = str(min_volume)+'/second'            
            
        
        super().transfer(source, dest, volume, one_source=one_source, one_tip=one_tip, 
              aspirate_speed=aspirate_speed, dispense_speed=dispense_speed, 
              aspirate_source=aspirate_source, dispense_target=dispense_target, 
              pre_buffer=pre_buffer, disposal_vol=disposal_vol, 
              transit_vol=transit_vol, blowout_buffer=blowout_buffer, 
              tip_type=tip_type, new_group=new_group,**mix_kwargs)
        
    