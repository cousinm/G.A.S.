from writer import Writer

# Define the tree reader class
class Builder():
    def __init__(self, configuration) -> None:
        """
        Create a tree builder
        """
        # Save configuration
        self.configuration = configuration['BUILDER']
        # Define a associated writer
        self.writer = Writer(configuration)
        # Reset
        self.reset()
    
    def reset(self) -> None:
        """
        Reset the snapshot to ID dict
        """
        # Link between snapshot and halo current referenced ID
        self.IDsAtSnapshot = {}
    
    def build(self, snapshots) -> bool:
        """
        Build tree indexes for a given set of halo dispatched into a set of snapshots
        """
        #
        # Run throught halo list
        for snapshot, haloList in snapshots.items():
            currentHaloID = -1
            previousHaloID = -1
            idsForSecondLoop = []
            for id, halo in haloList.items():
                #
                # Update halo indexes
                previousHaloID = currentHaloID
                currentHaloID = id
                #
                # Set default values (will be updated next)
                halo['IDs']['myFirstProg'] = -1
                halo['IDs']['myNextProg'] = -1
                halo['infos']['myLevel'] = 1  # I am a main halo per default
                halo['IDs']['myHost'] = -1 # And therefre I do not have an Host
                #
                # Set 'me'
                if snapshot in self.IDsAtSnapshot:
                    me = self.IDsAtSnapshot.get(snapshot) + 1
                else:
                    me = 1
                # Update
                self.IDsAtSnapshot[snapshot] = me
                # set my ID 'me'
                halo['IDs']['me'] = me
                #
                # Set 'myDesc' and 'myNextProg'
                # Set also 'myFirstProg' of my descendent if exist
                # Go to descendent if exist
                if halo['IDs']['desc_id'] > 0:
                    # Get the descendent and set myDesc
                    desc_id = halo['IDs']['desc_id']
                    myDesc = snapshots[snapshot+1][desc_id]['IDs']['me']
                    halo['IDs']['myDesc'] = myDesc
                    #
                    # Set myNextProg of the previous halo if descendent are similar
                    if previousHaloID > 0:
                        if desc_id == snapshots[snapshot][previousHaloID]['IDs']['desc_id']:
                            # I am the next progenitor of the previous halo
                            snapshots[snapshot][previousHaloID]['IDs']['myNextProg'] = me
                    #
                    # Set myFirstProg of my descendent if it is not already done
                    if snapshots[snapshot+1][desc_id]['IDs']['myFirstProg'] < 0:
                        # I am the first one
                        snapshots[snapshot+1][desc_id]['IDs']['myFirstProg'] = me
                else:
                    halo['IDs']['myDesc'] = -1
                #
                # Define sub-structure level and host link
                anIdForSecondLoop = self.build_substructure_link(snapshots, snapshot, id, halo)
                if anIdForSecondLoop > 0:
                    idsForSecondLoop.append(anIdForSecondLoop)
            #
            # Second
            #
            for id in idsForSecondLoop:
                anIdForSecondLoop= self.build_substructure_link(snapshots, snapshot, id, haloList[id])
                if anIdForSecondLoop > 0:
                    msg = 'halo {ih} can not be linked to its main halo'
                    raise IndexError(msg)

            #
        # Save snaphots
        self.writer.write(self.IDsAtSnapshot, snapshots)
        # Reset
        self.reset()
    
    def build_substructure_link(self, snapshots, snapshot, id, halo)-> int:
        """
        Build the substructure link if exist and if currently accessible
        """
        if halo['IDs']['upid'] > 0:
            host_id = halo['IDs']['upid']
            if host_id in snapshots[snapshot]:
                if 'me' in snapshots[snapshot][host_id]['IDs']:
                    myHost = snapshots[snapshot][host_id]['IDs']['me']
                    # I am a sub-structure, and my main halo is accesible 
                    halo['infos']['myLevel'] = 2
                    halo['IDs']['myHost'] = myHost
                    return -1
                else:
                    # My host is not currently build, the link will be created in a second loop
                    return id
            else:
                return -1
        else:
            return -1
