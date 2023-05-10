import h5py
from os import path
from numpy import array

# Define the tree reader class
class Writer():
    def __init__(self, configuration) -> None:
        """
        Create a tree writer based on hdf5 format
        """
        # Save configuration
        self.configuration = configuration['WRITER']
        # Get path for trees
        self.treePath = self.configuration['retreated_trees_path']
        #
        # Current tree file index
        self.currentTreeFileIndex = 1
        #
        # the hdf5 file
        self.h5f = None

    def write(self, IDsAtSnapshot, snapshots)-> None:
        """
        Save a set of trees in hdf5 format
        """
        # Create and open the file
        filename = path.join(self.treePath, 'tree_{:04d}.h5'.format(self.currentTreeFileIndex))
        self.h5f = h5py.File(filename, 'a')
        #
        # Create a "info" group containing size informations,
        # i.e number of halos at each snapshot
        self.write_infos(IDsAtSnapshot)
        #
        # Write snaphots
        # i.e save all halos at each snaphot in specific data-sets
        self.write_snapshots(snapshots)
        #
        # Close the current file
        self.h5f.close()
        # One more tree file has been generated, update counter
        self.currentTreeFileIndex += 1
    
    def write_infos(self, IDsAtSnapshot)-> None:
        """
        Save snapshot informations
        # Create a "info" group containing size informations,
        # i.e number of halos at each snapshot
        """
        infoGrp = self.h5f.create_group("infos")
        if 1 not in IDsAtSnapshot:
            # For completness reason, add snapshot 1 without any halos
            IDsAtSnapshot[1] = 0
        snapshotContains = array([nhalos for nhalos in list(IDsAtSnapshot.values())], dtype='i8')
        infoGrp.create_dataset("snapshot-contains", data=snapshotContains)
    
    def write_snapshots(self, snapshots)-> None:
        """
        Save all halo ID/properties for each snaphot in a dedicated data-set
        """
        # Group of snaphots
        snapshotsGrp = self.h5f.create_group("snapshots")
        for snapshot, haloList in snapshots.items():
            # Create a sub-group for this specific snaphot
            snapshotGrp = snapshotsGrp.create_group("snapshot_{:03d}".format(snapshot))
            for halo in haloList.values():
                # Create a group for each halo and
                id = halo['IDs']['me']
                haloGrp = snapshotGrp.create_group("halo_{:05d}".format(id))
                #
                # In this specific group create a set of 3 data-set to store:
                # - halos IDs
                # - halos properties
                # - halos infos
                IDs = array([halo['IDs'][key] for key in list(halo['IDs'].keys())], dtype='int8')
                props = array([halo['properties'][key] for key in list(halo['properties'].keys())], dtype='float32')
                infos = array([halo['infos'][key] for key in list(halo['infos'].keys())], dtype='int8')
                # Save
                haloGrp.create_dataset("IDs", data=IDs)
                haloGrp.create_dataset("properties", data=props)
                haloGrp.create_dataset("infos", data=infos)






