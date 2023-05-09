from os import path
from itertools import product
from builder import Builder

# Define the tree reader class
class Reader():
    def __init__(self, configuration) -> None:
        """
        Create a reader
        """
        # Create a builder
        self.builder = Builder(configuration)
        # Save configuration
        self.configuration = configuration['READER']
        # Get max number of halos expected
        self.nbMaxHalo = int(self.configuration['nb_halos_treated'])
        self.DnbMaxHalo = int(self.configuration['Dnb_halos_treated'])
        # Cosmological parameters
        self.cosmo = {}
        # Box size
        self.BoxSize = 0.
        # Number of tree contains in the current file
        self.nbTrees = 0
        # init
        self.reset()

    def reset(self)-> None:
        """
        Reset data
        """
        # Snapshots
        self.currentTreeID = -1
        self.snapshots = {}
        # Infos
        self.infos = {'nbTrees': 0, 'nbHalos': 0}
    
    def get_tree_file_list(self)-> None:
        """
        Get tree file list
        """
        treeFileList = self.configuration['tree_file_list']
        # Switch to specific case
        self.treeFileList = []
        if isinstance(treeFileList, str):
            if treeFileList in ['all']:
                # Build the complete list of tree files
                indexes = product([0, 1, 2, 3, 4], repeat=3)
                for t in list(indexes):
                    treeFile = 'tree_{}_{}_{}.dat'.format(t[0], t[1], t[2])
                    self.treeFileList.append(treeFile)
            else:
                indexes = treeFileList.split(', ')
                for ind in indexes:
                    treeFile = 'tree_{}.dat'.format(ind)
                    self.treeFileList.append(treeFile)
    
    def run(self)-> None :
        """
        Loop over tree file and load each one
        """
        # Get tree file list
        self.get_tree_file_list()
        #
        treePath = self.configuration['original_trees_path']
        for treeFile in self.treeFileList:
            filepath = path.join(treePath, treeFile)
            self.read(filepath)

    def read(self, filepath):
        """
        Read a treefile
        """
        # File path of the tree file
        self.filepath = filepath
        #
        with open(self.filepath) as myTreeFile:
            for il, line in enumerate(myTreeFile):
                if il == 0:
                    # Read line information labels
                    infos = line[1:].strip().split(' ')
                    #
                elif il in [1, 2]:
                    # Cosmological parameters and box size
                    continue
                    #
                elif line.startswith('#tree'):
                    # A new tree start
                    #
                    # BUILDER
                    #
                    if (self.infos['nbHalos'] > self.nbMaxHalo - self.DnbMaxHalo):
                        # Build and save the current set of trees
                        self.builder.build(self.snapshots)
                        # Reset current data
                        self.reset()
                    #
                    self.currentTreeID = int(line.split(' ')[1])  # The tree ID is the z=0 corresponding halo ID
                    # One more tree will be treated
                    self.infos['nbTrees'] += 1
                    #
                elif line.startswith('#'):
                    continue  # it is just a comment
                else:
                    data = line.strip().split()
                    if len(data) == 1:
                        # Save the number of tree contains in the file
                        self.nbTrees = int(line)
                    else:
                        # Load halo data
                        halo = {'IDs': {}, 'infos': {}, 'properties': {}}
                        for i, info in enumerate(infos):
                            key = info.split('(')[0]
                            if key in ['num_prog', 'Snap_num', 'mmp?', 'phantom']:
                                # Infos
                                halo['infos'][key] = int(data[i])
                            elif (key.find('id') < 0 and key.find('ID') < 0
                                or (key in ['Tidal_Force'])):
                                # Properties
                                halo['properties'][key] = float(data[i])
                            else:
                                # IDs
                                halo['IDs'][key] = int(data[i])
                        id = halo['IDs']['id']
                        snapshot = halo['infos']['Snap_num']
                        if snapshot not in self.snapshots:
                            self.snapshots[snapshot] = {id: halo}
                        else:
                            self.snapshots[snapshot][id] = halo
                        # One more halo load
                        # Update the number of total halo treated
                        self.infos['nbHalos'] += 1
            #
            # End of loading
            print('{} / {} have been loaded'.format(self.infos['nbTrees'], self.nbTrees))



