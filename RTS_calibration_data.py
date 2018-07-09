import numpy as np

class CalGains_data:

    def __init__(self):
        #NOTE: This class can either load a single tile, or all tiles. The parameter single_tile passes in this information, and the dict changes dependent on that.
        self.allx_obs_dict = {}
        self.ally_obs_dict = {}
        self.obs_list = []
        self.problem_obs_list = []
        self.NaN_list = []

    def load_observationsfromFile(self, filename, single_tile=True):

        f = open (filename, "r")

        for line in f:
            values = line.split(',')
            if (values[0] not in self.obs_list):
                self.obs_list.append(values[0])
                if (not single_tile):
                    self.allx_obs_dict[values[0]] = [None]*128
                    self.ally_obs_dict[values[0]] = [None]*128

            if values[2] == "X":
                if (single_tile):
                    self.allx_obs_dict[values[0]] = [float(i) for i in values[3:-1]]
                else:
                    self.allx_obs_dict[values[0]][int(values[1])] = [float(i) for i in values[3:-1]]
            else:
                if (single_tile):
                    self.ally_obs_dict[values[0]] = [float(i) for i in values[3:-1]]
                else:
                    self.ally_obs_dict[values[0]][int(values[1])] = [float(i) for i in values[3:-1]]
        self.obs_list.sort()
        f.close()
