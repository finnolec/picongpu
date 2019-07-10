"""
This file is part of the PIConGPU.

Copyright 2017-2019 PIConGPU contributors
Authors: Axel Huebl
License: GPLv3+
"""
from .base_reader import DataReader

import numpy as np
import os


class TransitionRadiationData(DataReader):
    """
    Data Reader for Transition Radiation Plugin.
    """

    def __init__(self, run_directory):
        """
        Parameters
        ----------
        simulation_directory : string
            path to the run directory of PIConGPU
            (the path before ``simOutput/``)
        """
        super().__init__(run_directory)

        self.data_file_prefix = "_transRad_"
        self.data_file_suffix = ".dat"
        self.data = None
        self.omegas = None
        self.thetas = None
        self.phis = None

    def get_data_path(self, species, iteration=None):
        """
        Return the path to the underlying data file.

        Parameters
        ----------
        species : string
            short name of the particle species, e.g. "e" for electrons
            (defined in ``speciesDefinition.param``)
        iteration: (unsignied) int or list of int
            The iteration at which to read the data.

        Returns
        -------
        A string with a path.
        """
        if species is None:
            raise ValueError("The species parameter can not be None!")
        # if iteration is None:
        #    raise ValueError("The iteration can't be None!")


        sim_output_dir = os.path.join(self.run_directory, "simOutput")
        if not os.path.isdir(sim_output_dir):
            raise IOError("The simOutput/ directory does not exist inside "
                          "path:\n  {}\n"
                          "Did you set the proper path to the run directory?\n"
                          "Did the simulation already run?"
                          .format(self.run_directory))

        if iteration is None:
            return sim_output_dir
        else:
            data_file_path = os.path.join(
                sim_output_dir,
                species + self.data_file_prefix + iteration +
                self.data_file_suffix
            )
            if not os.path.isfile(data_file_path):
                raise IOError("The file {} does not exist.\n"
                              "Did the simulation already run?"
                              .format(data_file_path))
            return data_file_path

    def _get_for_iteration(self, species, iteration=None, **kwargs):
        """
        Returns data for transition radiation visualization.

        Parameters
        ----------
        iteration : (unsigned)
            The iteration at which to read the data.
        :param species:
        :param kwargs:
        :return:
        """
        if species is None:
            raise ValueError("The species parameter can not be None!")
        if iteration is None:
            raise ValueError("The iteration can't be None!")

        # Do loading once and store them in ram
        if self.data is None:
            data_file_path = self.get_data_path(species, iteration)

            self.data = np.loadtxt(data_file_path)

            # Read values to automatically create theta, phi and omega arrays
            f = open(self.path + self.trOutputFile)
            tmp = f.readlines()[0]
            f.close()

            tmp = tmp.split("\t")
            if len(tmp) == 1:  # TODO is if necessary
                tmp = str(tmp).split(" ")
                tmp[9] = tmp[9].split("\\")[0]

            self.omegas = np.linspace(tmp[2], tmp[3], tmp[1])
            self.phis = np.linspace(tmp[5], tmp[6], tmp[4])
            self.thetas = np.linspace(tmp[8], tmp[9], tmp[7])

        return self.data

    def get_data(self, species, iteration=None, **kwargs):
        """
        :param observer:  shoudl
        :param n_observer:
        :return:
        """
        if iteration is None:
            raise ValueError("Can't return data for an unknown iteration!")
        if self.data == None:
            self.data = self._get_for_iteration(iteration, species)

        type = kwargs["type"]
        theta = kwargs["theta"]
        phi = kwargs["phi"]
        omega = kwargs["omega"]

        if type is "spectrum":
            # find phi and theta with maximum intensity if they are not given as parameters
            if theta is None and phi is None:
                maxIndex = np.argmax(self.data[, ::])
                theta = np.floor(maxIndex / len(self.phis))
                phi = maxIndex % len(self.phis)
            elif theta is None and phi is not None:
                theta = np.argmax(self.data[phi::len(self.phis), ::])
            elif theta is not None and phi is None:
                phi = np.argmax(self.data[theta * len(self.phis):(theta + 1) * len(self.phis):, ::])
            print("Spectrum is calculated for theta={:.4e} and phi={:.4e}".format(self.thetas[theta],
                                                                                  self.phis[phi]))
            return self.data[theta * len(self.phis) + phi, ::]

    def get_iterations(self, species, species_filter="all"):
        """
        Return an array of iterations with available data.

        Parameters
        ----------
        species : string
            short name of the particle species, e.g. "e" for electrons
            (defined in ``speciesDefinition.param``)
        species_filter: string
            name of the particle species filter, default is "all"
            (defined in ``particleFilters.param``)

        Returns
        -------
        An array with unsigned integers.
        """
        ctr_path = self.get_data_path(species)

        ctr_files = [f for f in os.listdir(ctr_path) if f.endswith(".dat")]
        iters = [int(f.split("_")[2].split(".")[0]) for f in ctr_files]
        return np.array(sorted(iters))
