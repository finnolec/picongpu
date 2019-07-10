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
        self.data_file_folder = "transRad/"
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
                self.data_file_folder + species + self.data_file_prefix +
                str(iteration) + self.data_file_suffix
            )
            if not os.path.isfile(data_file_path):
                raise IOError("The file {} does not exist.\n"
                              "Did the simulation already run?"
                              .format(data_file_path))
            return data_file_path

    def _get_for_iteration(self, iteration=None, **kwargs):
        """
        Returns data for transition radiation visualization as specified with "type" for
        transition_radiation_visualizer.py

        Parameters
        ----------
        iteration : (unsigned) int
            The iteration at which to read the data.

        :return:
        """
        species = kwargs["species"]
        if species is None:
            raise ValueError("The species parameter can not be None!")
        kwargs.pop("species")  # remove entry from kwargs because we pass it on

        if iteration is None:
            raise ValueError("The iteration can't be None!")

        # Do loading once and store them in ram
        if self.data is None:
            data_file_path = self.get_data_path(species=species, iteration=iteration)

            self.data = np.loadtxt(data_file_path)

            # Read values to automatically create theta, phi and omega arrays
            f = open(data_file_path)
            parameters = f.readlines()[0].split("\t")
            f.close()

            # Create discretized arrays or angles and frequency as they are discretized for the
            # calculation in PIConGPU. This is necessary for the labels for the axes.
            self.omegas = np.logspace(np.log10(float(parameters[2])), np.log10(float(parameters[3])),
                                      int(parameters[1]))
            self.phis = np.linspace(float(parameters[5]), float(parameters[6]), int(parameters[4]))
            self.thetas = np.linspace(float(parameters[8]), float(parameters[9]), int(parameters[7]))

        return self.get_data(species=species, iteration=iteration, **kwargs)

    def get_data(self, species, iteration=None, **kwargs):
        """
        Calculates data as specified with "type" for the plot from transition_radiation_visualizer.py.

        Parameters
        ----------
        kwargs: dictionary with further keyword arguments, valid are:
            species: string
                short name of the particle species, e.g. 'e' for electrons
                (defined in ``speciesDefinition.param``)
            iteration: int
                number of the iteration
            time: float
                simulation time.
                Only one of 'iteration' or 'time' should be passed!
            type: string
                name of figure type. valid figure types are:
                    'spectrum' - (default) plots transition radiation spectrum
                        at angles theta and phi over the frequency omega
            phi: int
                index of polar angle for a fixed value
            theta: int
                index of azimuth angle for a fixed value
            omega: int
                index of frequency for a fixed value, pointless in a spectrum

        Returns
        ----------
        A tuple of the x and y values for the plot as one dimensional arrays for normal figures
        (not heatmaps) and as colormeshes for heatmaps.
        """
        if iteration is None:
            raise ValueError("Can't return data for an unknown iteration!")
        if self.data is None:
            self.data = self._get_for_iteration(species, iteration, **kwargs)

        # Specify plot type
        type = kwargs["type"]

        # Load fixed values for observation angles from arguments, if given
        theta = kwargs["theta"]
        phi = kwargs["phi"]
        omega = kwargs["omega"]

        if type is "spectrum":
            # find phi and theta with maximum intensity if they are not given as parameters
            if theta is None and phi is None:
                maxIndex = np.argmax(self.data[:, 0])
                theta = int(np.floor(maxIndex / len(self.phis)))
                phi = maxIndex % len(self.phis)
            elif theta is None and phi is not None:
                theta = np.argmax(self.data[phi::len(self.phis), :])
            elif theta is not None and phi is None:
                phi = np.argmax(self.data[theta * len(self.phis):(theta + 1) * len(self.phis):, :])

            # return arrays prepared for the plot with transition_radiation_visualizer.py
            return self.omegas, self.data[theta * len(self.phis) + phi, :]

    def get_iterations(self, species):
        """
        Return an array of iterations with available data.

        Parameters
        ----------
        species : string
            short name of the particle species, e.g. "e" for electrons
            (defined in ``speciesDefinition.param``)

        Returns
        -------
        An array with unsigned integers.
        """
        ctr_path = self.get_data_path(species)

        ctr_files = [f for f in os.listdir(ctr_path) if f.endswith(".dat")]
        iters = [int(f.split("_")[2].split(".")[0]) for f in ctr_files]
        return np.array(sorted(iters))
