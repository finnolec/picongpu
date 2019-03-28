"""
This file is part of the PIConGPU.

Copyright 2017-2019 PIConGPU contributors
Authors: Sebastian Starke
License: GPLv3+
"""

import matplotlib.pyplot as plt
from warnings import warn

from picongpu.plugins.plot_mpl.utils import get_different_colors


class Visualizer(object):
    """
    Abstract base class for matplotlib visualizers that implements
    the visualization logic.
    Classes that derive from this class need to write their own implementations
    for the following functions in order to work:
        _create_data_reader(self, run_directory)
        _create_plt_obj(self, ax)
        _update_plt_obj(self)

    Note: When using classes derived from this within jupyter notebooks, use
    %matplotlib notebook mode.
    """

    def __init__(self, reader_cls, run_directories=None, ax=None):
        """
        Initialize the reader and data as member parameters.

        Parameters
        ----------
        run_directories : list of tuples of length 2
            or single tuple of length 2
            or list of strings or string.
            If tuples are specified, they have to be
            of the following form (sim_label, sim_path) and consist
            of strings with 'sim_label' being a short string
            used e.g. in plot legends and 'sim_path'
            leading to the run directory of PIConGPU
            (the path before ``simOutput/``).
            If only a string or list of strings is given,
            labels are generated by enumeration.
            If None, the user is responsible for providing run_directories
            later on via set_run_directories() before calling visualize().
        reader_cls: class
            handle of a PIConGPU Data reader class (not string!)
            which inherited from BaseReader in plugins.data.base_reader.
        ax: matplotlib.axes
        """
        self.reader_cls = reader_cls

        if ax is None:
            warn("No axes was given, using plt.gca() instead!")
            ax = plt.gca()
        self.ax = ax

        # they will only be set when run_directories are provided
        self.data_reader = None
        self.data = None
        self.plt_obj = None
        self.sim_labels = None
        self.colors = None

        if run_directories is not None:
            self.set_run_directories(run_directories)

    def _init_members(self, run_directories):
        # 0) clean the ax from previous data (own or other objects)
        self._clean_ax()

        # 1) get new labels
        self.sim_labels = [str(run_dir[0]) for run_dir in run_directories]

        # 2) reinstantiate the readers
        self._init_reader(run_directories)

        # 3) set data to None
        self._init_data(run_directories)

        # 4) set all plt_objects to None
        self._init_plt_obj(run_directories)

        # 5) set all colors or colormaps
        self._init_colors(run_directories)

    def set_run_directories(self, run_directories):
        """
        Prepare everything for a fresh start with new run_directories
        """
        run_directories = self._check_and_fix_run_dirs(run_directories)

        self._init_members(run_directories)

    def _clean_ax(self):
        """
        clear the ax and if available also clear possible colorbars
        that are still left from previous plots
        """
        for idx in range(len(self.ax.images)):
            if self.ax.images[idx].colorbar is not None:
                # this also removes the ax that the colorbar
                # was plotted in from the figure!
                self.ax.images[idx].colorbar.remove()
        self.ax.clear()

    def _init_reader(self, run_directories):
        self.data_reader = [
            self.reader_cls(run_dir[1]) for run_dir in run_directories]

    def _init_data(self, run_directories):
        self.data = [None] * len(run_directories)

    def _init_plt_obj(self, run_directories):
        self.plt_obj = [None] * len(run_directories)

    def _init_colors(self, run_directories):
        self.colors = get_different_colors(len(run_directories))

    def _check_and_fix_run_dirs(self, run_directories):
        """
        Check variable type for the run_directories and change
        to list of tuples if necessary.
        This can be overridden in derived classes to e.g. restrict
        to single simulation visualization.

        Returns
        -------
        a list of tuples, each of the form
        (simulation_label, path_to_simulation).
        """
        # silently convert str to list of length 1
        if not isinstance(run_directories, list):
            run_directories = [run_directories]

        if len(run_directories) < 1:
            warn("Empty run_directories list was passed!")
            return run_directories

        if isinstance(run_directories[0], str):
            warn("First element is str. Assuming the same for all "
                 "other elements. Will use enumeration for labeling!")
            run_directories = list(enumerate(run_directories))

        return run_directories

    def _remove_plt_obj(self, idx):
        """
        Removes the plt_obj at position idx from the current plot
        and sets it to None so that in a subsequent visualization call
        the plt_obj is created again.
        """
        # clear the previous data from the plot
        self.plt_obj[idx].remove()
        self.plt_obj[idx] = None

    def _create_plt_obj(self, idx):
        """
        Sets 'self.plt_obj' to an instance of a matplotlib.artist.Artist
        object (or derived classes) created by using 'self.ax'
        which can later be updated by feeding new data into it.
        Only called on the first call for visualization.
        """
        raise NotImplementedError

    def _update_plt_obj(self, idx):
        """
        Take the 'self.data' member, interpret it and feed it into the
        'self.plt_obj'.
        """
        raise NotImplementedError

    def visualize(self, **kwargs):
        """
        1. gathers the data for the selected kwargs
        2. removes plot elements for sources which have no data
        3. plot the data
        3.a Creates the 'plt_obj' if it does not exist
        3.b Updates the 'plt_obj' with the new data.
        4. adjusts the plot
        """
        if self.data_reader is None:
            raise RuntimeError(
                "No data readers are given. Did you provide run_directories"
                " during construction?"
                " If not, call set_run_directories() before calling"
                " visualize()!")

        # check that not multiple iteratons or timesteps are selected.
        # Since time is preferred over iteration if both are passed,
        # we check time first.
        # Note: the readers might be capable of dealing with multiple
        # values but visualization is not.
        if 'time' in kwargs:
            if isinstance(kwargs['time'], list):
                raise ValueError(
                    "This class only supports single timestep visualization!")
        if 'iteration' in kwargs:
            if isinstance(kwargs['iteration'], list):
                raise ValueError(
                    "This class only supports single iteration visualization!")

        self.collect_data(**kwargs)

        # first remove the stuff for which we don't have data from the ax
        self.remove_plots_without_data()

        # second plot everything for which we have data
        self.draw_data()

        self.adjust_plot(**kwargs)

    def collect_data(self, **kwargs):
        # get the data for the given args from all the readers
        # and catch errors when iteration is not there
        for i, reader in enumerate(self.data_reader):
            try:
                d = reader.get(**kwargs)
                self.data[i] = d
            except IndexError:
                # iteration is not there, so the data value
                # is going to be None
                self.data[i] = None

    def remove_plots_without_data(self):
        for i, (plt_obj, data) in enumerate(zip(self.plt_obj, self.data)):
            if data is None:
                # no data for the current combination of kwargs
                # so we should just omit the data from this input file
                warn("Data from {} is None! Will omit!".format(
                    self.sim_labels[i]))
                if plt_obj is not None:
                    self._remove_plt_obj(i)

    def draw_data(self):
        for i, (plt_obj, data) in enumerate(zip(self.plt_obj, self.data)):
            if data is not None:
                # we have data and only decide about creating/updating
                if plt_obj is None:
                    # create a new plot object for data from file i
                    self._create_plt_obj(i)
                else:
                    # update the plot object i by updating its data
                    self._update_plt_obj(i)

    def adjust_plot(self, **kwargs):
        """
        Executed after the plotting is done for adjusting legends etc...
        """
        pass

    def clear_cbar(self):
        """
        Clear colorbars if present. Should be implemented
        in derived classes that use colorbars.
        """
        pass
