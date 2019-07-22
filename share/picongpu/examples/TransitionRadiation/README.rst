TransitionRadiation : Transtion Radiation
=============================================

.. sectionauthor:: Finn-Ole Carstens <f.carstens (at) hzdr.de>

This example simulates an electron bunch moving in a 45° angle in the x-y plane and calculates the coherent and incoherent transition radiation created by this bunch.
The implemented transition radiation follows the studies from [Schroeder2004] and [Downer2018].
The transition radiation is computed for an infinitely large interface perpendicular to the y-axis of the simulation.
The result can be interpreted with the according python script `lib/python/picongpu/plugins/plot_mpl/transition_radiation_visualizer.py`.

References
----------

.. [Schroeder2004]
       Schroeder, C. B. and Esarey, E. and van Tilborg, J. and Leemans, W. P.
       *Theory of coherent transition radiation generated at a plasma-vacuum interface*,
       American Physical Society(2004),
       https://link.aps.org/doi/10.1103/PhysRevE.69.016501

.. [Downer2018]
       Downer, M. C. and Zgadzaj, R. and Debus, A. and Schramm, U. and Kaluza, M. C.
       *Diagnostics for plasma-based electron accelerators*,
       American Physical Society(2018),
       https://link.aps.org/doi/10.1103/RevModPhys.90.035002
