"""lant: ambient-noise surface-wave tomography for dense arrays, from continuous records to a 3-D Vs model.

Three steps, run one after the other and checked in between (docs/tutorial.md):
    lant cc     cross-correlation of every station pair and its stack
    lant disp   reference curve and two-station Rayleigh phase velocities
    lant inv    1-D starting model and 3-D direct inversion of the travel times (Fang et al. 2015)
and the resolution test  lant checker.
"""

__version__ = "3.0.0"
