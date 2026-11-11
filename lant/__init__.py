"""local_ANT: high-frequency ambient-noise surface-wave tomography for dense (nodal) arrays.

Raw records of any layout and format -> cleaned archive (parameters from the data) -> hourly cross-coherence ->
stack and station diagnosis -> F-J reference curve -> two-station phase velocities with quality
control -> DSurfTomo direct inversion for 3-D Vs (areal arrays) or narrow-strip inversion along a line -> checks and
figures.  Command line: `lant --help`; documentation: docs/tutorial.md.
"""

__version__ = "1.4.0"
