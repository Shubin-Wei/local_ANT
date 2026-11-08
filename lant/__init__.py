"""local_ANT: high-frequency ambient-noise surface-wave tomography for dense (nodal) arrays.

Raw records of any layout and format -> cleaned archive -> hourly cross-coherence -> station-day / station-hour
screening -> notch -> F-J and median reference curves -> two-station phase and group velocities with quality
control -> DSurfTomo direct inversion for 3-D Vs (areal arrays) or narrow-strip inversion along a line -> checks and
figures.  Command line: `lant --help`; documentation: docs/tutorial.md.
"""
__version__ = "1.0.0"
