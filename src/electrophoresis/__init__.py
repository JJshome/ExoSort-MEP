"""
DC Electrophoresis module for ExoSort-MEP.

This module provides functionality for sorting extracellular vesicles
using direct current (DC) electrophoresis based on their membrane potential.
"""

from .dc_sorting import EVSample, DCElectrophoresis

__all__ = ['EVSample', 'DCElectrophoresis']
