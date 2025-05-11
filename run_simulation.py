#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ExoSort-MEP Simulation Launcher

This script launches the web-based simulation environment for exploring the
electrophoresis and dielectrophoresis sorting of extracellular vesicles.
"""

import sys
import panel as pn
from deployment.simulation.ev_sorting_simulator import EVSortingSimulator


def main():
    """Launch the simulation app."""
    print("Starting ExoSort-MEP Simulation...")
    
    # Initialize Panel extensions
    pn.extension('bokeh')
    
    # Create simulator
    simulator = EVSortingSimulator()
    app = simulator.get_app()
    
    # Serve the app
    print("\nSimulation is running! Open a browser and navigate to:")
    print("    http://localhost:5006")
    print("\nPress Ctrl+C to exit.")
    try:
        pn.serve(app, port=5006, show=True)
    except KeyboardInterrupt:
        print("\nExiting simulation.")
        sys.exit(0)


if __name__ == "__main__":
    main()
