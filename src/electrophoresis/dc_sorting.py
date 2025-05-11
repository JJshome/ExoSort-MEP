#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
DC Electrophoresis-based Sorting of Extracellular Vesicles

This module implements the direct current (DC) electrophoresis method for 
sorting extracellular vesicles based on their membrane potential differences.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass


@dataclass
class EVSample:
    """Class for storing extracellular vesicle sample information."""
    name: str
    concentration: float  # particles/mL
    size_distribution: Dict[str, float]  # nm: percentage
    zeta_potential: float  # mV
    membrane_potential: float  # mV
    origin: str
    
    def __str__(self):
        return (f"EVSample: {self.name}\n"
                f"  Concentration: {self.concentration:.2e} particles/mL\n"
                f"  Avg Size: {self._avg_size():.2f} nm\n"
                f"  Zeta Potential: {self.zeta_potential:.2f} mV\n"
                f"  Membrane Potential: {self.membrane_potential:.2f} mV\n"
                f"  Origin: {self.origin}")
    
    def _avg_size(self) -> float:
        """Calculate average size from distribution."""
        sizes = np.array(list(self.size_distribution.keys()), dtype=float)
        percentages = np.array(list(self.size_distribution.values()))
        return np.sum(sizes * percentages) / np.sum(percentages)


class DCElectrophoresis:
    """
    Class for DC electrophoresis-based sorting of extracellular vesicles.
    
    This class handles the setup, running, and analysis of DC electrophoresis
    for sorting extracellular vesicles based on membrane potential differences.
    """
    
    def __init__(self, 
                 buffer_conductivity: float = 1.5,  # mS/cm
                 buffer_ph: float = 7.4,
                 voltage: float = 100.0,  # V
                 gel_percentage: float = 0.8,  # %
                 run_time: float = 60.0,  # minutes
                 collection_intervals: List[float] = None  # mm
                ):
        """
        Initialize the DC electrophoresis setup.
        
        Args:
            buffer_conductivity: Conductivity of the buffer (mS/cm)
            buffer_ph: pH of the buffer
            voltage: Applied voltage (V)
            gel_percentage: Agarose gel percentage (%)
            run_time: Total run time (minutes)
            collection_intervals: Distance intervals for collecting vesicles (mm)
        """
        self.buffer_conductivity = buffer_conductivity
        self.buffer_ph = buffer_ph
        self.voltage = voltage
        self.gel_percentage = gel_percentage
        self.run_time = run_time
        
        # Default collection intervals if none provided
        if collection_intervals is None:
            self.collection_intervals = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
        else:
            self.collection_intervals = collection_intervals
            
        # Storage for results
        self.samples = []
        self.results = None
        
    def add_sample(self, sample: EVSample) -> None:
        """
        Add an EV sample to be sorted.
        
        Args:
            sample: EVSample object containing sample information
        """
        self.samples.append(sample)
        print(f"Added sample: {sample.name}")
        
    def run_simulation(self) -> Dict[str, pd.DataFrame]:
        """
        Simulate the electrophoresis run with the added samples.
        
        Returns:
            Dictionary mapping sample names to DataFrames with sorting results
        """
        if not self.samples:
            raise ValueError("No samples added. Use add_sample() before running.")
            
        # Dictionary to store results for each sample
        results = {}
        
        for sample in self.samples:
            # Calculate mobility based on membrane potential
            # This is a simplified model; actual mobility would depend on many factors
            mobility = -0.2 * sample.membrane_potential  # µm·cm/V·s
            
            # Field strength
            field_strength = self.voltage / 100.0  # V/cm (assumes 10 cm gel)
            
            # Calculate velocity
            velocity = mobility * field_strength  # µm/s
            
            # Convert to mm/min for easier interpretation
            velocity_mm_min = velocity * 60 / 1000
            
            # Calculate distance traveled
            distance = velocity_mm_min * self.run_time  # mm
            
            # Create dataframe to store results
            df = pd.DataFrame({
                'Sample': [sample.name],
                'Membrane_Potential_mV': [sample.membrane_potential],
                'Mobility_µm_cm_per_Vs': [mobility],
                'Velocity_mm_per_min': [velocity_mm_min],
                'Distance_mm': [distance]
            })
            
            results[sample.name] = df
            
        self.results = results
        return results
    
    def collect_fractions(self) -> Dict[str, Dict[str, float]]:
        """
        Collect sorted vesicle fractions based on the simulation results.
        
        Returns:
            Dictionary mapping fractions to sample concentrations
        """
        if self.results is None:
            self.run_simulation()
            
        # Dictionary to store collected fractions
        # Format: {interval: {sample_name: concentration}}
        fractions = {}
        
        # Initialize intervals
        for interval in self.collection_intervals:
            fractions[interval] = {}
            
        # Assign samples to fractions based on distance traveled
        for sample_name, df in self.results.items():
            distance = df['Distance_mm'].values[0]
            
            # Find appropriate interval
            for i, interval in enumerate(sorted(self.collection_intervals)):
                if distance <= interval:
                    # Add sample to this fraction
                    sample = next(s for s in self.samples if s.name == sample_name)
                    fractions[interval][sample_name] = sample.concentration
                    break
                    
                # If it's the last interval and distance is greater, put in last fraction
                if i == len(self.collection_intervals) - 1:
                    sample = next(s for s in self.samples if s.name == sample_name)
                    fractions[interval][sample_name] = sample.concentration
                    
        return fractions
    
    def visualize_results(self, filename: Optional[str] = None) -> None:
        """
        Generate visualization of the electrophoresis results.
        
        Args:
            filename: Optional filename to save the plot
        """
        if self.results is None:
            self.run_simulation()
            
        # Set up figure
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))
        
        # Plot membrane potential vs. distance
        potentials = []
        distances = []
        sample_names = []
        
        for sample_name, df in self.results.items():
            potentials.append(df['Membrane_Potential_mV'].values[0])
            distances.append(df['Distance_mm'].values[0])
            sample_names.append(sample_name)
            
        # Plot distance vs. membrane potential
        scatter = ax1.scatter(potentials, distances, c=potentials, cmap='viridis', 
                             s=100, alpha=0.7, edgecolors='k')
        
        # Add labels for points
        for i, name in enumerate(sample_names):
            ax1.annotate(name, (potentials[i], distances[i]), 
                       xytext=(5, 5), textcoords='offset points')
            
        ax1.set_xlabel('Membrane Potential (mV)')
        ax1.set_ylabel('Distance Traveled (mm)')
        ax1.set_title('Membrane Potential vs. Distance')
        ax1.grid(True, linestyle='--', alpha=0.7)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax1)
        cbar.set_label('Membrane Potential (mV)')
        
        # Plot gel visualization
        ax2.set_xlim(0, 10)  # gel width
        ax2.set_ylim(0, max(distances) * 1.2)  # gel length + margin
        
        # Draw gel
        gel_rect = plt.Rectangle((0, 0), 10, max(distances) * 1.2, 
                               fc='lightblue', ec='skyblue', alpha=0.3)
        ax2.add_patch(gel_rect)
        
        # Draw wells
        for i in range(len(sample_names)):
            well = plt.Rectangle((3 + i, 0), 0.8, 2, fc='black', ec='black')
            ax2.add_patch(well)
            
        # Draw bands
        for i, (name, dist) in enumerate(zip(sample_names, distances)):
            band_color = plt.cm.viridis(plt.Normalize(-100, 0)(potentials[i]))
            band = plt.Rectangle((3 + i, dist - 1), 0.8, 2, fc=band_color, ec='black')
            ax2.add_patch(band)
            ax2.text(3 + i + 0.4, dist + 2, name, ha='center')
            
        # Draw collection intervals
        for interval in self.collection_intervals:
            ax2.axhline(y=interval, color='red', linestyle='--', alpha=0.7)
            ax2.text(0.5, interval + 0.5, f"{interval} mm", color='red')
            
        ax2.set_title('Gel Visualization')
        ax2.set_xticks([])
        ax2.set_yticks([])
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()


# Example usage:
if __name__ == "__main__":
    # Create sample EV samples with different membrane potentials
    samples = [
        EVSample(
            name="Hyperpolarized EVs",
            concentration=1e9,
            size_distribution={80: 0.1, 100: 0.3, 120: 0.4, 150: 0.2},
            zeta_potential=-30.5,
            membrane_potential=-90.0,  # Highly negative (hyperpolarized)
            origin="Neural cells"
        ),
        EVSample(
            name="Normal EVs",
            concentration=5e9,
            size_distribution={90: 0.2, 110: 0.5, 130: 0.3},
            zeta_potential=-20.2,
            membrane_potential=-70.0,  # Normal negative potential
            origin="Fibroblasts"
        ),
        EVSample(
            name="Depolarized EVs",
            concentration=3e9,
            size_distribution={100: 0.4, 120: 0.3, 140: 0.3},
            zeta_potential=-15.5,
            membrane_potential=-30.0,  # Less negative (depolarized)
            origin="Cancer cells"
        ),
        EVSample(
            name="Positive EVs",
            concentration=2e8,
            size_distribution={70: 0.2, 90: 0.5, 110: 0.3},
            zeta_potential=5.5,
            membrane_potential=10.0,  # Positive potential (rare)
            origin="Activated immune cells"
        )
    ]
    
    # Setup electrophoresis
    ep = DCElectrophoresis(
        buffer_conductivity=1.5,
        buffer_ph=7.4,
        voltage=120.0,
        gel_percentage=0.8,
        run_time=90.0,
        collection_intervals=[5.0, 10.0, 15.0, 20.0, 25.0]
    )
    
    # Add samples
    for sample in samples:
        ep.add_sample(sample)
    
    # Run simulation
    results = ep.run_simulation()
    
    # Print results
    for sample_name, df in results.items():
        print(f"\nResults for {sample_name}:")
        print(df.to_string(index=False))
    
    # Collect fractions
    fractions = ep.collect_fractions()
    print("\nCollected fractions:")
    for interval, samples in fractions.items():
        print(f"Fraction at {interval} mm:")
        for sample_name, conc in samples.items():
            print(f"  {sample_name}: {conc:.2e} particles/mL")
    
    # Visualize results
    ep.visualize_results("electrophoresis_results.png")
