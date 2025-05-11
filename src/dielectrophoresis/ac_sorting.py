#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
AC Dielectrophoresis-based Sorting of Extracellular Vesicles

This module implements the alternating current (AC) dielectrophoresis method for 
sorting extracellular vesicles based on their membrane potential differences.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
from scipy.constants import epsilon_0  # Vacuum permittivity
import cmath


@dataclass
class EVProperties:
    """Class for storing extracellular vesicle dielectric properties."""
    name: str
    diameter: float  # nm
    membrane_thickness: float  # nm
    membrane_permittivity: float  # relative permittivity (unitless)
    membrane_conductivity: float  # S/m
    interior_permittivity: float  # relative permittivity (unitless)
    interior_conductivity: float  # S/m
    membrane_potential: float  # mV
    zeta_potential: float  # mV
    origin: str
    
    def __str__(self):
        return (f"EVProperties: {self.name}\n"
                f"  Diameter: {self.diameter:.2f} nm\n"
                f"  Membrane Thickness: {self.membrane_thickness:.2f} nm\n"
                f"  Membrane Potential: {self.membrane_potential:.2f} mV\n"
                f"  Zeta Potential: {self.zeta_potential:.2f} mV\n"
                f"  Origin: {self.origin}")
    
    @property
    def radius(self) -> float:
        """Returns the radius in meters."""
        return self.diameter * 1e-9 / 2
    
    @property
    def membrane_thickness_m(self) -> float:
        """Returns the membrane thickness in meters."""
        return self.membrane_thickness * 1e-9


class ACDielectrophoresis:
    """
    Class for AC dielectrophoresis-based sorting of extracellular vesicles.
    
    This class handles the setup, running, and analysis of AC dielectrophoresis
    for sorting extracellular vesicles based on membrane potential differences.
    """
    
    def __init__(self, 
                 medium_permittivity: float = 78.5,  # relative permittivity of water
                 medium_conductivity: float = 1e-3,  # S/m
                 frequencies: List[float] = None,  # Hz
                 voltage: float = 10.0,  # V
                 electrode_distance: float = 100e-6,  # m
                 flow_rate: float = 1.0  # µL/min
                ):
        """
        Initialize the AC dielectrophoresis setup.
        
        Args:
            medium_permittivity: Relative permittivity of the suspension medium
            medium_conductivity: Conductivity of the suspension medium (S/m)
            frequencies: List of frequencies to analyze (Hz)
            voltage: Applied voltage (V)
            electrode_distance: Distance between electrodes (m)
            flow_rate: Flow rate through the microfluidic channel (µL/min)
        """
        self.medium_permittivity = medium_permittivity
        self.medium_conductivity = medium_conductivity
        
        # Default frequencies if none provided (logspace from 10^4 to 10^9 Hz)
        if frequencies is None:
            self.frequencies = np.logspace(4, 9, 100)
        else:
            self.frequencies = np.array(frequencies)
            
        self.voltage = voltage
        self.electrode_distance = electrode_distance
        self.flow_rate = flow_rate
        
        # Storage for vesicles and results
        self.vesicles = []
        self.results = None
        
    def add_vesicle(self, vesicle: EVProperties) -> None:
        """
        Add an EV to be analyzed.
        
        Args:
            vesicle: EVProperties object containing vesicle properties
        """
        self.vesicles.append(vesicle)
        print(f"Added vesicle: {vesicle.name}")
        
    def _complex_permittivity(self, permittivity: float, conductivity: float, 
                             frequency: float) -> complex:
        """
        Calculate complex permittivity.
        
        Args:
            permittivity: Relative permittivity
            conductivity: Conductivity (S/m)
            frequency: Frequency (Hz)
            
        Returns:
            Complex permittivity
        """
        omega = 2 * np.pi * frequency
        return permittivity * epsilon_0 - 1j * conductivity / omega
    
    def _complex_permittivity_medium(self, frequency: float) -> complex:
        """
        Calculate complex permittivity of the medium.
        
        Args:
            frequency: Frequency (Hz)
            
        Returns:
            Complex permittivity of medium
        """
        return self._complex_permittivity(
            self.medium_permittivity, 
            self.medium_conductivity, 
            frequency
        )
    
    def _complex_permittivity_vesicle(self, vesicle: EVProperties, 
                                     frequency: float) -> complex:
        """
        Calculate complex permittivity of a vesicle using single-shell model.
        
        Args:
            vesicle: EVProperties object
            frequency: Frequency (Hz)
            
        Returns:
            Effective complex permittivity of the vesicle
        """
        # Complex permittivity of membrane
        e_membrane = self._complex_permittivity(
            vesicle.membrane_permittivity,
            vesicle.membrane_conductivity,
            frequency
        )
        
        # Complex permittivity of interior
        e_interior = self._complex_permittivity(
            vesicle.interior_permittivity,
            vesicle.interior_conductivity,
            frequency
        )
        
        # Factor for membrane potential contribution
        # This is a simplified model where membrane potential affects conductivity
        # Actual mechanism is more complex and would require detailed modeling
        potential_factor = 1.0 + abs(vesicle.membrane_potential) / 100.0
        if vesicle.membrane_potential < 0:
            # Hyperpolarized membrane has increased effective conductivity
            e_membrane = e_membrane * potential_factor
        else:
            # Depolarized membrane has decreased effective conductivity
            e_membrane = e_membrane / potential_factor
            
        # Calculate radii
        r1 = vesicle.radius - vesicle.membrane_thickness_m  # Inner radius
        r2 = vesicle.radius  # Outer radius
        
        # Volume fraction of inner sphere
        gamma = (r1 / r2) ** 3
        
        # Calculate Clausius-Mossotti factor components
        e_2 = e_membrane
        e_1 = e_interior
        e_3 = self._complex_permittivity_medium(frequency)
        
        # First-order effective permittivity (inner sphere + membrane)
        numerator = 2 * e_2 + e_1 - 2 * gamma * (e_2 - e_1)
        denominator = 2 * e_2 + e_1 + gamma * (e_2 - e_1)
        e_effective = e_2 * (numerator / denominator)
        
        return e_effective
    
    def _clausius_mossotti_factor(self, vesicle: EVProperties, 
                                 frequency: float) -> complex:
        """
        Calculate the Clausius-Mossotti factor for a vesicle.
        
        Args:
            vesicle: EVProperties object
            frequency: Frequency (Hz)
            
        Returns:
            Complex Clausius-Mossotti factor
        """
        e_vesicle = self._complex_permittivity_vesicle(vesicle, frequency)
        e_medium = self._complex_permittivity_medium(frequency)
        
        return (e_vesicle - e_medium) / (e_vesicle + 2 * e_medium)
    
    def _dep_force(self, vesicle: EVProperties, frequency: float, 
                  gradient_e_squared: float = 1e12) -> float:
        """
        Calculate the DEP force on a vesicle.
        
        Args:
            vesicle: EVProperties object
            frequency: Frequency (Hz)
            gradient_e_squared: Gradient of the electric field squared (V²/m³)
            
        Returns:
            DEP force (N)
        """
        # Volume of vesicle
        volume = 4/3 * np.pi * vesicle.radius**3
        
        # Clausius-Mossotti factor
        cm_factor = self._clausius_mossotti_factor(vesicle, frequency)
        
        # Medium permittivity
        e_medium = self.medium_permittivity * epsilon_0
        
        # DEP force
        force = 2 * np.pi * e_medium * volume * cm_factor.real * gradient_e_squared
        
        return force
    
    def run_simulation(self) -> Dict[str, pd.DataFrame]:
        """
        Simulate the DEP behavior of added vesicles across frequencies.
        
        Returns:
            Dictionary mapping vesicle names to DataFrames with frequency response
        """
        if not self.vesicles:
            raise ValueError("No vesicles added. Use add_vesicle() before running.")
            
        # Gradient of electric field squared (approx for parallel electrodes)
        gradient_e_squared = (self.voltage / self.electrode_distance)**2 / self.electrode_distance
        
        # Dictionary to store results for each vesicle
        results = {}
        
        for vesicle in self.vesicles:
            # Calculate DEP force and CM factor for each frequency
            cm_factors_real = []
            cm_factors_imag = []
            dep_forces = []
            
            for freq in self.frequencies:
                cm_factor = self._clausius_mossotti_factor(vesicle, freq)
                cm_factors_real.append(cm_factor.real)
                cm_factors_imag.append(cm_factor.imag)
                
                dep_force = self._dep_force(vesicle, freq, gradient_e_squared)
                dep_forces.append(dep_force)
                
            # Create dataframe to store results
            df = pd.DataFrame({
                'Frequency_Hz': self.frequencies,
                'CM_Factor_Real': cm_factors_real,
                'CM_Factor_Imag': cm_factors_imag,
                'DEP_Force_N': dep_forces
            })
            
            # Add vesicle properties to dataframe for reference
            df['Vesicle_Name'] = vesicle.name
            df['Membrane_Potential_mV'] = vesicle.membrane_potential
            df['Diameter_nm'] = vesicle.diameter
            
            results[vesicle.name] = df
            
        self.results = results
        return results
    
    def find_crossover_frequencies(self) -> Dict[str, List[float]]:
        """
        Find crossover frequencies where the real part of CM factor changes sign.
        
        Returns:
            Dictionary mapping vesicle names to lists of crossover frequencies
        """
        if self.results is None:
            self.run_simulation()
            
        crossovers = {}
        
        for vesicle_name, df in self.results.items():
            # Find where real part of CM factor changes sign
            signs = np.sign(df['CM_Factor_Real'])
            sign_changes = np.where(np.diff(signs) != 0)[0]
            
            crossover_freqs = []
            for idx in sign_changes:
                # Linear interpolation to find more precise crossover point
                f1 = df['Frequency_Hz'].iloc[idx]
                f2 = df['Frequency_Hz'].iloc[idx + 1]
                cm1 = df['CM_Factor_Real'].iloc[idx]
                cm2 = df['CM_Factor_Real'].iloc[idx + 1]
                
                # Interpolate to find where CM factor = 0
                crossover = f1 - cm1 * (f2 - f1) / (cm2 - cm1)
                crossover_freqs.append(crossover)
                
            crossovers[vesicle_name] = crossover_freqs
            
        return crossovers
    
    def get_sorting_frequencies(self) -> List[float]:
        """
        Identify optimal frequencies for sorting vesicles based on CM factor differences.
        
        Returns:
            List of recommended frequencies for sorting
        """
        if self.results is None:
            self.run_simulation()
            
        # Get crossover frequencies for all vesicles
        crossovers = self.find_crossover_frequencies()
        
        # Combine all crossover frequencies
        all_crossovers = []
        for freqs in crossovers.values():
            all_crossovers.extend(freqs)
            
        # If there are no crossover frequencies, find frequencies with maximum CM factor difference
        if not all_crossovers:
            # Sample frequencies to test
            test_freqs = np.logspace(4, 9, 20)
            max_diff = 0
            best_freq = test_freqs[0]
            
            for freq in test_freqs:
                # Find closest frequency in our simulation results
                idx = np.abs(self.frequencies - freq).argmin()
                freq_actual = self.frequencies[idx]
                
                # Get CM factors for all vesicles at this frequency
                cm_factors = [df['CM_Factor_Real'].iloc[idx] for df in self.results.values()]
                
                # Calculate maximum difference between any two vesicles
                cm_diff = max(cm_factors) - min(cm_factors)
                
                if cm_diff > max_diff:
                    max_diff = cm_diff
                    best_freq = freq_actual
                    
            return [best_freq]
            
        # Return sorted unique crossover frequencies
        return sorted(list(set(all_crossovers)))
    
    def predict_sorting(self, frequency: float) -> Dict[str, float]:
        """
        Predict sorting outcomes at a specific frequency.
        
        Args:
            frequency: Frequency to analyze (Hz)
            
        Returns:
            Dictionary mapping vesicle names to relative DEP forces
        """
        if self.results is None:
            self.run_simulation()
            
        # Find closest frequency in our simulation results
        idx = np.abs(self.frequencies - frequency).argmin()
        freq_actual = self.frequencies[idx]
        
        # Get DEP forces for all vesicles at this frequency
        sorting = {}
        for vesicle_name, df in self.results.items():
            force = df['DEP_Force_N'].iloc[idx]
            cm_factor = df['CM_Factor_Real'].iloc[idx]
            
            sorting[vesicle_name] = {
                'DEP_Force_N': force,
                'CM_Factor': cm_factor,
                'Direction': 'Positive DEP (to high field)' if cm_factor > 0 else 'Negative DEP (to low field)'
            }
            
        return sorting
    
    def visualize_results(self, filename: Optional[str] = None, 
                         highlight_frequencies: List[float] = None) -> None:
        """
        Generate visualization of the DEP results.
        
        Args:
            filename: Optional filename to save the plot
            highlight_frequencies: List of frequencies to highlight on the plot
        """
        if self.results is None:
            self.run_simulation()
            
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
        
        # Color map for vesicles based on membrane potential
        potentials = [v.membrane_potential for v in self.vesicles]
        norm = plt.Normalize(min(potentials), max(potentials))
        cmap = plt.cm.viridis
        
        # Plot real part of CM factor vs frequency
        for i, (vesicle_name, df) in enumerate(self.results.items()):
            vesicle = next(v for v in self.vesicles if v.name == vesicle_name)
            color = cmap(norm(vesicle.membrane_potential))
            
            ax1.semilogx(df['Frequency_Hz'], df['CM_Factor_Real'], 
                       label=f"{vesicle_name} ({vesicle.membrane_potential} mV)", 
                       color=color, linewidth=2)
            
        # Add zero line and highlight frequencies
        ax1.axhline(y=0, color='r', linestyle='--', alpha=0.7)
        
        if highlight_frequencies:
            for freq in highlight_frequencies:
                ax1.axvline(x=freq, color='k', linestyle=':', alpha=0.7)
                ax1.text(freq, ax1.get_ylim()[1] * 0.9, f"{freq:.2e} Hz", 
                        rotation=90, va='top')
                
        ax1.set_ylabel('Real Part of CM Factor')
        ax1.set_title('Clausius-Mossotti Factor vs Frequency')
        ax1.grid(True, which='both', linestyle='--', alpha=0.5)
        ax1.legend(loc='best')
        
        # Plot DEP force vs frequency
        for i, (vesicle_name, df) in enumerate(self.results.items()):
            vesicle = next(v for v in self.vesicles if v.name == vesicle_name)
            color = cmap(norm(vesicle.membrane_potential))
            
            ax2.loglog(df['Frequency_Hz'], np.abs(df['DEP_Force_N']), 
                     label=f"{vesicle_name} ({vesicle.membrane_potential} mV)", 
                     color=color, linewidth=2)
            
        # Highlight frequencies
        if highlight_frequencies:
            for freq in highlight_frequencies:
                ax2.axvline(x=freq, color='k', linestyle=':', alpha=0.7)
                
        ax2.set_xlabel('Frequency (Hz)')
        ax2.set_ylabel('|DEP Force| (N)')
        ax2.set_title('DEP Force Magnitude vs Frequency')
        ax2.grid(True, which='both', linestyle='--', alpha=0.5)
        
        # Add colorbar for membrane potential
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=[ax1, ax2], orientation='vertical', pad=0.01)
        cbar.set_label('Membrane Potential (mV)')
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()
        
    def visualize_3d_sorting(self, frequencies: List[float], 
                            show_vectors: bool = True,
                            filename: Optional[str] = None) -> None:
        """
        Generate 3D visualization of DEP-based sorting at given frequencies.
        
        Args:
            frequencies: List of frequencies to visualize
            show_vectors: Whether to show force vectors
            filename: Optional filename to save the plot
        """
        if self.results is None:
            self.run_simulation()
            
        if len(frequencies) != 2 and len(frequencies) != 3:
            raise ValueError("Provide either 2 or 3 frequencies for visualization")
            
        # Find indices of closest frequencies in our results
        indices = []
        actual_freqs = []
        for freq in frequencies:
            idx = np.abs(self.frequencies - freq).argmin()
            indices.append(idx)
            actual_freqs.append(self.frequencies[idx])
            
        # Get CM factors for all vesicles at selected frequencies
        vesicle_colors = []
        cm_values = []
        names = []
        
        for vesicle_name, df in self.results.items():
            vesicle = next(v for v in self.vesicles if v.name == vesicle_name)
            
            # Get CM factors at selected frequencies
            cm_vals = [df['CM_Factor_Real'].iloc[idx] for idx in indices]
            
            # Add 0 for z-axis if only 2 frequencies provided
            if len(frequencies) == 2:
                cm_vals.append(0)
                
            cm_values.append(cm_vals)
            vesicle_colors.append(vesicle.membrane_potential)
            names.append(vesicle_name)
            
        # Create 3D plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Convert to numpy array for easier manipulation
        cm_values = np.array(cm_values)
        
        # Color map based on membrane potential
        norm = plt.Normalize(min(vesicle_colors), max(vesicle_colors))
        cmap = plt.cm.viridis
        colors = [cmap(norm(p)) for p in vesicle_colors]
        
        # Plot points
        if len(frequencies) == 2:
            # 2D plot in 3D space (z=0)
            ax.scatter(cm_values[:, 0], cm_values[:, 1], cm_values[:, 2], 
                     c=colors, s=100, marker='o', edgecolors='k')
            
            # Show vectors from origin if requested
            if show_vectors:
                for i, (x, y, z) in enumerate(cm_values):
                    ax.plot([0, x], [0, y], [0, z], color=colors[i], alpha=0.7)
                    
            ax.set_zlim(-1, 1)
        else:
            # 3D plot
            ax.scatter(cm_values[:, 0], cm_values[:, 1], cm_values[:, 2], 
                     c=colors, s=100, marker='o', edgecolors='k')
            
            # Show vectors from origin if requested
            if show_vectors:
                for i, (x, y, z) in enumerate(cm_values):
                    ax.plot([0, x], [0, y], [0, z], color=colors[i], alpha=0.7)
        
        # Add labels for points
        for i, name in enumerate(names):
            ax.text(cm_values[i, 0], cm_values[i, 1], cm_values[i, 2], 
                  name, size=10, zorder=1)
            
        # Set labels and title
        ax.set_xlabel(f'CM Factor at {actual_freqs[0]:.2e} Hz')
        ax.set_ylabel(f'CM Factor at {actual_freqs[1]:.2e} Hz')
        
        if len(frequencies) == 3:
            ax.set_zlabel(f'CM Factor at {actual_freqs[2]:.2e} Hz')
            title = '3D Sorting Space Based on CM Factors at Three Frequencies'
        else:
            ax.set_zlabel('No third frequency')
            title = '2D Sorting Space Based on CM Factors at Two Frequencies'
            
        ax.set_title(title)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.1)
        cbar.set_label('Membrane Potential (mV)')
        
        # Add grid
        ax.grid(True)
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()


# Example usage:
if __name__ == "__main__":
    # Create sample EV vesicles with different membrane potentials
    vesicles = [
        EVProperties(
            name="Hyperpolarized EVs",
            diameter=100.0,  # nm
            membrane_thickness=5.0,  # nm
            membrane_permittivity=10.0,
            membrane_conductivity=1e-7,  # S/m
            interior_permittivity=50.0,
            interior_conductivity=0.5,  # S/m
            membrane_potential=-90.0,  # mV (hyperpolarized)
            zeta_potential=-30.0,  # mV
            origin="Neural cells"
        ),
        EVProperties(
            name="Normal EVs",
            diameter=120.0,  # nm
            membrane_thickness=5.0,  # nm
            membrane_permittivity=8.0,
            membrane_conductivity=1e-7,  # S/m
            interior_permittivity=60.0,
            interior_conductivity=0.6,  # S/m
            membrane_potential=-70.0,  # mV (normal)
            zeta_potential=-25.0,  # mV
            origin="Fibroblasts"
        ),
        EVProperties(
            name="Depolarized EVs",
            diameter=130.0,  # nm
            membrane_thickness=4.5,  # nm
            membrane_permittivity=7.0,
            membrane_conductivity=2e-7,  # S/m
            interior_permittivity=55.0,
            interior_conductivity=0.7,  # S/m
            membrane_potential=-30.0,  # mV (depolarized)
            zeta_potential=-15.0,  # mV
            origin="Cancer cells"
        ),
        EVProperties(
            name="Positive EVs",
            diameter=90.0,  # nm
            membrane_thickness=4.0,  # nm
            membrane_permittivity=9.0,
            membrane_conductivity=3e-7,  # S/m
            interior_permittivity=45.0,
            interior_conductivity=0.4,  # S/m
            membrane_potential=10.0,  # mV (positive, rare)
            zeta_potential=5.0,  # mV
            origin="Activated immune cells"
        )
    ]
    
    # Setup DEP
    dep = ACDielectrophoresis(
        medium_permittivity=78.5,
        medium_conductivity=0.05,  # Higher conductivity buffer
        frequencies=np.logspace(3, 9, 200),  # 1 kHz to 1 GHz
        voltage=10.0,
        electrode_distance=50e-6,  # 50 µm
        flow_rate=0.5  # µL/min
    )
    
    # Add vesicles
    for vesicle in vesicles:
        dep.add_vesicle(vesicle)
    
    # Run simulation
    results = dep.run_simulation()
    
    # Get sorting frequencies
    sorting_freqs = dep.get_sorting_frequencies()
    print("\nRecommended sorting frequencies:")
    for freq in sorting_freqs:
        print(f"  {freq:.2e} Hz")
    
    # Predict sorting at first recommended frequency
    if sorting_freqs:
        freq = sorting_freqs[0]
        print(f"\nPredicted sorting at {freq:.2e} Hz:")
        sorting = dep.predict_sorting(freq)
        
        for vesicle_name, info in sorting.items():
            print(f"  {vesicle_name}:")
            print(f"    CM Factor: {info['CM_Factor']:.4f}")
            print(f"    DEP Force: {info['DEP_Force_N']:.4e} N")
            print(f"    Direction: {info['Direction']}")
    
    # Visualize results
    dep.visualize_results("dep_frequency_response.png", highlight_frequencies=sorting_freqs)
    
    # 3D visualization with 2 frequencies
    if len(sorting_freqs) >= 2:
        dep.visualize_3d_sorting(sorting_freqs[:2], show_vectors=True, 
                               filename="dep_2d_sorting.png")
    
    # 3D visualization with 3 frequencies if available
    if len(sorting_freqs) >= 3:
        dep.visualize_3d_sorting(sorting_freqs[:3], show_vectors=True, 
                               filename="dep_3d_sorting.png")
