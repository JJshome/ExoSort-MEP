#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extracellular Vesicle Data Handling Utilities

This module provides utility functions for loading, saving, and converting 
extracellular vesicle data across various formats used in the ExoSort-MEP system.
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Union, Optional, Tuple, Any


def load_ev_data(file_path: str) -> pd.DataFrame:
    """
    Load extracellular vesicle data from various file formats.
    
    Args:
        file_path: Path to the data file (CSV, Excel, JSON, TXT)
        
    Returns:
        DataFrame containing the EV data
    
    Raises:
        ValueError: If file format is not supported
    """
    file_ext = os.path.splitext(file_path)[1].lower()
    
    if file_ext == '.csv':
        return pd.read_csv(file_path)
    elif file_ext in ['.xls', '.xlsx']:
        return pd.read_excel(file_path)
    elif file_ext == '.json':
        with open(file_path, 'r') as f:
            data = json.load(f)
        return pd.DataFrame(data)
    elif file_ext == '.txt':
        # Try to determine delimiter automatically
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
        
        if '\t' in first_line:
            return pd.read_csv(file_path, sep='\t')
        elif ',' in first_line:
            return pd.read_csv(file_path, sep=',')
        elif ';' in first_line:
            return pd.read_csv(file_path, sep=';')
        else:
            return pd.read_csv(file_path, delim_whitespace=True)
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")


def save_ev_data(data: pd.DataFrame, file_path: str, **kwargs) -> None:
    """
    Save extracellular vesicle data to various file formats.
    
    Args:
        data: DataFrame containing EV data
        file_path: Path where the file will be saved
        **kwargs: Additional parameters for the specific save function
    
    Raises:
        ValueError: If file format is not supported
    """
    file_ext = os.path.splitext(file_path)[1].lower()
    
    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)
    
    if file_ext == '.csv':
        data.to_csv(file_path, index=False, **kwargs)
    elif file_ext in ['.xls', '.xlsx']:
        data.to_excel(file_path, index=False, **kwargs)
    elif file_ext == '.json':
        with open(file_path, 'w') as f:
            json.dump(data.to_dict('records'), f, indent=2)
    elif file_ext == '.txt':
        data.to_csv(file_path, sep='\t', index=False, **kwargs)
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")


def convert_nta_to_standard(
    nta_file: str, 
    output_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Convert NTA (Nanoparticle Tracking Analysis) data to standard format.
    
    Args:
        nta_file: Path to the NTA data file
        output_file: Optional path to save the converted data
        
    Returns:
        DataFrame in standard format
    """
    # Load NTA data based on file extension
    file_ext = os.path.splitext(nta_file)[1].lower()
    
    if file_ext == '.csv':
        df = pd.read_csv(nta_file)
    elif file_ext in ['.xls', '.xlsx']:
        df = pd.read_excel(nta_file)
    else:
        raise ValueError(f"Unsupported NTA file format: {file_ext}")
    
    # Check if file is from Malvern NanoSight
    if 'Size (nm)' in df.columns and 'Concentration (particles/ml)' in df.columns:
        # Already in standard format
        std_df = df
    
    # Check if file is from ZetaView
    elif 'Diameter (nm)' in df.columns and 'Particles' in df.columns:
        std_df = pd.DataFrame({
            'Size (nm)': df['Diameter (nm)'],
            'Concentration (particles/ml)': df['Particles'],
        })
        # Add other columns if available
        if 'Mean Intensity' in df.columns:
            std_df['Intensity'] = df['Mean Intensity']
            
    # Check if file is from IZON qNano
    elif 'Particle Diameter (nm)' in df.columns and 'Concentration (particles/mL)' in df.columns:
        std_df = pd.DataFrame({
            'Size (nm)': df['Particle Diameter (nm)'],
            'Concentration (particles/ml)': df['Concentration (particles/mL)'],
        })
        
    # Unknown format - try to guess based on column names
    else:
        # Look for size column
        size_cols = [col for col in df.columns if any(term in col.lower() 
                                                     for term in ['size', 'diameter', 'nm'])]
        # Look for concentration column
        conc_cols = [col for col in df.columns if any(term in col.lower() 
                                                     for term in ['concentration', 'particles', 'count'])]
        
        if size_cols and conc_cols:
            std_df = pd.DataFrame({
                'Size (nm)': df[size_cols[0]],
                'Concentration (particles/ml)': df[conc_cols[0]],
            })
        else:
            raise ValueError("Could not determine the format of the NTA data file")
    
    # Save to output file if specified
    if output_file:
        save_ev_data(std_df, output_file)
        
    return std_df


def convert_zeta_to_standard(
    zeta_file: str, 
    output_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Convert zeta potential data to standard format.
    
    Args:
        zeta_file: Path to the zeta potential data file
        output_file: Optional path to save the converted data
        
    Returns:
        DataFrame in standard format
    """
    # Load zeta data based on file extension
    file_ext = os.path.splitext(zeta_file)[1].lower()
    
    if file_ext == '.csv':
        df = pd.read_csv(zeta_file)
    elif file_ext in ['.xls', '.xlsx']:
        df = pd.read_excel(zeta_file)
    else:
        raise ValueError(f"Unsupported zeta file format: {file_ext}")
    
    # Check if file is from Malvern Zetasizer
    if 'Zeta Potential (mV)' in df.columns:
        # Already in standard format
        std_df = df
    
    # Check if file is from ZetaView
    elif 'Zeta' in df.columns:
        std_df = pd.DataFrame({
            'Zeta Potential (mV)': df['Zeta'],
        })
        
    # Unknown format - try to guess based on column names
    else:
        # Look for zeta column
        zeta_cols = [col for col in df.columns if any(term in col.lower() 
                                                    for term in ['zeta', 'potential', 'mv'])]
        
        if zeta_cols:
            std_df = pd.DataFrame({
                'Zeta Potential (mV)': df[zeta_cols[0]],
            })
        else:
            raise ValueError("Could not determine the format of the zeta potential data file")
    
    # Save to output file if specified
    if output_file:
        save_ev_data(std_df, output_file)
        
    return std_df


def calculate_size_distribution(
    data: pd.DataFrame,
    size_column: str = 'Size (nm)',
    concentration_column: str = 'Concentration (particles/ml)',
    bins: Optional[List[float]] = None,
    normalize: bool = True
) -> Dict[str, np.ndarray]:
    """
    Calculate size distribution from particle data.
    
    Args:
        data: DataFrame containing size and concentration data
        size_column: Name of the column containing size values
        concentration_column: Name of the column containing concentration values
        bins: Optional list of bin edges for size distribution
        normalize: Whether to normalize the distribution (default: True)
        
    Returns:
        Dictionary with 'bins', 'counts', and 'distribution' keys
    """
    # Extract size and concentration data
    sizes = data[size_column].values
    concentrations = data[concentration_column].values
    
    # Create bins if not provided
    if bins is None:
        min_size = np.floor(np.min(sizes) / 10) * 10
        max_size = np.ceil(np.max(sizes) / 10) * 10
        bins = np.arange(min_size, max_size + 10, 10)  # 10 nm bin width
    
    # Calculate histogram
    hist, bin_edges = np.histogram(sizes, bins=bins, weights=concentrations)
    
    # Normalize if requested
    if normalize:
        hist = hist / np.sum(hist)
    
    # Return as dictionary
    return {
        'bins': bin_edges,
        'bin_centers': (bin_edges[:-1] + bin_edges[1:]) / 2,
        'counts': hist,
        'distribution': hist / np.sum(hist) if not normalize else hist
    }


def merge_ev_data(
    file_paths: List[str],
    output_file: Optional[str] = None
) -> pd.DataFrame:
    """
    Merge multiple EV data files into a single DataFrame.
    
    Args:
        file_paths: List of paths to EV data files
        output_file: Optional path to save the merged data
        
    Returns:
        Merged DataFrame
    """
    dfs = []
    
    for file_path in file_paths:
        df = load_ev_data(file_path)
        
        # Add source file info
        df['Source_File'] = os.path.basename(file_path)
        
        dfs.append(df)
    
    # Merge dataframes
    merged_df = pd.concat(dfs, ignore_index=True)
    
    # Save to output file if specified
    if output_file:
        save_ev_data(merged_df, output_file)
    
    return merged_df


def calculate_membrane_potential(
    zeta_potential: float,
    size_nm: float,
    ionic_strength_mM: float = 150.0,
    temperature_C: float = 25.0,
    model: str = 'gouy_chapman'
) -> float:
    """
    Estimate membrane potential from zeta potential.
    
    Args:
        zeta_potential: Zeta potential in mV
        size_nm: Vesicle diameter in nm
        ionic_strength_mM: Ionic strength in mM
        temperature_C: Temperature in Celsius
        model: Model to use ('gouy_chapman' or 'basic')
        
    Returns:
        Estimated membrane potential in mV
    """
    # Constants
    R = 8.314  # Gas constant, J/(K·mol)
    F = 96485  # Faraday constant, C/mol
    
    # Convert temperature to Kelvin
    T = temperature_C + 273.15
    
    # Calculate Debye length (nm)
    debye_length = 0.304 / np.sqrt(ionic_strength_mM)
    
    if model == 'gouy_chapman':
        # Calculate surface potential
        surface_potential = zeta_potential * np.exp(size_nm / (2 * debye_length))
        
        # Estimate membrane potential
        membrane_potential = surface_potential - 40  # Empirical offset
    else:
        # Basic approximation
        membrane_potential = zeta_potential * 1.5
    
    return membrane_potential


def convert_dep_crossover_to_properties(
    crossover_frequency_Hz: float,
    medium_conductivity_Sm: float,
    vesicle_size_nm: float
) -> Dict[str, float]:
    """
    Estimate EV properties from DEP crossover frequency.
    
    Args:
        crossover_frequency_Hz: Crossover frequency in Hz
        medium_conductivity_Sm: Medium conductivity in S/m
        vesicle_size_nm: Vesicle diameter in nm
        
    Returns:
        Dictionary of estimated properties
    """
    # Constants
    epsilon_0 = 8.85e-12  # Vacuum permittivity, F/m
    medium_permittivity = 78.5  # Relative permittivity of water
    
    # Calculate crossover frequency in rad/s
    omega_crossover = 2 * np.pi * crossover_frequency_Hz
    
    # Convert size to meters
    radius_m = vesicle_size_nm * 1e-9 / 2
    
    # Estimate membrane capacitance (F/m²)
    Cmem = np.sqrt(2) * medium_conductivity_Sm / (radius_m * omega_crossover)
    
    # Estimate membrane thickness (assuming typical permittivity)
    membrane_permittivity = 8.0
    membrane_thickness_nm = 1000 * membrane_permittivity * epsilon_0 / Cmem
    
    # Estimate effective conductivity at crossover
    effective_conductivity = medium_conductivity_Sm
    
    # Estimate membrane potential (very approximate)
    potential_factor = crossover_frequency_Hz / 1e6
    membrane_potential = -70 * potential_factor
    
    return {
        'Membrane_Capacitance_F_m2': Cmem,
        'Membrane_Thickness_nm': membrane_thickness_nm,
        'Effective_Conductivity_S_m': effective_conductivity,
        'Estimated_Membrane_Potential_mV': membrane_potential
    }


# Example usage
if __name__ == "__main__":
    # Example of loading and saving EV data
    try:
        # This is just an example and won't run without actual files
        data = load_ev_data("example_data.csv")
        print(f"Loaded data with {len(data)} rows")
        
        # Calculate size distribution
        dist = calculate_size_distribution(data)
        print(f"Size distribution calculated with {len(dist['bins'])-1} bins")
        
        # Estimate membrane potential from zeta potential
        zeta = -25.0  # Example zeta potential in mV
        size = 100.0  # Example vesicle size in nm
        membrane_pot = calculate_membrane_potential(zeta, size)
        print(f"Estimated membrane potential: {membrane_pot:.2f} mV")
        
        # Save processed data
        save_ev_data(data, "processed_data.csv")
        print("Data saved successfully")
        
    except FileNotFoundError:
        print("This is an example script. Please provide actual file paths to use.")
