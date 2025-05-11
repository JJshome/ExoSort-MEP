#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Characterization of Extracellular Vesicles

This module provides tools for analyzing and characterizing extracellular vesicles
based on various properties, with a focus on membrane potential differences.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Optional, Union
from dataclasses import dataclass
from scipy.stats import gaussian_kde
import matplotlib.patches as mpatches


@dataclass
class EVFraction:
    """Class for storing extracellular vesicle fraction information."""
    name: str
    identifier: str
    membrane_potential_range: Tuple[float, float]  # mV
    concentration: float  # particles/mL
    size_distribution: Dict[float, float]  # nm: frequency
    protein_content: Dict[str, float]  # protein: abundance
    nucleic_acid_content: Dict[str, float]  # RNA/DNA type: concentration (ng/mL)
    lipid_content: Dict[str, float]  # lipid type: percentage
    source: str
    sorting_method: str
    sorting_parameters: Dict[str, float]
    
    def __str__(self):
        return (f"EVFraction: {self.name} ({self.identifier})\n"
                f"  Membrane Potential Range: {self.membrane_potential_range[0]} to {self.membrane_potential_range[1]} mV\n"
                f"  Concentration: {self.concentration:.2e} particles/mL\n"
                f"  Size Range: {min(self.size_distribution.keys()):.1f} to {max(self.size_distribution.keys()):.1f} nm\n"
                f"  Source: {self.source}\n"
                f"  Sorting Method: {self.sorting_method}")
    
    @property
    def mean_size(self) -> float:
        """Calculate the mean vesicle size from the distribution."""
        sizes = np.array(list(self.size_distribution.keys()))
        frequencies = np.array(list(self.size_distribution.values()))
        return np.sum(sizes * frequencies) / np.sum(frequencies)
    
    @property
    def median_size(self) -> float:
        """Calculate the median vesicle size from the distribution."""
        sizes = np.array(list(self.size_distribution.keys()))
        frequencies = np.array(list(self.size_distribution.values()))
        
        # Calculate cumulative distribution
        cum_frequencies = np.cumsum(frequencies)
        cum_frequencies /= cum_frequencies[-1]  # Normalize
        
        # Find median size (where cumulative frequency crosses 0.5)
        median_idx = np.searchsorted(cum_frequencies, 0.5)
        
        return sizes[median_idx]
    
    @property
    def mode_size(self) -> float:
        """Find the most frequent vesicle size."""
        sizes = np.array(list(self.size_distribution.keys()))
        frequencies = np.array(list(self.size_distribution.values()))
        return sizes[np.argmax(frequencies)]
    
    @property
    def size_range(self) -> Tuple[float, float]:
        """Return the range of vesicle sizes."""
        sizes = list(self.size_distribution.keys())
        return (min(sizes), max(sizes))
    
    def size_percentile(self, percentile: float) -> float:
        """
        Calculate the vesicle size at a given percentile.
        
        Args:
            percentile: Percentile (0-100)
            
        Returns:
            Size at specified percentile
        """
        if percentile < 0 or percentile > 100:
            raise ValueError("Percentile must be between 0 and 100")
            
        sizes = np.array(list(self.size_distribution.keys()))
        frequencies = np.array(list(self.size_distribution.values()))
        
        # Calculate cumulative distribution
        cum_frequencies = np.cumsum(frequencies)
        cum_frequencies /= cum_frequencies[-1]  # Normalize
        
        # Find size at percentile
        perc_idx = np.searchsorted(cum_frequencies, percentile / 100.0)
        
        return sizes[min(perc_idx, len(sizes) - 1)]


class EVAnalyzer:
    """
    Class for analyzing and characterizing extracellular vesicle fractions.
    
    This class provides tools for comparing and visualizing different EV fractions,
    particularly those sorted by membrane potential difference.
    """
    
    def __init__(self):
        """Initialize the EV analyzer."""
        self.fractions = []
        
    def add_fraction(self, fraction: EVFraction) -> None:
        """
        Add an EV fraction to the analyzer.
        
        Args:
            fraction: EVFraction object to add
        """
        self.fractions.append(fraction)
        print(f"Added fraction: {fraction.name}")
        
    def get_fraction(self, identifier: str) -> Optional[EVFraction]:
        """
        Retrieve a fraction by its identifier.
        
        Args:
            identifier: Unique identifier for the fraction
            
        Returns:
            Matching EVFraction or None if not found
        """
        for fraction in self.fractions:
            if fraction.identifier == identifier:
                return fraction
        return None
    
    def compare_size_distributions(self, 
                                  identifiers: Optional[List[str]] = None,
                                  kde: bool = True,
                                  filename: Optional[str] = None) -> None:
        """
        Compare size distributions of EV fractions.
        
        Args:
            identifiers: List of fraction identifiers to compare (None for all)
            kde: Whether to use kernel density estimation for smoothing
            filename: Optional filename to save the plot
        """
        if not self.fractions:
            raise ValueError("No fractions added. Use add_fraction() first.")
            
        # Select fractions to compare
        if identifiers:
            fractions = [f for f in self.fractions if f.identifier in identifiers]
            if not fractions:
                raise ValueError(f"No fractions found with the given identifiers: {identifiers}")
        else:
            fractions = self.fractions
            
        # Create figure
        plt.figure(figsize=(10, 6))
        
        # Color map based on membrane potential
        potentials = [np.mean(f.membrane_potential_range) for f in fractions]
        norm = plt.Normalize(min(potentials), max(potentials))
        cmap = plt.cm.viridis
        
        # Plot each distribution
        for i, fraction in enumerate(fractions):
            sizes = np.array(list(fraction.size_distribution.keys()))
            frequencies = np.array(list(fraction.size_distribution.values()))
            mean_potential = np.mean(fraction.membrane_potential_range)
            color = cmap(norm(mean_potential))
            
            if kde:
                # Use KDE for smoothing
                x_grid = np.linspace(min(sizes) * 0.8, max(sizes) * 1.2, 1000)
                
                # Calculate weights based on frequencies
                weights = frequencies / np.sum(frequencies)
                
                # Generate samples for KDE
                samples = []
                for size, weight in zip(sizes, weights):
                    n_samples = int(weight * 10000)
                    samples.extend([size] * n_samples)
                
                samples = np.array(samples)
                kde_model = gaussian_kde(samples)
                density = kde_model(x_grid)
                
                plt.plot(x_grid, density, label=f"{fraction.name} ({mean_potential:.1f} mV)", 
                       color=color, linewidth=2)
            else:
                # Plot raw histogram
                plt.plot(sizes, frequencies, label=f"{fraction.name} ({mean_potential:.1f} mV)", 
                       color=color, marker='o', linestyle='-', markersize=4)
                
            # Add mean, median, mode markers
            plt.axvline(x=fraction.mean_size, color=color, linestyle='--', alpha=0.7)
            plt.axvline(x=fraction.median_size, color=color, linestyle=':', alpha=0.7)
            plt.axvline(x=fraction.mode_size, color=color, linestyle='-.', alpha=0.7)
            
        # Add legend items for statistics
        handles, labels = plt.gca().get_legend_handles_labels()
        handles.extend([
            mpatches.Patch(color='gray', linestyle='--', label='Mean'),
            mpatches.Patch(color='gray', linestyle=':', label='Median'),
            mpatches.Patch(color='gray', linestyle='-.', label='Mode')
        ])
        
        # Add labels and title
        plt.xlabel('Vesicle Size (nm)')
        plt.ylabel('Density' if kde else 'Frequency')
        plt.title('Size Distribution of EV Fractions')
        plt.legend(handles=handles, loc='best')
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Add colorbar for membrane potential
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, orientation='vertical', pad=0.01)
        cbar.set_label('Mean Membrane Potential (mV)')
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()
    
    def compare_protein_content(self, 
                               identifiers: Optional[List[str]] = None,
                               top_n: int = 10,
                               filename: Optional[str] = None) -> None:
        """
        Compare protein content of EV fractions.
        
        Args:
            identifiers: List of fraction identifiers to compare (None for all)
            top_n: Number of top proteins to display
            filename: Optional filename to save the plot
        """
        if not self.fractions:
            raise ValueError("No fractions added. Use add_fraction() first.")
            
        # Select fractions to compare
        if identifiers:
            fractions = [f for f in self.fractions if f.identifier in identifiers]
            if not fractions:
                raise ValueError(f"No fractions found with the given identifiers: {identifiers}")
        else:
            fractions = self.fractions
            
        # Combine all proteins across fractions
        all_proteins = set()
        for fraction in fractions:
            all_proteins.update(fraction.protein_content.keys())
            
        # Get the top N proteins by average abundance
        protein_avg_abundance = {}
        for protein in all_proteins:
            abundances = [f.protein_content.get(protein, 0) for f in fractions]
            protein_avg_abundance[protein] = np.mean(abundances)
            
        top_proteins = sorted(protein_avg_abundance.items(), 
                            key=lambda x: x[1], reverse=True)[:top_n]
        top_protein_names = [p[0] for p in top_proteins]
        
        # Create dataframe for plotting
        data = []
        for fraction in fractions:
            mean_potential = np.mean(fraction.membrane_potential_range)
            for protein in top_protein_names:
                abundance = fraction.protein_content.get(protein, 0)
                data.append({
                    'Fraction': fraction.name,
                    'Membrane Potential (mV)': mean_potential,
                    'Protein': protein,
                    'Abundance': abundance
                })
                
        df = pd.DataFrame(data)
        
        # Create figure
        plt.figure(figsize=(12, 8))
        
        # Use seaborn for better heatmap
        pivot_df = df.pivot(index='Protein', columns='Fraction', values='Abundance')
        
        # Add membrane potential to column names
        potentials = [np.mean(f.membrane_potential_range) for f in fractions]
        column_labels = [f"{f.name}\n({p:.1f} mV)" for f, p in zip(fractions, potentials)]
        pivot_df.columns = column_labels
        
        # Sort proteins by average abundance
        pivot_df = pivot_df.reindex(top_protein_names)
        
        # Plot heatmap
        sns.heatmap(pivot_df, annot=True, cmap='viridis', fmt='.2f', 
                  linewidths=0.5, cbar_kws={'label': 'Protein Abundance'})
        
        plt.title('Protein Content Comparison Across EV Fractions')
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()
    
    def compare_nucleic_acid_content(self, 
                                    identifiers: Optional[List[str]] = None,
                                    filename: Optional[str] = None) -> None:
        """
        Compare nucleic acid content of EV fractions.
        
        Args:
            identifiers: List of fraction identifiers to compare (None for all)
            filename: Optional filename to save the plot
        """
        if not self.fractions:
            raise ValueError("No fractions added. Use add_fraction() first.")
            
        # Select fractions to compare
        if identifiers:
            fractions = [f for f in self.fractions if f.identifier in identifiers]
            if not fractions:
                raise ValueError(f"No fractions found with the given identifiers: {identifiers}")
        else:
            fractions = self.fractions
            
        # Combine all nucleic acid types across fractions
        all_types = set()
        for fraction in fractions:
            all_types.update(fraction.nucleic_acid_content.keys())
            
        # Create dataframe for plotting
        data = []
        for fraction in fractions:
            mean_potential = np.mean(fraction.membrane_potential_range)
            for acid_type in all_types:
                concentration = fraction.nucleic_acid_content.get(acid_type, 0)
                data.append({
                    'Fraction': fraction.name,
                    'Membrane Potential (mV)': mean_potential,
                    'Nucleic Acid Type': acid_type,
                    'Concentration (ng/mL)': concentration
                })
                
        df = pd.DataFrame(data)
        
        # Create figure
        plt.figure(figsize=(12, 6))
        
        # Color map based on membrane potential
        potentials = [np.mean(f.membrane_potential_range) for f in fractions]
        norm = plt.Normalize(min(potentials), max(potentials))
        cmap = plt.cm.viridis
        colors = [cmap(norm(p)) for p in potentials]
        
        # Bar positions
        n_fractions = len(fractions)
        n_types = len(all_types)
        width = 0.8 / n_fractions
        acid_types = sorted(all_types)
        
        # Plot grouped bars
        for i, fraction in enumerate(fractions):
            mean_potential = np.mean(fraction.membrane_potential_range)
            values = [fraction.nucleic_acid_content.get(acid, 0) for acid in acid_types]
            
            x = np.arange(len(acid_types))
            offset = width * i - width * n_fractions / 2 + width / 2
            
            plt.bar(x + offset, values, width, 
                  label=f"{fraction.name} ({mean_potential:.1f} mV)", 
                  color=colors[i], edgecolor='black', linewidth=0.5)
            
        # Add labels and title
        plt.xlabel('Nucleic Acid Type')
        plt.ylabel('Concentration (ng/mL)')
        plt.title('Nucleic Acid Content Comparison Across EV Fractions')
        plt.xticks(np.arange(len(acid_types)), acid_types, rotation=45, ha='right')
        plt.legend(loc='best')
        plt.grid(True, axis='y', linestyle='--', alpha=0.7)
        
        # Add colorbar for membrane potential
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, orientation='vertical', pad=0.01)
        cbar.set_label('Mean Membrane Potential (mV)')
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()
    
    def compare_lipid_content(self, 
                             identifiers: Optional[List[str]] = None,
                             filename: Optional[str] = None) -> None:
        """
        Compare lipid content of EV fractions.
        
        Args:
            identifiers: List of fraction identifiers to compare (None for all)
            filename: Optional filename to save the plot
        """
        if not self.fractions:
            raise ValueError("No fractions added. Use add_fraction() first.")
            
        # Select fractions to compare
        if identifiers:
            fractions = [f for f in self.fractions if f.identifier in identifiers]
            if not fractions:
                raise ValueError(f"No fractions found with the given identifiers: {identifiers}")
        else:
            fractions = self.fractions
            
        # Get all lipid types
        all_lipids = set()
        for fraction in fractions:
            all_lipids.update(fraction.lipid_content.keys())
            
        # Create figure with subplots - one pie chart per fraction
        n_fractions = len(fractions)
        fig, axes = plt.subplots(1, n_fractions, figsize=(5 * n_fractions, 6))
        
        # Handle case of single fraction
        if n_fractions == 1:
            axes = [axes]
            
        # Plot pie charts
        for i, (fraction, ax) in enumerate(zip(fractions, axes)):
            mean_potential = np.mean(fraction.membrane_potential_range)
            
            # Get lipid percentages
            lipids = sorted(all_lipids)
            values = [fraction.lipid_content.get(lipid, 0) for lipid in lipids]
            
            # Use a colormap based on values
            cmap = plt.cm.viridis
            colors = [cmap(i/len(lipids)) for i in range(len(lipids))]
            
            # Create pie chart
            wedges, texts, autotexts = ax.pie(
                values, 
                labels=None,
                autopct='%1.1f%%', 
                startangle=90, 
                colors=colors
            )
            
            # Customize text
            for autotext in autotexts:
                autotext.set_fontsize(8)
                autotext.set_color('white')
                
            # Add title
            ax.set_title(f"{fraction.name}\n({mean_potential:.1f} mV)")
            
        # Add a single legend for all subplots
        handles = [mpatches.Patch(color=cmap(i/len(lipids)), label=lipid) 
                 for i, lipid in enumerate(lipids)]
        fig.legend(handles=handles, loc='lower center', ncol=min(5, len(lipids)), 
                 bbox_to_anchor=(0.5, 0.05))
        
        plt.suptitle('Lipid Content Comparison Across EV Fractions', fontsize=16)
        plt.tight_layout(rect=[0, 0.1, 1, 0.95])  # Make room for legend
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()
    
    def membrane_potential_correlation(self, 
                                      property_name: str,
                                      filename: Optional[str] = None) -> None:
        """
        Analyze correlation between membrane potential and a specific property.
        
        Args:
            property_name: Name of the property to correlate ('size', 'protein', 'nucleic_acid', 'lipid')
            filename: Optional filename to save the plot
        """
        if not self.fractions:
            raise ValueError("No fractions added. Use add_fraction() first.")
            
        potentials = [np.mean(f.membrane_potential_range) for f in self.fractions]
        
        if property_name == 'size':
            # Correlate with mean vesicle size
            values = [f.mean_size for f in self.fractions]
            y_label = 'Mean Vesicle Size (nm)'
            title = 'Correlation Between Membrane Potential and Vesicle Size'
            
        elif property_name == 'protein':
            # Calculate total protein content
            values = [sum(f.protein_content.values()) for f in self.fractions]
            y_label = 'Total Protein Content'
            title = 'Correlation Between Membrane Potential and Protein Content'
            
        elif property_name == 'nucleic_acid':
            # Calculate total nucleic acid content
            values = [sum(f.nucleic_acid_content.values()) for f in self.fractions]
            y_label = 'Total Nucleic Acid Content (ng/mL)'
            title = 'Correlation Between Membrane Potential and Nucleic Acid Content'
            
        elif property_name == 'lipid':
            # We need a specific lipid type to correlate
            # Use first lipid type that's present in all fractions
            all_lipids = set.intersection(*[set(f.lipid_content.keys()) for f in self.fractions])
            
            if not all_lipids:
                raise ValueError("No common lipid type found across all fractions")
                
            lipid_type = sorted(all_lipids)[0]
            values = [f.lipid_content[lipid_type] for f in self.fractions]
            y_label = f'{lipid_type} Content (%)'
            title = f'Correlation Between Membrane Potential and {lipid_type} Content'
            
        else:
            raise ValueError(f"Unknown property name: {property_name}")
            
        # Create figure
        plt.figure(figsize=(10, 6))
        
        # Create scatter plot
        plt.scatter(potentials, values, c=potentials, cmap='viridis', 
                  s=100, edgecolors='k', alpha=0.7)
        
        # Add fraction names as labels
        for i, fraction in enumerate(self.fractions):
            plt.annotate(fraction.name, (potentials[i], values[i]), 
                       xytext=(5, 5), textcoords='offset points')
            
        # Calculate and plot regression line if we have enough points
        if len(self.fractions) > 2:
            z = np.polyfit(potentials, values, 1)
            p = np.poly1d(z)
            x_line = np.linspace(min(potentials), max(potentials), 100)
            plt.plot(x_line, p(x_line), 'r--', alpha=0.7)
            
            # Calculate correlation coefficient
            corr = np.corrcoef(potentials, values)[0, 1]
            plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', 
                   transform=plt.gca().transAxes, fontsize=12,
                   bbox=dict(facecolor='white', alpha=0.7))
            
        # Add labels and title
        plt.xlabel('Membrane Potential (mV)')
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(min(potentials), max(potentials)))
        sm.set_array([])
        cbar = plt.colorbar(sm, orientation='vertical', pad=0.01)
        cbar.set_label('Membrane Potential (mV)')
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"Figure saved to {filename}")
            
        plt.show()
    
    def generate_summary_report(self, 
                               identifiers: Optional[List[str]] = None,
                               filename: Optional[str] = None) -> pd.DataFrame:
        """
        Generate a summary report of EV fractions.
        
        Args:
            identifiers: List of fraction identifiers to include (None for all)
            filename: Optional filename to save the report
            
        Returns:
            DataFrame with summary statistics
        """
        if not self.fractions:
            raise ValueError("No fractions added. Use add_fraction() first.")
            
        # Select fractions to include
        if identifiers:
            fractions = [f for f in self.fractions if f.identifier in identifiers]
            if not fractions:
                raise ValueError(f"No fractions found with the given identifiers: {identifiers}")
        else:
            fractions = self.fractions
            
        # Create summary dataframe
        data = []
        for fraction in fractions:
            # Basic properties
            row = {
                'Fraction Name': fraction.name,
                'Identifier': fraction.identifier,
                'Min Membrane Potential (mV)': fraction.membrane_potential_range[0],
                'Max Membrane Potential (mV)': fraction.membrane_potential_range[1],
                'Mean Membrane Potential (mV)': np.mean(fraction.membrane_potential_range),
                'Concentration (particles/mL)': fraction.concentration,
                'Mean Size (nm)': fraction.mean_size,
                'Median Size (nm)': fraction.median_size,
                'Mode Size (nm)': fraction.mode_size,
                'Min Size (nm)': fraction.size_range[0],
                'Max Size (nm)': fraction.size_range[1],
                'Source': fraction.source,
                'Sorting Method': fraction.sorting_method
            }
            
            # Add sorting parameters
            for param, value in fraction.sorting_parameters.items():
                row[f'Sorting Parameter: {param}'] = value
                
            # Add top proteins
            top_proteins = sorted(fraction.protein_content.items(), 
                               key=lambda x: x[1], reverse=True)[:3]
            for i, (protein, abundance) in enumerate(top_proteins):
                row[f'Top Protein {i+1}'] = protein
                row[f'Top Protein {i+1} Abundance'] = abundance
                
            # Add nucleic acid totals
            row['Total Nucleic Acid (ng/mL)'] = sum(fraction.nucleic_acid_content.values())
            
            # Add data to list
            data.append(row)
            
        # Create DataFrame
        df = pd.DataFrame(data)
        
        # Save to file if specified
        if filename:
            df.to_csv(filename, index=False)
            print(f"Summary report saved to {filename}")
            
        return df


# Example usage:
if __name__ == "__main__":
    # Create sample EV fractions with different membrane potentials
    fractions = [
        EVFraction(
            name="Hyperpolarized EVs",
            identifier="hyper_ev_001",
            membrane_potential_range=(-100, -80),
            concentration=1e9,
            size_distribution={
                80: 0.05, 90: 0.15, 100: 0.30, 110: 0.25, 120: 0.15, 130: 0.05, 140: 0.05
            },
            protein_content={
                "CD9": 0.85, "CD63": 0.92, "CD81": 0.90, 
                "TSG101": 0.75, "Alix": 0.80, "Flotillin-1": 0.65,
                "HSP70": 0.70, "Syntenin-1": 0.60, "Annexin A2": 0.55,
                "GAPDH": 0.40, "Actin": 0.45
            },
            nucleic_acid_content={
                "miRNA": 25.0, "mRNA": 15.5, "lncRNA": 8.2, "DNA": 2.5
            },
            lipid_content={
                "Phosphatidylcholine": 45.0, "Phosphatidylserine": 15.0,
                "Phosphatidylethanolamine": 20.0, "Cholesterol": 15.0,
                "Sphingomyelin": 5.0
            },
            source="Neural cells",
            sorting_method="DC Electrophoresis",
            sorting_parameters={
                "Voltage": 120.0, "Buffer pH": 7.4, "Run Time (min)": 60.0
            }
        ),
        EVFraction(
            name="Normal EVs",
            identifier="normal_ev_001",
            membrane_potential_range=(-75, -55),
            concentration=3e9,
            size_distribution={
                80: 0.10, 90: 0.20, 100: 0.25, 110: 0.20, 120: 0.15, 130: 0.05, 140: 0.05
            },
            protein_content={
                "CD9": 0.80, "CD63": 0.85, "CD81": 0.82, 
                "TSG101": 0.70, "Alix": 0.75, "Flotillin-1": 0.60,
                "HSP70": 0.65, "Syntenin-1": 0.50, "Annexin A2": 0.60,
                "GAPDH": 0.45, "Actin": 0.50
            },
            nucleic_acid_content={
                "miRNA": 20.0, "mRNA": 18.0, "lncRNA": 7.5, "DNA": 3.0
            },
            lipid_content={
                "Phosphatidylcholine": 40.0, "Phosphatidylserine": 12.0,
                "Phosphatidylethanolamine": 22.0, "Cholesterol": 18.0,
                "Sphingomyelin": 8.0
            },
            source="Fibroblasts",
            sorting_method="DC Electrophoresis",
            sorting_parameters={
                "Voltage": 120.0, "Buffer pH": 7.4, "Run Time (min)": 60.0
            }
        ),
        EVFraction(
            name="Depolarized EVs",
            identifier="depol_ev_001",
            membrane_potential_range=(-50, -20),
            concentration=2e9,
            size_distribution={
                80: 0.10, 90: 0.15, 100: 0.20, 110: 0.25, 120: 0.20, 130: 0.05, 140: 0.05
            },
            protein_content={
                "CD9": 0.70, "CD63": 0.75, "CD81": 0.72, 
                "TSG101": 0.60, "Alix": 0.65, "Flotillin-1": 0.55,
                "HSP70": 0.60, "Syntenin-1": 0.45, "Annexin A2": 0.65,
                "GAPDH": 0.50, "Actin": 0.55
            },
            nucleic_acid_content={
                "miRNA": 15.0, "mRNA": 20.0, "lncRNA": 6.5, "DNA": 4.0
            },
            lipid_content={
                "Phosphatidylcholine": 35.0, "Phosphatidylserine": 10.0,
                "Phosphatidylethanolamine": 25.0, "Cholesterol": 20.0,
                "Sphingomyelin": 10.0
            },
            source="Cancer cells",
            sorting_method="DC Electrophoresis",
            sorting_parameters={
                "Voltage": 120.0, "Buffer pH": 7.4, "Run Time (min)": 60.0
            }
        ),
        EVFraction(
            name="Positive EVs",
            identifier="pos_ev_001",
            membrane_potential_range=(5, 15),
            concentration=5e8,
            size_distribution={
                80: 0.15, 90: 0.20, 100: 0.15, 110: 0.15, 120: 0.15, 130: 0.10, 140: 0.10
            },
            protein_content={
                "CD9": 0.60, "CD63": 0.65, "CD81": 0.62, 
                "TSG101": 0.50, "Alix": 0.55, "Flotillin-1": 0.50,
                "HSP70": 0.55, "Syntenin-1": 0.40, "Annexin A2": 0.70,
                "GAPDH": 0.55, "Actin": 0.60
            },
            nucleic_acid_content={
                "miRNA": 10.0, "mRNA": 25.0, "lncRNA": 5.0, "DNA": 6.0
            },
            lipid_content={
                "Phosphatidylcholine": 30.0, "Phosphatidylserine": 8.0,
                "Phosphatidylethanolamine": 28.0, "Cholesterol": 22.0,
                "Sphingomyelin": 12.0
            },
            source="Activated immune cells",
            sorting_method="DC Electrophoresis",
            sorting_parameters={
                "Voltage": 120.0, "Buffer pH": 7.4, "Run Time (min)": 60.0
            }
        )
    ]
    
    # Create analyzer
    analyzer = EVAnalyzer()
    
    # Add fractions
    for fraction in fractions:
        analyzer.add_fraction(fraction)
    
    # Compare size distributions
    analyzer.compare_size_distributions(kde=True, filename="size_distributions.png")
    
    # Compare protein content
    analyzer.compare_protein_content(top_n=8, filename="protein_content.png")
    
    # Compare nucleic acid content
    analyzer.compare_nucleic_acid_content(filename="nucleic_acid_content.png")
    
    # Compare lipid content
    analyzer.compare_lipid_content(filename="lipid_content.png")
    
    # Analyze correlations
    for prop in ['size', 'protein', 'nucleic_acid', 'lipid']:
        analyzer.membrane_potential_correlation(prop, filename=f"{prop}_correlation.png")
    
    # Generate summary report
    report = analyzer.generate_summary_report(filename="ev_fractions_summary.csv")
    print("\nSummary Report:")
    print(report.to_string())
