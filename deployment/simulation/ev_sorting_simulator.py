#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ExoSort-MEP Simulation Tool

This module provides a web-based simulation environment for exploring the
electrophoresis and dielectrophoresis sorting of extracellular vesicles based on
membrane potential differences.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import panel as pn
import holoviews as hv
from holoviews import opts
import param
from bokeh.models import HoverTool
import io
import base64

# Initialize Panel and HoloViews
pn.extension('bokeh')
hv.extension('bokeh')

# Constants
DEFAULT_VESICLE_TYPES = [
    {"name": "Hyperpolarized EVs", "membrane_potential": -90.0, "size": 100, "concentration": 1e9},
    {"name": "Normal EVs", "membrane_potential": -70.0, "size": 120, "concentration": 5e9},
    {"name": "Depolarized EVs", "membrane_potential": -30.0, "size": 130, "concentration": 3e9},
    {"name": "Positive EVs", "membrane_potential": 10.0, "size": 90, "concentration": 2e8}
]

# Color mapping based on membrane potential
def get_color(potential):
    """Get color based on membrane potential."""
    if potential <= -80:
        return "#1976D2"  # Blue for hyperpolarized
    elif potential <= -50:
        return "#7B1FA2"  # Purple for normal negative
    elif potential <= 0:
        return "#E53935"  # Red for depolarized
    else:
        return "#FF9800"  # Orange for positive


class EVSortingSimulator(param.Parameterized):
    """Interactive simulator for EV sorting based on membrane potential."""
    
    # DC Electrophoresis Parameters
    voltage = param.Number(120.0, bounds=(10, 300), step=10, doc="Applied voltage (V)")
    buffer_ph = param.Number(7.4, bounds=(5.0, 9.0), step=0.1, doc="Buffer pH")
    buffer_conductivity = param.Number(1.5, bounds=(0.1, 5.0), step=0.1, doc="Buffer conductivity (mS/cm)")
    run_time = param.Number(60.0, bounds=(10, 240), step=10, doc="Run time (minutes)")
    
    # AC Dielectrophoresis Parameters
    frequency = param.Number(1000000.0, bounds=(1000, 1000000000), step=1000, doc="AC Frequency (Hz)")
    ac_voltage = param.Number(10.0, bounds=(1, 30), step=1, doc="AC Voltage (V)")
    medium_conductivity = param.Number(0.05, bounds=(0.001, 0.5), step=0.001, doc="Medium conductivity (S/m)")
    
    # Simulation type
    simulation_type = param.Selector(objects=["DC Electrophoresis", "AC Dielectrophoresis"], 
                                   default="DC Electrophoresis")
    
    # Button to run simulation
    run_simulation_button = param.Action(lambda self: self.run_simulation(), label="Run Simulation")
    
    def __init__(self, **params):
        super().__init__(**params)
        self.vesicle_types = DEFAULT_VESICLE_TYPES.copy()
        self.simulation_results = None
        self._create_ui()
        
    def _create_ui(self):
        """Create the user interface for the simulator."""
        # Main layout
        self.main_layout = pn.Column(
            pn.pane.Markdown("# ExoSort-MEP Simulation Tool", style={'text-align': 'center'}),
            pn.pane.Markdown("## Extracellular Vesicle Sorting by Membrane Potential", style={'text-align': 'center'}),
            
            pn.Row(
                # Left column - Parameters and controls
                pn.Column(
                    pn.pane.Markdown("### Simulation Parameters"),
                    pn.Param(self.param.simulation_type, name="Simulation Type"),
                    
                    pn.pane.Markdown("#### DC Electrophoresis Parameters", 
                                   visible=self.simulation_type == "DC Electrophoresis"),
                    pn.Param(self.param.voltage, name="Voltage (V)",
                           visible=self.simulation_type == "DC Electrophoresis"),
                    pn.Param(self.param.buffer_ph, name="Buffer pH",
                           visible=self.simulation_type == "DC Electrophoresis"),
                    pn.Param(self.param.buffer_conductivity, name="Buffer Conductivity (mS/cm)",
                           visible=self.simulation_type == "DC Electrophoresis"),
                    pn.Param(self.param.run_time, name="Run Time (minutes)",
                           visible=self.simulation_type == "DC Electrophoresis"),
                    
                    pn.pane.Markdown("#### AC Dielectrophoresis Parameters", 
                                   visible=self.simulation_type == "AC Dielectrophoresis"),
                    pn.Param(self.param.frequency, name="Frequency (Hz)",
                           visible=self.simulation_type == "AC Dielectrophoresis"),
                    pn.Param(self.param.ac_voltage, name="AC Voltage (V)",
                           visible=self.simulation_type == "AC Dielectrophoresis"),
                    pn.Param(self.param.medium_conductivity, name="Medium Conductivity (S/m)",
                           visible=self.simulation_type == "AC Dielectrophoresis"),
                    
                    pn.pane.Markdown("### Vesicle Properties"),
                    self._create_vesicle_table(),
                    
                    pn.Param(self.param.run_simulation_button, name="Run Simulation"),
                    width=400
                ),
                
                # Right column - Visualization
                pn.Column(
                    pn.pane.Markdown("### Simulation Results"),
                    pn.pane.Markdown("Click 'Run Simulation' to see results."),
                    name="results_column",
                    width=600
                )
            )
        )
        
        # Update UI visibility when simulation type changes
        self.param.watch(self._update_ui_visibility, "simulation_type")
        
    def _create_vesicle_table(self):
        """Create an editable table for vesicle properties."""
        vesicle_df = pd.DataFrame(self.vesicle_types)
        
        def edit_vesicles(event):
            self.vesicle_types = event.new.to_dict('records')
            
        table = pn.widgets.Tabulator(
            vesicle_df, 
            pagination='local', 
            page_size=10,
            widths={'name': 150, 'membrane_potential': 100, 'size': 80, 'concentration': 120},
            editors={
                'name': 'text',
                'membrane_potential': 'number',
                'size': 'number',
                'concentration': 'number'
            },
            formatters={
                'membrane_potential': {'precision': 1},
                'size': {'precision': 0},
                'concentration': {'type': 'scientific', 'precision': 2}
            },
            titles={
                'name': 'Name',
                'membrane_potential': 'Membrane Potential (mV)',
                'size': 'Size (nm)',
                'concentration': 'Concentration (particles/mL)'
            },
            height=200
        )
        
        table.param.watch(edit_vesicles, 'value')
        
        return table
    
    def _update_ui_visibility(self, *events):
        """Update UI element visibility based on simulation type."""
        for widget in self.main_layout[1][0]:
            if isinstance(widget, pn.Param) or isinstance(widget, pn.pane.Markdown):
                if "DC Electrophoresis" in widget.name and self.simulation_type == "DC Electrophoresis":
                    widget.visible = True
                elif "AC Dielectrophoresis" in widget.name and self.simulation_type == "AC Dielectrophoresis":
                    widget.visible = True
                elif "DC Electrophoresis" in widget.name and self.simulation_type != "DC Electrophoresis":
                    widget.visible = False
                elif "AC Dielectrophoresis" in widget.name and self.simulation_type != "AC Dielectrophoresis":
                    widget.visible = False
    
    def run_simulation(self):
        """Run the selected simulation type."""
        if self.simulation_type == "DC Electrophoresis":
            self._run_dc_simulation()
        else:
            self._run_ac_simulation()
            
        # Update results display
        self._update_results_display()
    
    def _run_dc_simulation(self):
        """Run DC electrophoresis simulation."""
        # Calculate mobility based on membrane potential for each vesicle type
        results = []
        
        for vesicle in self.vesicle_types:
            # Mobility calculation (simplified model)
            # More negative membrane potential -> higher mobility towards positive pole
            mobility = -0.2 * vesicle["membrane_potential"]  # µm·cm/V·s
            
            # Field strength
            field_strength = self.voltage / 100.0  # V/cm (assumes 10 cm gel)
            
            # Calculate velocity
            velocity = mobility * field_strength  # µm/s
            
            # Convert to mm/min for easier interpretation
            velocity_mm_min = velocity * 60 / 1000
            
            # Calculate distance traveled
            distance = velocity_mm_min * self.run_time  # mm
            
            # Determine collection fraction (simplified)
            if distance < 5:
                fraction = "Fraction 1 (0-5 mm)"
            elif distance < 10:
                fraction = "Fraction 2 (5-10 mm)"
            elif distance < 15:
                fraction = "Fraction 3 (10-15 mm)"
            elif distance < 20:
                fraction = "Fraction 4 (15-20 mm)"
            else:
                fraction = "Fraction 5 (>20 mm)"
            
            # Add to results
            results.append({
                "Name": vesicle["name"],
                "Membrane_Potential_mV": vesicle["membrane_potential"],
                "Size_nm": vesicle["size"],
                "Concentration": vesicle["concentration"],
                "Mobility": mobility,
                "Velocity_mm_min": velocity_mm_min,
                "Distance_mm": distance,
                "Fraction": fraction,
                "Color": get_color(vesicle["membrane_potential"])
            })
        
        self.simulation_results = pd.DataFrame(results)
    
    def _run_ac_simulation(self):
        """Run AC dielectrophoresis simulation."""
        # Calculate CM factor and DEP force based on membrane potential for each vesicle type
        results = []
        
        for vesicle in self.vesicle_types:
            # Simplified Clausius-Mossotti factor calculation
            # This is a simplified model for educational purposes
            # Real CM factor depends on complex permittivity of vesicle and medium
            potential = vesicle["membrane_potential"]
            
            # Frequency-dependent scaling factor
            freq_factor = 1.0 / (1.0 + (self.frequency / 1e6)**0.5)
            
            # Effect of membrane potential on CM factor
            # More negative potential -> stronger polarization
            if potential < 0:
                cm_factor = -0.5 * (abs(potential) / 100.0) * freq_factor
            else:
                cm_factor = 0.3 * (potential / 100.0) * freq_factor
            
            # For very high frequencies, CM factor tends toward 0.5
            if self.frequency > 1e8:
                cm_factor = 0.5 - (0.5 - cm_factor) * (1e9 - self.frequency) / 1e9
            
            # Calculate DEP force (simplified)
            # Force depends on volume, CM factor, and gradient of electric field squared
            volume = 4/3 * np.pi * (vesicle["size"] * 1e-9 / 2)**3  # m³
            gradient_e_squared = (self.ac_voltage / 50e-6)**2 / 50e-6  # V²/m³ (assumed electrode gap of 50 µm)
            
            # Medium permittivity
            medium_permittivity = 78.5 * 8.85e-12  # F/m
            
            # DEP force
            dep_force = 2 * np.pi * medium_permittivity * volume * cm_factor * gradient_e_squared  # N
            
            # Determine collection fraction based on CM factor
            if cm_factor < -0.3:
                fraction = "Strong Negative DEP"
                direction = "Low field regions"
            elif cm_factor < 0:
                fraction = "Weak Negative DEP"
                direction = "Low field regions"
            elif cm_factor < 0.3:
                fraction = "Weak Positive DEP"
                direction = "High field regions"
            else:
                fraction = "Strong Positive DEP"
                direction = "High field regions"
            
            # Add to results
            results.append({
                "Name": vesicle["name"],
                "Membrane_Potential_mV": vesicle["membrane_potential"],
                "Size_nm": vesicle["size"],
                "Concentration": vesicle["concentration"],
                "CM_Factor": cm_factor,
                "DEP_Force_N": dep_force,
                "Fraction": fraction,
                "Direction": direction,
                "Color": get_color(vesicle["membrane_potential"])
            })
        
        self.simulation_results = pd.DataFrame(results)
    
    def _update_results_display(self):
        """Update the results display with visualizations."""
        if self.simulation_results is None:
            return
        
        # Clear previous results
        results_column = self.main_layout[1][1]
        results_column.clear()
        results_column.append(pn.pane.Markdown("### Simulation Results"))
        
        # Add results table
        results_table = pn.widgets.Tabulator(
            self.simulation_results.drop(columns=['Color']), 
            pagination='local',
            page_size=10,
            formatters={
                'Membrane_Potential_mV': {'precision': 1},
                'Size_nm': {'precision': 0},
                'Concentration': {'type': 'scientific', 'precision': 2},
                'Mobility': {'precision': 4} if 'Mobility' in self.simulation_results.columns else None,
                'Velocity_mm_min': {'precision': 4} if 'Velocity_mm_min' in self.simulation_results.columns else None,
                'Distance_mm': {'precision': 2} if 'Distance_mm' in self.simulation_results.columns else None,
                'CM_Factor': {'precision': 4} if 'CM_Factor' in self.simulation_results.columns else None,
                'DEP_Force_N': {'type': 'scientific', 'precision': 2} if 'DEP_Force_N' in self.simulation_results.columns else None
            }
        )
        results_column.append(results_table)
        
        # Add visualizations based on simulation type
        if self.simulation_type == "DC Electrophoresis":
            results_column.append(self._create_dc_visualization())
        else:
            results_column.append(self._create_ac_visualization())
    
    def _create_dc_visualization(self):
        """Create visualization for DC electrophoresis results."""
        df = self.simulation_results
        
        # Create scatter plot of membrane potential vs distance
        scatter = hv.Points(
            df, kdims=['Membrane_Potential_mV', 'Distance_mm'], 
            vdims=['Name', 'Size_nm', 'Fraction', 'Color']
        ).opts(
            color='Color',
            size=10,
            marker='o',
            padding=0.1,
            width=500,
            height=400,
            xlabel='Membrane Potential (mV)',
            ylabel='Distance (mm)',
            title='Membrane Potential vs. Distance Traveled',
            tools=['hover'],
            toolbar='above'
        )
        
        # Create bar chart of vesicle distribution by fraction
        fraction_counts = df.groupby('Fraction')['Concentration'].sum().reset_index()
        fraction_counts['Concentration_normalized'] = fraction_counts['Concentration'] / fraction_counts['Concentration'].sum() * 100
        
        bars = hv.Bars(
            fraction_counts, kdims=['Fraction'], vdims=['Concentration_normalized']
        ).opts(
            color='Fraction',
            cmap='Category10',
            width=500,
            height=300,
            xlabel='Fraction',
            ylabel='Percentage (%)',
            title='Vesicle Distribution by Fraction',
            xrotation=45,
            tools=['hover'],
            toolbar='above'
        )
        
        # Create gel visualization
        gel_width = 10  # cm
        gel_height = 25  # mm (max distance + margin)
        gel_padding = 2  # mm
        
        # Create gel rectangle
        gel_rect = hv.Rectangles([(0, 0, gel_width, gel_height)]).opts(
            color='lightblue', alpha=0.3, line_width=1, line_color='blue',
            xlabel='Width (cm)', ylabel='Distance (mm)',
            xlim=(-1, gel_width + 1), ylim=(-1, gel_height + 1),
            title='Gel Electrophoresis Visualization',
            width=500, height=400
        )
        
        # Create wells at the top of the gel
        wells = []
        well_width = 0.8  # cm
        well_height = 2  # mm
        n_wells = len(df)
        well_spacing = gel_width / (n_wells + 1)
        
        for i in range(n_wells):
            well_x = (i + 1) * well_spacing - well_width / 2
            wells.append((well_x, 0, well_x + well_width, well_height))
        
        well_rects = hv.Rectangles(wells).opts(color='black')
        
        # Create bands based on distances
        bands = []
        for i, row in df.iterrows():
            band_x = (i + 1) * well_spacing - well_width / 2
            band_y = row['Distance_mm'] - 1  # Offset by 1mm
            bands.append((band_x, band_y, band_x + well_width, band_y + 2))
        
        # Get colors for bands
        band_colors = df['Color'].tolist()
        
        band_rects = hv.Rectangles(bands).opts(color=band_colors)
        
        # Add collection intervals as horizontal lines
        intervals = [5, 10, 15, 20]
        interval_lines = []
        
        for interval in intervals:
            interval_lines.append(hv.HLine(interval).opts(
                color='red', line_dash='dashed', line_width=1
            ))
        
        # Combine gel components
        gel_vis = gel_rect * well_rects * band_rects
        for line in interval_lines:
            gel_vis = gel_vis * line
            
        # Add labels for wells and bands
        labels = []
        for i, row in df.iterrows():
            # Well label
            well_x = (i + 1) * well_spacing
            labels.append(hv.Text(well_x, -0.5, f"Well {i+1}"))
            
            # Band label
            band_x = (i + 1) * well_spacing
            band_y = row['Distance_mm'] + 0.5
            labels.append(hv.Text(band_x, band_y, row['Name']))
        
        # Add interval labels
        for interval in intervals:
            labels.append(hv.Text(0.5, interval + 0.5, f"{interval} mm", halign='left'))
        
        label_overlay = hv.Overlay(labels).opts(text_font_size='8pt')
        gel_vis = gel_vis * label_overlay
        
        # Return visualizations in a layout
        return pn.Column(
            scatter,
            pn.Row(bars, gel_vis)
        )
    
    def _create_ac_visualization(self):
        """Create visualization for AC dielectrophoresis results."""
        df = self.simulation_results
        
        # Create scatter plot of membrane potential vs CM factor
        scatter_cm = hv.Points(
            df, kdims=['Membrane_Potential_mV', 'CM_Factor'], 
            vdims=['Name', 'Size_nm', 'Fraction', 'Direction', 'Color']
        ).opts(
            color='Color',
            size=10,
            marker='o',
            padding=0.1,
            width=500,
            height=350,
            xlabel='Membrane Potential (mV)',
            ylabel='Clausius-Mossotti Factor',
            title='Membrane Potential vs. CM Factor',
            tools=['hover'],
            toolbar='above'
        )
        
        # Add a horizontal line at CM factor = 0
        cm_zero_line = hv.HLine(0).opts(color='red', line_dash='dashed')
        scatter_cm = scatter_cm * cm_zero_line
        
        # Create bar chart of DEP force by vesicle type
        bars_dep = hv.Bars(
            df, kdims=['Name'], vdims=['DEP_Force_N']
        ).opts(
            color='Color',
            width=500,
            height=300,
            xlabel='Vesicle Type',
            ylabel='DEP Force (N)',
            title='DEP Force by Vesicle Type',
            xrotation=45,
            tools=['hover'],
            toolbar='above'
        )
        
        # Create diagram of DEP sorting
        # Create electric field lines visual
        x = np.linspace(0, 10, 100)
        y = np.linspace(0, 10, 100)
        
        # Create electric field gradient representation
        x_grid, y_grid = np.meshgrid(x, y)
        center_x, center_y = 5, 5
        distance = np.sqrt((x_grid - center_x)**2 + (y_grid - center_y)**2)
        field_strength = 1 / (distance + 0.5)
        
        # Normalize field strength
        field_strength = (field_strength - field_strength.min()) / (field_strength.max() - field_strength.min())
        
        # Convert to DataFrame for HoloViews
        field_df = pd.DataFrame({
            'x': x_grid.flatten(),
            'y': y_grid.flatten(),
            'strength': field_strength.flatten()
        })
        
        # Create heatmap of field gradient
        field_heatmap = hv.Points(
            field_df, kdims=['x', 'y'], vdims=['strength']
        ).opts(
            color='strength',
            cmap='viridis',
            colorbar=True,
            width=500,
            height=400,
            xlabel='',
            ylabel='',
            title='DEP Sorting Visualization',
            tools=['hover'],
            toolbar='above',
            clim=(0, 1)
        )
        
        # Add electrodes
        electrodes = []
        electrodes.append(hv.Rectangles([(4, 0, 6, 0.5)]).opts(color='black'))  # Bottom electrode
        electrodes.append(hv.Rectangles([(4, 9.5, 6, 10)]).opts(color='black'))  # Top electrode
        electrodes.append(hv.Rectangles([(0, 4, 0.5, 6)]).opts(color='black'))  # Left electrode
        electrodes.append(hv.Rectangles([(9.5, 4, 10, 6)]).opts(color='black'))  # Right electrode
        
        electrode_overlay = hv.Overlay(electrodes)
        
        # Add vesicle positions based on CM factor
        vesicle_positions = []
        
        for _, row in df.iterrows():
            cm_factor = row['CM_Factor']
            name = row['Name']
            color = row['Color']
            
            # Position based on CM factor
            if cm_factor < 0:  # Negative DEP (low field regions)
                # Corners (low field)
                pos_x = 1 + np.random.rand() * 2
                pos_y = 1 + np.random.rand() * 2
                if np.random.rand() > 0.5:
                    pos_x = 9 - np.random.rand() * 2
                if np.random.rand() > 0.5:
                    pos_y = 9 - np.random.rand() * 2
            else:  # Positive DEP (high field regions)
                # Edges (high field)
                if np.random.rand() > 0.5:
                    pos_x = 5 + np.random.rand() * 0.5 - 0.25
                    pos_y = (0.5 if np.random.rand() > 0.5 else 9.5) + np.random.rand() * 0.5 - 0.25
                else:
                    pos_x = (0.5 if np.random.rand() > 0.5 else 9.5) + np.random.rand() * 0.5 - 0.25
                    pos_y = 5 + np.random.rand() * 0.5 - 0.25
            
            vesicle_positions.append((pos_x, pos_y, name, color, cm_factor))
        
        vesicle_df = pd.DataFrame(vesicle_positions, columns=['x', 'y', 'Name', 'Color', 'CM_Factor'])
        
        vesicle_points = hv.Points(
            vesicle_df, kdims=['x', 'y'], vdims=['Name', 'Color', 'CM_Factor']
        ).opts(
            color='Color',
            size=15,
            marker='o',
            tools=['hover']
        )
        
        # Combine field visualization with electrodes and vesicles
        dep_vis = field_heatmap * electrode_overlay * vesicle_points
        
        # Create prediction table
        prediction_df = df[['Name', 'Membrane_Potential_mV', 'CM_Factor', 'Direction']].copy()
        prediction_table = pn.widgets.Tabulator(
            prediction_df,
            pagination='local',
            page_size=10,
            formatters={
                'Membrane_Potential_mV': {'precision': 1},
                'CM_Factor': {'precision': 4}
            }
        )
        
        # Return visualizations in a layout
        return pn.Column(
            scatter_cm,
            pn.Row(bars_dep, dep_vis),
            pn.pane.Markdown("### DEP Force Predictions"),
            prediction_table
        )
    
    def get_app(self):
        """Get the Panel app for the simulator."""
        return self.main_layout
    
    def save_results(self, filename=None):
        """
        Save the simulation results to a CSV file.
        
        Args:
            filename: Optional filename for the CSV file (default: auto-generated)
        """
        if self.simulation_results is None:
            print("No simulation results to save. Run a simulation first.")
            return
            
        if filename is None:
            sim_type = "DC" if self.simulation_type == "DC Electrophoresis" else "AC"
            filename = f"evsort_{sim_type}_simulation_results.csv"
            
        # Save to CSV
        self.simulation_results.drop(columns=['Color']).to_csv(filename, index=False)
        print(f"Results saved to {filename}")


def main():
    """Run the simulator as a standalone app."""
    simulator = EVSortingSimulator()
    app = simulator.get_app()
    
    # Serve the app
    pn.serve(app, port=5006)


if __name__ == "__main__":
    main()
