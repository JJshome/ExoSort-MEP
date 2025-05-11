# Protocol for Extracellular Vesicle Sorting Based on Membrane Potential

This protocol describes the methods for sorting extracellular vesicles (EVs) based on their membrane potential differences using both DC electrophoresis and AC dielectrophoresis techniques.

## 1. DC Electrophoresis Protocol

### 1.1 Materials and Equipment

- **EV Sample Preparation**
  - Isolated extracellular vesicles (≥ 1×10⁹ particles/mL)
  - PBS buffer (pH 7.4)
  - Ultracentrifuge (100,000 × g capability)
  - 0.22 μm filters

- **Electrophoresis Setup**
  - Horizontal electrophoresis chamber
  - Agarose (0.5-1.0%)
  - TAE buffer (Tris-acetate-EDTA, pH 7.4)
  - DC power supply (0-300V)
  - Micropipettes and tips

- **Collection Equipment**
  - Gel cutter or razor blade
  - Microcentrifuge tubes
  - Gel extraction buffer
  - Electroelution equipment (optional)
  - Eppendorf tubes

### 1.2 EV Sample Preparation

1. **Isolation of EVs**
   - Obtain EVs from cell culture supernatant or biological fluids using standard isolation techniques (differential ultracentrifugation, size exclusion chromatography, or commercial kits).

2. **Quality Control**
   - Measure concentration and size distribution using nanoparticle tracking analysis (NTA) or dynamic light scattering (DLS).
   - Confirm EV markers (CD63, CD9, CD81) by Western blot or flow cytometry.
   - Measure initial zeta potential using zetasizer.

3. **Sample Preparation for Electrophoresis**
   - Resuspend EVs in TAE buffer at concentration of 1×10¹⁰ particles/mL.
   - Add 10% glycerol to prevent sample diffusion.
   - If needed, stain with lipophilic dyes for visualization (e.g., DiO, DiI).

### 1.3 Electrophoresis Procedure

1. **Gel Preparation**
   - Prepare 0.8% agarose gel in TAE buffer.
   - Cast gel with large loading wells (50-100 μL capacity).
   - Allow gel to solidify for 30 minutes at room temperature.

2. **Sample Loading**
   - Place gel in electrophoresis chamber filled with TAE buffer.
   - Load 50-100 μL of prepared EV sample into wells.
   - Include appropriate size/charge markers in separate lanes.

3. **Electrophoresis Run**
   - Connect power supply to electrophoresis chamber.
   - Set voltage to 100-120V (or 5-10 V/cm).
   - Run for 60-90 minutes at room temperature.
   - Monitor run progress visually if samples were stained.

4. **Fractionation and Collection**
   - Mark positions on gel at 5mm intervals from the well.
   - Cut gel at these marked positions using a sterile gel cutter.
   - Place each gel section in a separate tube containing 500 μL extraction buffer.
   - Incubate at 4°C overnight with gentle shaking to extract EVs.
   - Alternatively, use electroelution to recover EVs from gel sections.

5. **Recovery of EVs from Gel Fragments**
   - Centrifuge tubes at 10,000 × g for 10 minutes.
   - Collect supernatant containing EVs.
   - Filter through 0.22 μm filter to remove gel particles.
   - Concentrate EVs using ultracentrifugation (100,000 × g for 70 minutes).
   - Resuspend EV pellet in PBS or desired buffer.

### 1.4 Quality Control and Characterization

1. Measure concentration of collected fractions using NTA or DLS.
2. Verify size distribution of EVs in each fraction.
3. Measure zeta potential of each fraction.
4. Perform transmission electron microscopy (TEM) for morphological analysis.
5. Analyze protein content using Western blot or mass spectrometry.
6. Quantify nucleic acid content using qPCR or RNA-seq.

## 2. AC Dielectrophoresis Protocol

### 2.1 Materials and Equipment

- **EV Sample Preparation**
  - Isolated extracellular vesicles (≥ 1×10⁹ particles/mL)
  - Low conductivity buffer (typically 0.05-0.1 S/m)
  - Sucrose for density adjustment (optional)
  - 0.22 μm filters

- **Dielectrophoresis Setup**
  - Microfluidic DEP device with electrode arrays
  - Function generator (AC, 10 kHz - 50 MHz range)
  - Amplifier (optional)
  - Oscilloscope for signal verification
  - Syringe pump or pressure control system
  - Microscope for visualization (fluorescence capability preferred)

- **Collection Equipment**
  - Outlet tubing
  - Microcentrifuge tubes
  - Micropipettes and tips

### 2.2 EV Sample Preparation

1. **Buffer Preparation**
   - Prepare low conductivity buffer (e.g., 8.5% sucrose + 0.3% dextrose in deionized water).
   - Adjust conductivity to 0.05-0.1 S/m by adding minimal amounts of PBS.
   - Filter buffer through 0.22 μm filter.
   - Adjust pH to 7.4.

2. **Sample Preparation**
   - Resuspend EVs in prepared buffer at concentration of 1×10⁹ particles/mL.
   - If needed, stain with lipophilic fluorescent dyes for visualization.
   - Filter through 0.22 μm filter to remove aggregates.

### 2.3 DEP Device Preparation

1. **Device Assembly**
   - Clean microfluidic channels and electrode surfaces.
   - Assemble microfluidic device according to manufacturer's instructions.
   - Connect input and output tubing.
   - If using custom device, ensure proper sealing to prevent leakage.

2. **System Priming**
   - Prime all tubing and channels with buffer solution.
   - Ensure no air bubbles are present in the system.
   - Establish stable flow conditions before introducing sample.

3. **Electrical Connections**
   - Connect function generator to device electrodes.
   - Verify connections and signal parameters using oscilloscope.
   - Set initial frequency to 1 MHz for testing.

### 2.4 DEP Sorting Procedure

1. **Initial System Test**
   - Apply AC signal (1-10 V peak-to-peak) at 1 MHz.
   - Verify proper device operation with buffer only.
   - Check for any electrode damage or abnormal heating.

2. **Optimization**
   - Determine optimal frequency for separation.
   - Test frequency sweeps from 100 kHz to 10 MHz.
   - Identify crossover frequencies where DEP force changes direction.
   - Optimize voltage amplitude (typically 5-20 Vpp).

3. **Sample Processing**
   - Load EV sample into syringe or reservoir.
   - Set flow rate to 0.5-2 μL/min (optimize based on device dimensions).
   - Apply optimized AC signal (frequency and voltage).
   - Collect fractions from different outlet channels.

4. **Multiple Frequency Operation (Optional)**
   - For more refined separation, implement frequency shifting.
   - Start with high frequency (e.g., 5 MHz) for initial separation.
   - Switch to lower frequency (e.g., 500 kHz) for secondary separation.
   - Collect resulting sub-fractions separately.

### 2.5 Fraction Collection and Processing

1. **Collection**
   - Collect fractions from different outlet channels in separate tubes.
   - Keep samples on ice during collection.
   - Monitor separation using microscopy if EVs are fluorescently labeled.

2. **Post-Processing**
   - Concentrate collected fractions by ultracentrifugation (100,000 × g for 70 minutes).
   - Resuspend EV pellets in PBS or desired buffer for downstream applications.
   - Store at 4°C for short-term or -80°C for long-term storage.

### 2.6 Quality Control and Characterization

1. Measure concentration of collected fractions using NTA or DLS.
2. Verify size distribution of EVs in each fraction.
3. Measure membrane potential using potentiometric dyes (e.g., DiBAC4(3)).
4. Perform transmission electron microscopy (TEM) for morphological analysis.
5. Analyze protein content using Western blot or mass spectrometry.
6. Assess functionality using appropriate biological assays.

## 3. Troubleshooting

### 3.1 DC Electrophoresis Issues

| Problem | Possible Causes | Solutions |
|---------|-----------------|-----------|
| Poor EV migration | Buffer conductivity too high | Dilute buffer or use lower ionic strength buffer |
| | Voltage too low | Increase voltage (up to 10 V/cm) |
| | Gel percentage too high | Reduce agarose percentage to 0.5-0.6% |
| EV aggregation | Sample concentration too high | Dilute sample before loading |
| | Buffer incompatibility | Optimize buffer conditions |
| Low recovery | EVs stuck in gel matrix | Use electroelution for recovery |
| | Loss during filtration | Use low-protein-binding filters |
| Poor separation | Run time too short | Increase electrophoresis time |
| | Voltage too high | Reduce voltage to prevent heating |

### 3.2 AC Dielectrophoresis Issues

| Problem | Possible Causes | Solutions |
|---------|-----------------|-----------|
| Weak DEP response | Buffer conductivity too high | Reduce conductivity to 0.05 S/m |
| | Frequency not optimal | Perform frequency sweep to find optimal frequency |
| | Voltage too low | Increase voltage (within safe limits) |
| Electrode damage | Voltage too high | Reduce applied voltage |
| | Electrolysis | Use higher frequencies (>100 kHz) |
| | Buffer conductivity too high | Use lower conductivity buffer |
| Poor separation | Flow rate too high | Reduce flow rate |
| | Inadequate electrode design | Optimize electrode geometry |
| | Improper frequency | Identify crossover frequencies and optimize |
| EV adherence to surfaces | Surface properties | Add 0.1% BSA to buffer to reduce adhesion |
| | Excessive DEP force | Reduce voltage |

## 4. Notes and Recommendations

1. **Sample Quality**
   - Use freshly isolated EVs whenever possible.
   - Minimize freeze-thaw cycles which can alter membrane properties.
   - Verify EV integrity before and after sorting.

2. **Optimization**
   - Both methods require optimization for specific EV populations.
   - Perform pilot runs with small sample volumes before processing valuable samples.
   - Document all parameters for reproducibility.

3. **Controls**
   - Include synthetic liposomes of known sizes and charges as controls.
   - Process buffer-only controls to identify potential contamination.
   - Use EVs with characterized membrane potentials as references.

4. **Biological Validation**
   - Verify that the sorting process does not compromise EV biological activity.
   - Compare sorted fractions in functional assays.
   - Assess target cell uptake and biological effects for each fraction.

## 5. Safety Considerations

1. **Electrical Safety**
   - Follow all electrical safety protocols when working with high voltage equipment.
   - Ensure proper grounding of all equipment.
   - Turn off power supplies before making any connections or adjustments.

2. **Biological Safety**
   - Handle all biological materials according to appropriate biosafety level.
   - Use appropriate personal protective equipment.
   - Properly dispose of biological waste according to institutional guidelines.

3. **Chemical Safety**
   - Handle all chemicals according to safety data sheets.
   - Use appropriate ventilation when working with volatile reagents.
   - Wear appropriate protective equipment.

## 6. References

1. Tayebi M, Yang D, Collins DJ, Ai Y. Deterministic Sorting of Submicrometer Particles and Extracellular Vesicles Using a Combined Electric and Acoustic Field. Nano Lett. 2021;21(16):6935-6941.

2. Nakano M, Nakabayashi R, Koyama R, Inaba M, Suehiro J. Dielectrophoresis Low and High Crossover Frequencies of Cancerous Exosomes. Electrophoresis. 2025;46(4):623-631.

3. Park M, Lee CH, Noh H, et al. High-precision extracellular-vesicle isolation-analysis integrated platform for rapid cancer diagnosis directly from blood plasma. Biosens Bioelectron. 2025;225:116863.

4. Lan M, Wu D, Cheng C, et al. Small Extracellular Vesicles Detection by Dielectrophoresis-Based Microfluidic Chip Filled with Transparent Antibody-Conjugated Microbeads for Breast Cancer Diagnosis. Anal Chem. 2025;97(6):3142-3151.

5. Wang D, Yang S, Wang N, et al. A Novel Microfluidic Strategy for Efficient Exosome Separation via Thermally Oxidized Non-Uniform Deterministic Lateral Displacement (DLD) Arrays and Dielectrophoresis (DEP) Synergy. Biosensors. 2024;14(4):174.

---

**Protocol Version:** 1.0  
**Last Updated:** May 2025  
**Patent Status:** Patent Pending
