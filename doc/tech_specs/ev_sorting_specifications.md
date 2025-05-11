# Technical Specifications for ExoSort-MEP System

This document outlines the technical specifications for the ExoSort-MEP system, which sorts extracellular vesicles (EVs) based on their membrane potential differences using both DC electrophoresis and AC dielectrophoresis methods.

## 1. System Overview

ExoSort-MEP is a comprehensive platform for the separation and classification of extracellular vesicles based on their membrane potential differences. The system integrates two primary sorting technologies:

1. **DC Electrophoresis** - Utilizes direct current to separate vesicles based on charge magnitude
2. **AC Dielectrophoresis** - Employs alternating current to create non-uniform electric fields for separation

The platform includes hardware components, control software, and analysis tools for characterizing the sorted populations.

## 2. DC Electrophoresis Module Specifications

### 2.1 Hardware Components

| Component | Specification | Notes |
|-----------|---------------|-------|
| **Electrophoresis Chamber** | High-resolution horizontal format | Specialized for nanoscale vesicles |
| **Power Supply** | 0-300V DC, 0-500mA, precision control ±0.5V | Digital display with real-time monitoring |
| **Cooling System** | Integrated Peltier, 4-25°C range | Prevents sample degradation |
| **Gel System** | Ultra-low porosity agarose (0.5-1.0%) | Specialized for EV applications |
| **Collection System** | 5 fraction collectors | Positioned at 5mm intervals |
| **Buffer Circulation** | Programmable flow rate (0-10 mL/min) | Maintains pH and conductivity |

### 2.2 Performance Parameters

| Parameter | Range | Resolution | Accuracy |
|-----------|-------|------------|----------|
| Operating Voltage | 10-300V | 0.1V | ±0.5V |
| Operating Current | 0-500mA | 0.1mA | ±0.5mA |
| Running Time | 1-240 minutes | 1 minute | ±0.5 minute |
| Temperature Control | 4-25°C | 0.1°C | ±0.5°C |
| pH Range | 5.0-9.0 | 0.1 pH unit | ±0.1 pH unit |
| Buffer Conductivity | 0.1-10 mS/cm | 0.01 mS/cm | ±0.05 mS/cm |

### 2.3 Sample Requirements

- **Volume**: 50-500 μL per run
- **Concentration**: 10⁸-10¹⁰ particles/mL
- **Buffer Compatibility**: TAE, TBE, PBS (with adjustment)
- **Sample Preparation**: Minimal (resuspension in appropriate buffer)

### 2.4 Expected Separation Performance

| Vesicle Type | Membrane Potential Range | Expected Migration Distance (60 min @ 120V) |
|--------------|--------------------------|-------------------------------------------|
| Hyperpolarized | -100 to -80 mV | 20-25 mm |
| Normal | -80 to -50 mV | 15-20 mm |
| Depolarized | -50 to -20 mV | 10-15 mm |
| Near-neutral | -20 to 0 mV | 5-10 mm |
| Positive | >0 mV | 0-5 mm |

## 3. AC Dielectrophoresis Module Specifications

### 3.1 Hardware Components

| Component | Specification | Notes |
|-----------|---------------|-------|
| **Microfluidic Chip** | Quadrupole electrode array | Field gradient optimization |
| **Signal Generator** | 10 Hz - 50 MHz, 0-20 Vpp | Precision frequency control |
| **Amplifier** | 20W, 0-30 Vpp output | Low distortion (<0.1% THD) |
| **Flow System** | 0.1-10 μL/min, programmable | Precision syringe pumps |
| **Collection System** | 4 outlet channels | Automated fraction collection |
| **Visualization** | Integrated fluorescence microscopy | Real-time monitoring |
| **Temperature Sensor** | Integrated thermistor array | Prevents sample heating |

### 3.2 Performance Parameters

| Parameter | Range | Resolution | Accuracy |
|-----------|-------|------------|----------|
| Frequency Range | 10 Hz - 50 MHz | 1 Hz | ±0.01% |
| Voltage Amplitude | 0-30 Vpp | 0.1 V | ±0.1 V |
| Flow Rate | 0.1-10 μL/min | 0.1 μL/min | ±0.05 μL/min |
| Operating Temperature | 4-37°C | 0.1°C | ±0.2°C |
| Medium Conductivity | 0.01-0.5 S/m | 0.001 S/m | ±0.005 S/m |
| Channel Dimensions | 200μm width, 50μm depth | - | ±2 μm |

### 3.3 Sample Requirements

- **Volume**: 20-200 μL per run
- **Concentration**: 10⁷-10⁹ particles/mL
- **Medium Compatibility**: Low conductivity buffers (0.01-0.1 S/m)
- **Sample Preparation**: Resuspension in DEP buffer required

### 3.4 Expected Separation Performance

| Vesicle Type | Membrane Potential Range | CM Factor at 1 MHz | DEP Behavior |
|--------------|--------------------------|-------------------|-------------|
| Hyperpolarized | -100 to -80 mV | -0.4 to -0.5 | Strong negative DEP |
| Normal | -80 to -50 mV | -0.2 to -0.4 | Moderate negative DEP |
| Depolarized | -50 to -20 mV | -0.1 to -0.2 | Weak negative DEP |
| Near-neutral | -20 to 0 mV | -0.1 to +0.1 | Crossover region |
| Positive | >0 mV | +0.1 to +0.3 | Positive DEP |

### 3.5 Crossover Frequency Ranges

| Vesicle Type | First Crossover | Second Crossover |
|--------------|-----------------|------------------|
| Hyperpolarized | 50-100 kHz | 10-20 MHz |
| Normal | 100-200 kHz | 5-10 MHz |
| Depolarized | 200-500 kHz | 2-5 MHz |
| Positive | 500-1000 kHz | 1-2 MHz |

## 4. Control System Specifications

### 4.1 Hardware Platform

- **Processor**: Embedded ARM Cortex-A53, quad-core 1.5 GHz
- **Memory**: 4 GB RAM, 64 GB storage
- **Display**: 10.1" capacitive touchscreen, 1920×1200 resolution
- **Connectivity**: Ethernet, Wi-Fi, USB 3.0, Bluetooth 5.0
- **Data Acquisition**: 16-bit ADC/DAC, 1 MS/s sampling rate
- **Power**: 110-240 VAC, 50/60 Hz, 200W max

### 4.2 Software Features

- **Operating System**: Real-time Linux
- **User Interface**: Intuitive touchscreen interface
- **Method Development**: Customizable protocol creation and storage
- **Data Analysis**: Real-time visualization and processing
- **Export Formats**: CSV, Excel, PDF, Images (PNG, TIFF)
- **Security**: User authentication and audit trails
- **Remote Access**: Optional network control and monitoring

### 4.3 Automation Capabilities

- **Protocol Library**: Pre-programmed methods for common vesicle types
- **Sequence Processing**: Batch processing of multiple samples
- **Adaptive Control**: Real-time parameter adjustment based on feedback
- **Fail-safe Features**: Error detection and safe shutdown procedures
- **Method Optimization**: Automated parameter optimization

## 5. Analysis Module Specifications

### 5.1 Integrated Analytics

| Capability | Specification | Notes |
|------------|---------------|-------|
| **Size Analysis** | 30-1000 nm range | Integration with NTA data |
| **Concentration Measurement** | 10⁶-10¹² particles/mL | Direct or dilution-based |
| **Membrane Potential Estimation** | -100 to +50 mV range | Based on migration pattern |
| **Zeta Potential Correlation** | -50 to +50 mV range | Calculated from mobility |
| **Protein Content Estimation** | Semi-quantitative | Based on standard curves |
| **Nucleic Acid Detection** | Presence/absence | Fluorescence-based detection |

### 5.2 Visualization and Reporting

- **Real-time Monitoring**: Live tracking of separation process
- **Migration Profiles**: Distance vs. time visualization
- **Fraction Analysis**: Composition reports by collection zone
- **Comparative Analysis**: Multi-sample overlay and comparison
- **Quality Metrics**: Reproducibility and separation efficiency calculations
- **Automated Reports**: Customizable templates with key parameters

## 6. Environmental and Operational Requirements

### 6.1 Operating Conditions

- **Temperature**: 15-30°C (instrument); 4-25°C (samples)
- **Humidity**: 20-80% non-condensing
- **Altitude**: Up to 2000m
- **Electrical**: 110-240V AC, 50/60 Hz, 10A max
- **Space Requirements**: Benchtop, 80 × 60 × 50 cm (W×D×H)
- **Weight**: Main unit 25 kg, accessories 5 kg

### 6.2 Regulatory Compliance

- **Laboratory Safety**: IEC 61010-1
- **EMC Standards**: IEC 61326-1
- **RoHS Compliant**: Yes
- **Quality System**: ISO 9001 manufacturing

### 6.3 Maintenance Requirements

- **Calibration Interval**: 12 months
- **Electrode Replacement**: 500 operating hours (AC-DEP)
- **Filter Replacement**: Monthly or every 100 runs
- **Software Updates**: Quarterly (web-based)
- **Preventive Maintenance**: Annual service recommended

## 7. Performance Validation

### 7.1 Calibration Standards

- **Size Standards**: NIST-traceable polystyrene beads (50, 100, 200, 500 nm)
- **Charge Standards**: Liposomes with defined surface charge
- **Membrane Potential Standards**: Synthetic vesicles with ionophores

### 7.2 System Validation Tests

| Test | Acceptance Criteria | Method |
|------|---------------------|--------|
| Separation Resolution | ≥ 90% pure fractions | NTA analysis of collected fractions |
| Run-to-Run Reproducibility | CV < 10% | Repeated runs of standard samples |
| Sample Recovery | > 80% recovery | Pre- and post-separation quantification |
| Cross-Contamination | < 5% between fractions | Fluorescent tracer testing |
| Temperature Control | ± 0.5°C of setpoint | Thermal mapping during operation |
| Flow Precision | ± 2% of setpoint | Volumetric collection and verification |

## 8. Consumables and Accessories

### 8.1 Standard Consumables

| Item | Specifications | Replacement Frequency |
|------|---------------|------------------------|
| **Running Buffers** | DC: TAE buffer kit; AC: DEP buffer kit | Each run |
| **Agarose** | Ultra-low porosity, EV-grade | Each DC run |
| **Collection Tubes** | 0.5 mL, low-protein binding | Each run |
| **Microfluidic Chips** | Disposable, pre-sterilized | Every 5 AC runs |
| **Filter Sets** | 0.22 μm, low protein binding | Each sample |
| **Calibration Kits** | Size and charge standards | Monthly |

### 8.2 Optional Accessories

- **Advanced Visualization Package**: High-resolution fluorescence imaging
- **Automated Sample Loader**: 24-sample capacity
- **Fraction Collector Upgrade**: 96-well format compatibility
- **Temperature Control Extension**: 4°C operation capability
- **Data Analysis Software Package**: Advanced statistics and multi-parameter analysis
- **Remote Monitoring System**: Cloud-based access and notifications

## 9. Integration Capabilities

### 9.1 Upstream Technologies

- **Sample Preparation**: Direct integration with EV isolation systems
- **Quality Control**: Compatibility with NTA, DLS, and flow cytometry data
- **Biomarker Screening**: Integration with proteomics and genomics workflows

### 9.2 Downstream Applications

- **Therapeutic Loading**: Direct transfer to drug loading platforms
- **Surface Modification**: Integration with EV engineering systems
- **Formulation**: Connection to stability and storage solutions
- **Assay Systems**: Direct transfer to functional testing platforms

### 9.3 Data Ecosystem

- **LIMS Integration**: Standard API for laboratory information systems
- **Cloud Storage**: Secure backup and remote access options
- **Collaborative Platform**: Multi-user access and annotation capabilities
- **Machine Learning**: Predictive maintenance and performance optimization

## 10. System Limitations and Considerations

### 10.1 Sample Limitations

- **Minimum Sample Volume**: 50 μL (DC), 20 μL (AC)
- **Minimum Concentration**: 10⁷ particles/mL
- **Sample Viscosity**: Must be < 5 cP
- **Particle Size Range**: 30-1000 nm optimal
- **Buffer Compatibility**: See detailed buffer guide in user manual

### 10.2 Performance Boundaries

- **Smallest Detectable Potential Difference**: ~10 mV
- **Maximum Throughput**: 4 samples/hour (DC), 2 samples/hour (AC)
- **Resolution Limit**: Cannot separately resolve potentials within ±5 mV
- **Continuous Operation**: Maximum 12 hours before maintenance
- **Sample Recovery**: Typically 60-90% depending on vesicle properties

---

**Document Version:** 1.0  
**Last Updated:** May 2025  
**Patent Status:** Patent Pending
