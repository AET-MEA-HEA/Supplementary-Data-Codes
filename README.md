# Supplementary Data Codes 

**Three-dimensional atomic positions and local chemical order of medium-and high-entropy alloys**

Saman Moniri<sup>1*</sup>, Yao Yang<sup>1*</sup>, Jun Ding<sup>2</sup>, Yakun Yuan<sup>1</sup>, Jihan Zhou<sup>1</sup>, Long Yang<sup>1</sup>, Fan Zhu<sup>1</sup>, Yuxuan Liao<sup>1</sup>, Yonggang Yao<sup>3</sup>, Liangbing Hu<sup>3</sup>, Peter Ercius<sup>4</sup> & Jianwei Miao<sup>1†</sup>    

*<sup>1</sup>Department of Physics & Astronomy and California NanoSystems Institute, University of California, Los Angeles, CA 90095, USA.* 
*<sup>2</sup>Center for Alloy Innovation and Design, State Key Laboratory for Mechanical Behavior of Materials, Xi’an Jiaotong University, Xi’an, China.* 
*<sup>3</sup>Department of Materials Science and Engineering, University of Maryland, College Park, MD, USA.*      
*<sup>4</sup>National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA.*   
**These authors contributed equally to this work.*     
*†Correspondence and requests for materials should be addressed to J.M. (miao@physics.ucla.edu).*  

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)

# Overview

Medium- and high-entropy alloys (M/HEAs) mix multiple principal elements with near-equiatomic composition and represent a paradigm-shift strategy for designing new materials for metallurgy, catalysis, and other fields. One of the core hypotheses of M/HEAs is lattice distortion, which has been investigated by different numerical and experimental techniques. However, determining the three-dimensional (3D) lattice distortion in M/HEAs remains a challenge. Additionally, the presumed random elemental mixing in M/HEAs has been questioned by x-ray/neuron study, atomistic simulations, energy dispersive spectroscopy (EDS), and electron diffraction, which suggest the existence of local chemical order in M/HEAs. However, the 3D local chemical order has eluded direct experimental observation since the EDS integrates the composition of atomic columns along the zone axes and diffuse electron reflections may originate from planar defects instead of local chemical order. 

Here, we determine the 3D atomic positions of M/HEA nanoparticles using atomic electron tomography, and quantitatively characterize the local lattice distortion, strain tensor, twin boundaries, dislocation cores, and chemical short-range order (CSRO). We find that the HEAs have larger local lattice distortion and more heterogeneous strain than the MEAs and strain is correlated with CSRO. We also observe CSRO-mediated twinning in the MEAs, that is, twinning occurs in energetically unfavoured CSRO regions but not in energetically favoured CSRO ones, which represents the first experimental observation of correlating local chemical order with structural defects in any material. We expect that this work will not only expand our fundamental understanding of this important class of materials, but also could provide the foundation for tailoring M/HEA properties through engineering lattice distortion and local chemical order. 

# System Requirements

We recommend a computer with 16G DRAM, standard i7 4-core CPU, and a GPU to run most data analysis source codes. But for the 3D reconstruction of the experimental data with RESIRE, atomic tracing and the calculation of properties such as local displacement, strain tensor and chemical short-range order, we recommend a computer with large memory (256G DRAM, 16-core CPU and 1 GPU).

## Software Requirements

### OS Requirements

This package has been tested on the following Operating System:

Linux: CentOS 6 2.6.32    
Windows: Windows 11 22621.1702    
Mac OSX: We have not tested it on a Mac yet, but it should in principle work.     

### Matlab Version Requirements

This package has been tested with `Matlab` R2022a. All the codes have to run in their own folders. We recommend the use of `Matlab` version R2021a or higher to test the data and source codes.

# Repositary Contents

### 1. Experiment Data

Folder: [Measured_data](./1_Measured_data)

This folder contains the experimental projections after denoising and alignment as well as their corresponding angles for all the M/HEA nanoparticles.

### 2. Reconstructed 3D Volume

Folder: [Final_reconstruction_volume](./2_Final_reconstruction_volume)

This folder contains the 3D volume of all the M/HEA nanoparticles reconstructed from experimental projections and angles in [Measured_data](./1_Measured_data). 
The reconstruction codes can be found in [RESIRE_package](https://github.com/AET-MetallicGlass/Supplementary-Data-Codes/tree/master/2_RESIRE_package).

### 3. Experimental Atomic Model

Folder: [Final_coordinates](./3_Final_coordinates)

This folder contains the final 3D atomic models and chemical species (i.e. Ni Pd Pt or Type 1, 2, 3) of the M/HEA nanoparticles. 
The tracing codes can be found in [Tracing_codes](https://github.com/AET-MetallicGlass/Supplementary-Data-Codes/tree/master/4_Tracing_and_classification).

### 4. Post Data Analysis — The experimentally measured structural and chemical properties of the M/HEA nanoparticles

Folder: [Data_analysis_properties](./4_Data_analysis_properties)

Run the code `Main_1_calculate_local_displacement.m` to calculate the bond orientation order parameter and the average nearby Pt bond length for all the atoms in the M/HEA nanoparticles.   
Run the code `Main_2_plot_twin_op.m` to plot the twin order parameters for all the atoms in the M/HEA nanoparticles.    
Run the code `Main_3_strain_map.m` to plot the strain map of the M/HEA nanoparticles.    
Run the code `Main_4_calc_csro.m` to calculate the chemical short-range order of the M/HEA nanoparticles.    

