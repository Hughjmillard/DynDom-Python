# DynDom-Py

A Python implementation of the DynDom algorithm for analyzing protein domain movements between conformational states.

## Overview

DynDom-Py is a Python implementation of the DynDom algorithm for analyzing protein conformational changes in terms of rigid-body domain movements. The software determines protein domains, hinge axes, and amino acid residues involved in hinge bending through fully automated analysis.

Given two conformational states of the same protein (e.g., from X-ray crystallography, NMR, or molecular dynamics simulations), DynDom-Py identifies how the protein moves by treating different regions as quasi-rigid bodies. This approach transforms complex conformational changes into an easily understood view of domain movements connected by flexible hinges.

The analysis reveals:
- **Dynamic domains**: Regions that move as rigid bodies
- **Hinge regions**: Flexible areas that allow domain movement  
- **Screw axes**: Mathematical description of how domains rotate and translate relative to each other

## Features

- **Automated Domain Detection**: Identifies rigid domains using rotation vector clustering
- **Hinge Analysis**: Locates flexible regions and mechanical hinges between domains
- **Screw Axis Calculation**: Determines interdomain rotation axes, angles, and translations
- **Motion Validation**: Applies geometric constraints and interdomain/intradomain motion ratios
- **Comprehensive Visualization**: Generates PyMOL scripts with domain coloring and 3D motion arrows
- **Flexible Input**: Supports any two conformational states in PDB format
- **Detailed Output**: Provides motion statistics, hinge residue identification, and structural files
- **Quality Control**: Validates domain connectivity and motion significance

## Installation

### Requirements

- Python 3.7+
- Required packages:
  ```
  numpy
  scipy
  scikit-learn
  gemmi
  ```

### Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/dyndom-py.git
   cd dyndom-py
   ```

2. Install dependencies:
   ```bash
   pip install numpy scipy scikit-learn gemmi
   ```

3. Ensure PDB files are available in your input directory or the software will automatically download them from RCSB.

## Usage

### Basic Usage

1. **Prepare input files** in the `data/` directory:

   **command.txt**:
   ```
   input_path=data
   output_path=output
   filename1=file-name-here
   chain1id=chain-id-here
   filename2=file-name-here
   chain2id=chain-id-here
   ```

   **param.txt**:
   ```
   window=5
   domain=20
   ratio=1.0
   k_means_n_init=1
   k_means_max_iter=500
   atoms=backbone
   ```

2. **Run the analysis**:
   ```bash
   python main.py
   ```

### Parameters

#### Command Parameters (`command.txt`)
- `input_path`: Directory containing PDB files
- `output_path`: Directory for output files
- `filename1`, `filename2`: PDB IDs or filenames (without .pdb extension)
- `chain1id`, `chain2id`: Chain identifiers

#### Analysis Parameters (`param.txt`)
- `window`: Sliding window size for local motion analysis (default: 5)
- `domain`: Minimum domain size in residues (default: 20)
- `ratio`: Minimum ratio of interdomain to intradomain motion (default: 1.0)
- `atoms`: Atoms to use for analysis (`backbone` for N,CA,C or `ca` for CA only)
- `k_means_n_init`: Number of K-means initializations (default: 1)
- `k_means_max_iter`: Maximum K-means iterations (default: 500)

## Output Files

The software generates several output files in the specified output directory:

- **`{protein1}_{chain1}_{protein2}_{chain2}.pdb`**: Superimposed structures
- **`{protein1}_{chain1}_{protein2}_{chain2}.pml`**: PyMOL visualization script
- **`{protein1}_{chain1}_{protein2}_{chain2}_arrows.pdb`**: Motion arrows for PyMOL
- **`{protein1}_{chain1}_{protein2}_{chain2}.w5_info`**: Detailed analysis results
- **`{protein1}_rot_vecs.pdb`**: Rotation vectors for each residue

### Visualization

Load the generated `.pml` file in PyMOL to visualize:
- **Domain coloring**: Different domains shown in different colors
- **Fixed domain**: The biggest or most connected domain will be shown in blue
- **Hinge regions**: Flexible regions highlighted in green
- **Motion arrows**: 3D arrows showing screw axes and rotation directions
- **Arrow colouring**: Shaft shows fixed domain colour, tip shows dynamic domain


## Algorithm Overview

DynDom-Py implements a 7-step automated workflow for domain motion analysis:

1. **Global Structure Alignment**: Performs whole-protein best-fit superposition of the two conformational states
2. **Local Motion Detection**: Slides a window along the protein backbone to determine rotation vectors for short main-chain segments
3. **Rotation Vector Clustering**: Uses K-means clustering to identify groups of segments with similar rotational behavior
4. **Dynamic Domain Construction**: Builds domains from clustered segments, ensuring spatial connectivity (≤4.0 Å)
5. **Domain Validation**: Applies minimum size criteria and tests interdomain vs. intradomain motion ratios
6. **Screw Axis Calculation**: Determines interdomain screw axes based on Chasles' theorem for rigid body motion
7. **Hinge Residue Identification**: Locates residues involved in interdomain bending using statistical analysis of rotation distributions

### Key Concepts

- **Dynamic Domains**: Groups of residues that move as quasi-rigid bodies
- **Bending Residues**: Residues at domain boundaries with rotations outside the main distribution (P < 0.2)

## Methodology

The algorithm is based on the original DynDom method with several key steps:

## Citation

If you use DynDom-Py in your research, please cite:

```
[Your citation will be generated when published with DOI]
```

Original DynDom algorithm:
```
Hayward, S. and Berendsen, H.J.C. (1998) Systematic analysis of domain motions in proteins from conformational change: new results on citrate synthase and T4 lysozyme. Proteins 30, 144-154.
```


## Acknowledgments

- Original DynDom algorithm by Steven Hayward and Herman Berendsen
- GEMMI library for crystallographic computations
- PyMOL for molecular visualization
