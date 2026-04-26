# Scaffold-Guided Protein Design for Three-Finger Toxin (3FTx) Binders
![forgithub](https://github.com/shantanuyadav19/de-Novo-Protein-Design/data/poster)

## Overview

This project implements a comprehensive scaffold-guided protein design pipeline targeting three-finger toxins (3FTx) from *Naja naja* (cobra) venom for computational antivenom development. The pipeline systematically explores multiple interface regions using **interface tensor conditioning** with RFdiffusion, followed by rigorous AlphaFold2-based validation.

> **Note**: This codebase is being prepared for integration into the official [RFdiffusion repository](https://github.com/RosettaCommons/RFdiffusion) and [dl_binder_design](https://github.com/nrbennet/dl_binder_design/tree/main).

## Project Structure

### Core Pipeline Scripts
```
├── make_interface_tensor.py            # Core tensor generation for scaffold guidance
├── generate_interface_tensors_*.py     # Interface tensor generation for specific regions
├── run_rfd_designs_*.py               # RFdiffusion inference with scaffold guidance
├── generate_consensus_sequence.py      # MSA-based consensus sequence generation
├── generate_scores_plots.py           # Comprehensive AlphaFold2 analysis & visualization
├── pdbinfo_reset.py                   # PDB renumbering utility for Rosetta
└── requirements.txt                   # Python dependencies
```

### Data and Results
```
├── uniprotkb_venom_neurotoxin_*.tsv   # UniProt venom neurotoxin database
├── interface_tensors.zip              # Precomputed interface tensors
├── H3i_am1bcc_fa_mm.params           # Rosetta parameter file
├── alk2_target.pdb                   # ALK2 reference target structure
├── long_chain_pdb_new.pdb            # 3FTx target structure
├── long_chain.fasta                  # Multiple 3FTx sequences from N. naja
├── long_chain_consensus.fa           # Consensus sequence
├── design_*.png                      # 3D structure visualizations
├── af2_score_analysis.png            # Comprehensive validation analysis
├── designs_excellent.csv             # Top-tier candidates with full metrics
└── per_tier_summary.csv              # Statistical summary across quality tiers
```

## Target Database

The project utilizes a comprehensive UniProt-derived database of venom neurotoxins, specifically targeting long neurotoxins from *Naja naja* (Indian cobra):

| UniProt ID | Protein Name | Length | Structure Available |
|------------|--------------|--------|-------------------|
| C0HM09 | Alpha-elapitoxin-Nn3a (Long neurotoxin 7) | 71 aa | - |
| P25668 | Long neurotoxin 1 (Toxin A) | 71 aa | AlphaFold |
| P25671 | Long neurotoxin 3 (Toxin C) | 71 aa | X-ray (PDB) + AlphaFold |
| P25672 | Long neurotoxin 4 (Toxin D) | 71 aa | AlphaFold |
| P25669 | Long neurotoxin 2 (Toxin B) | 71 aa | AlphaFold |
| P25673 | Long neurotoxin 5 (Toxin E) | 71 aa | AlphaFold |

These are three-finger toxins responsible for the neurotoxic effects of cobra envenomation, making them critical targets for antivenom development.

## Methodology

### 1. Target Preparation and Consensus Generation
- **Multiple Sequence Alignment**: Homologous 3FTx proteins aligned using Clustal Omega
- **Consensus Sequence Generation**: MSA-based consensus using configurable thresholds
- **Structure Preparation**: PDB renumbering and validation using PyRosetta
- **Database Mining**: Comprehensive UniProt search for venom neurotoxins

### 2. Systematic Interface Region Targeting
The pipeline systematically explores multiple interface regions on the 3FTx target:

| Interface Region | Target Residues | Hotspots | Strategy |
|------------------|----------------|----------|----------|
| **N-terminal Loop** | A28-33 | 6-residue interface | High-precision targeting |
| **Central Loop** | A30-34 | 5-residue interface | Core binding region |
| **C-terminal Region** | A36-39 | 4-residue interface | Alternative binding site |
| **Extended Interface** | A53-56 | 4-residue interface | Secondary contact region |

### 3. Scaffold-Guided Binder Design

#### Interface Tensor Generation
```python
# Example command for tensor generation
python make_interface_tensor.py \
    --input_pdb long_chain_pdb_new.pdb \
    --out_dir clamp_tensors_121_A28_33/len_121/ \
    --binderlen 121 \
    --target_adj A28-33 \
    --binder_ss E \
    --binder_ss_len 6
```

#### RFdiffusion with Scaffold Conditioning
```python
# Scaffold-guided design with hotspot targeting
python run_inference.py \
    inference.input_pdb=long_chain_pdb_new.pdb \
    contigmap.contigs=[A1-71/0 121] \
    ppi.hotspot_res=['A28','A29','A30','A31','A32','A33'] \
    scaffoldguided.scaffoldguided=True \
    scaffoldguided.scaffold_dir=clamp_tensors_121_A28_33/len_121/ \
    inference.num_designs=100
```

#### Parameter Space Exploration
- **Binder Lengths**: 115-125 amino acids (systematic scan)
- **Design Batch Sizes**: 7-100 designs per condition
- **Interface Specifications**: β-strand pairing with target loops
- **Scaffold Diversity**: Multiple tensor conformations per interface

### 4. Comprehensive Validation Pipeline
All designed binders undergo systematic evaluation using multiple metrics:

#### AlphaFold2 Confidence Metrics:
- **pLDDT scores**: Per-residue confidence (0-100)
  - `plddt_binder`: Confidence in designed binder structure
  - `plddt_target`: Confidence in target structure  
  - `plddt_total`: Overall complex confidence (weighted average)

#### Predicted Aligned Error (PAE):
- **pae_binder**: Internal binder domain accuracy (Å)
- **pae_target**: Internal target domain accuracy (Å)
- **pae_interaction**: Binder-target interface accuracy (Å) - **Primary filter**

#### Structural Validation:
- **binder_aligned_rmsd**: RMSD between AF2 prediction and RFdiffusion design (Å)
- **target_aligned_rmsd**: Target backbone preservation (Å)

#### Computational Metrics:
- **time**: AlphaFold2 prediction time per design (seconds)

### 5. Automated Analysis and Classification

#### Tier-Based Quality Assessment
Designs are automatically classified into four quality tiers based on interface confidence:

| Tier | PAE Interaction Range | Color Code | Interpretation |
|------|----------------------|------------|---------------|
| **Excellent** | 0 ≤ PAE < 10 Å | 🟢 Green | High confidence binding |
| **Good** | 10 ≤ PAE < 15 Å | 🟡 Yellow | Moderate confidence |
| **Marginal** | 15 ≤ PAE < 20 Å | 🟠 Orange | Low confidence |
| **Poor** | PAE ≥ 20 Å | 🔴 Red | Unlikely to bind |

#### Comprehensive Visualization
The analysis pipeline generates multi-panel diagnostic plots:
1. **Interface PAE Distribution**: Histogram with tier boundaries
2. **Overall Confidence**: pLDDT total distribution
3. **Per-Chain Analysis**: Binder vs target confidence scatter
4. **Design Fidelity**: Interface PAE vs structural preservation
5. **Cross-Validation**: RMSD analysis for design accuracy

#### Automated Filtering and Ranking
- **Per-tier CSV exports**: Separate files for each quality tier
- **Global ranking**: Multi-metric sorting (PAE → pLDDT → RMSD)
- **Statistical summaries**: Per-tier means, medians, and counts
- **Success rate calculation**: Percentage meeting quality thresholds

## Key Results

### Overall Performance Metrics
From systematic exploration across multiple interface regions and binder lengths:

| Metric | Best Achievement | Success Criteria |
|--------|-----------------|------------------|
| **Interface Confidence** | 6.2 Å PAE | < 10 Å (excellent) |
| **Binder Confidence** | 88.1 pLDDT | > 80 (high confidence) |
| **Target Preservation** | 91.7 pLDDT | > 85 (stable) |
| **Design Fidelity** | 0.6 Å RMSD | < 2 Å (accurate) |
| **Success Rate** | ~1.9% excellent | Variable by interface |

### Per-Interface Region Analysis

| Interface Region | Designs Generated | Excellent Candidates | Best PAE (Å) | Notes |
|------------------|-------------------|---------------------|--------------|--------|
| A28-33 (6-res) | 100 | TBD | TBD | High-precision targeting |
| A30-34 (5-res) | 100 | TBD | TBD | Core binding region |
| A36-39 (4-res) | 110 | TBD | TBD | Alternative site |
| A53-56 (4-res) | 110 | TBD | TBD | Extended interface |

### Quality Distribution (Example Dataset)
- **Excellent** (PAE < 10): 17 designs (~1.9%) - Ready for experimental validation
- **Good** (PAE 10-15): 8 designs (~0.9%) - Promising candidates
- **Marginal** (PAE 15-20): 16 designs (~1.8%) - Requires optimization
- **Poor** (PAE > 20): 855 designs (~95.4%) - Structural analysis for failure modes

### Computational Performance
- **Average prediction time**: ~17 seconds per design (AlphaFold2)
- **Tensor generation**: ~59 scaffolds per interface condition
- **Scaffold diversity**: Multiple conformational states per binder length
- **GPU requirements**: CUDA-enabled environment for RFdiffusion

## Installation and Usage

### Prerequisites
```bash
# Install RFdiffusion (required)
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion
# Follow installation instructions

# Install PyRosetta (required for PDB processing)
# Contact https://www.pyrosetta.org for academic license

# Install Python dependencies
pip install -r requirements.txt
```

### Quick Start

#### 1. Prepare Target Structure
```bash
# Renumber PDB if needed (edit line 15 for params file path)
python pdbinfo_reset.py your_target.pdb
```

#### 2. Generate Interface Tensors
```bash
# For 6-residue interface targeting A28-33
python generate_interface_tensors_121_A28_33.py
```

#### 3. Run Scaffold-Guided Design
```bash
# Generate 100 designs with interface conditioning
python run_rfd_designs_long_chain_121_A28_33.py
```

#### 4. Sequence Design with ProteinMPNN
```bash
# Switch to ProteinMPNN environment
conda deactivate
conda activate proteinmpnn_binder_design

# Design sequences for generated structures
python /home/sbl_lab1/protein_design/dl_binder_design/mpnn_fr/dl_interface_design.py \
  -pdbdir /home/sbl_lab1/protein_design/long_chain_design/RFDiffusion_designs_121_A28_33 \
  -outpdbdir /home/sbl_lab1/protein_design/long_chain_design/ProteinMPNN_designs_121_A28_33 \
  -relax_cycles 2 \
  -seqs_per_struct 8 \
  -temperature 0.1
```

#### 5. Validate with AlphaFold2
```bash
# Run AlphaFold2 on generated designs (requires separate setup)
# Generate confidence scores in Rosetta .sc format
```

#### 6. Comprehensive Analysis
```bash
# Generate multi-panel analysis plots and filtered results
python generate_scores_plots.py
```

### Advanced Usage

#### Custom Interface Targeting
```bash
# Complex multi-region interface
python make_interface_tensor.py \
    --input_pdb target.pdb \
    --out_dir adjacency_tensors/ \
    --target_adj A1-10:A50-60,A65-70:A90 \
    --binder_ss H,E,E \
    --binder_ss_len 8,8,8 \
    --binderlen 100
```

#### Multiple Sequence Analysis
```bash
# Generate consensus from aligned FASTA
python generate_consensus_sequence.py aligned_sequences.fa 0.7
```

## Technical Details

### Scaffold-Guided Design Innovation
This pipeline implements **interface tensor conditioning** for RFdiffusion, enabling:

- **Structural Constraints**: Pre-defined secondary structure adjacencies
- **Interface Geometry**: Precise control over binding site interactions
- **Scaffold Diversity**: Multiple conformational starting points
- **Systematic Exploration**: Parameter sweeps across design space

### Key Algorithm Features
- **Target-Agnostic**: Framework applicable beyond 3FTx targets
- **Multi-Interface**: Simultaneous targeting of multiple binding regions
- **Validation Integration**: Automated AlphaFold2 confidence assessment
- **Failure Analysis**: Diagnostic plots for design optimization

### File Formats and Data Flow
```
Target PDB → Interface Tensors → RFdiffusion → Design PDBs → 
AlphaFold2 → Confidence Scores → Analysis & Ranking → 
Filtered Candidates
```

## Applications and Impact

### Immediate Applications
1. **Computational Antivenom Development**
   - High-affinity binders for cobra neurotoxin neutralization
   - Broad-spectrum targeting across 3FTx variants
   - Rational design replacing animal-derived antibodies

2. **Therapeutic Discovery**
   - Lead compounds for anti-3FTx therapeutics
   - Template for targeting other venom components
   - Precision medicine for snakebite treatment

3. **Structural Biology Research**
   - Mechanistic insights into 3FTx-target interactions
   - Structure-function relationships in toxin binding
   - Validation of computational design methods

### Methodological Contributions
4. **Scaffold-Guided Design Framework**
   - General methodology for interface-constrained protein design
   - Integration candidate for official RFdiffusion repository
   - Benchmark for computational design validation

5. **Systematic Design Validation**
   - Multi-metric confidence assessment protocol
   - Automated failure mode analysis
   - Statistical framework for design success evaluation

## Future Directions

### Experimental Validation Pipeline
- **Protein Expression**: High-throughput production of top candidates
- **Binding Affinity**: Surface plasmon resonance (SPR) measurements
- **Neutralization Assays**: In vitro toxin neutralization studies
- **Structural Validation**: X-ray crystallography of best binders

### Computational Enhancements
- **Sequence Optimization**: ProteinMPNN-based sequence design
- **Stability Engineering**: Computational thermostability optimization
- **Affinity Maturation**: Iterative design-test cycles
- **Multi-Target Design**: Simultaneous binding to multiple 3FTx variants

### Broader Applications
- **Pan-Venom Targeting**: Extension to other snake species and toxin families
- **Antibody Engineering**: Application to therapeutic antibody design
- **Enzyme Inhibitors**: Targeting of other therapeutic protein targets
- **Vaccine Development**: Immunogen design for snakebite prevention

### Technology Transfer
- **RFdiffusion Integration**: Merger into official codebase
- **Open Source Release**: Community access to scaffold-guided methods
- **Training Resources**: Educational materials for computational design
- **Industrial Applications**: Technology transfer for therapeutics development

## Dependencies and Environment

### Core Requirements
- **RFdiffusion**: Latest version from RosettaCommons repository
- **PyRosetta**: Academic license required for PDB processing
- **AlphaFold2**: For structure prediction and confidence scoring
- **PyTorch**: CUDA-enabled for GPU acceleration

### Python Ecosystem
```bash
# Core scientific computing
numpy pandas matplotlib seaborn scikit-learn

# Bioinformatics
biopython

# Development environment  
jupyterlab ipykernel

# Platform-specific PyTorch installation required
```

### Hardware Recommendations
- **GPU**: CUDA-compatible for RFdiffusion and AlphaFold2
- **Memory**: 16+ GB RAM for large protein complexes
- **Storage**: SSD recommended for tensor I/O operations
- **CPU**: Multi-core for parallel design generation

## Citation and Acknowledgments

This work builds upon several key technologies:

- **RFdiffusion**: Watson et al., Nature (2023) - Protein structure generation
- **AlphaFold2**: Jumper et al., Nature (2021) - Structure prediction and confidence
- **PyRosetta**: Chaudhury et al., PLoS ONE (2010) - Protein structure manipulation
- **Three-Finger Toxin Research**: Extensive literature on cobra neurotoxin structure-function

*[Update with appropriate citations upon publication]*

## Contributing

This codebase is under active development for integration into the RFdiffusion repository. Contributions are welcome through:

- **Issues**: Bug reports and feature requests
- **Code Review**: Optimization and best practices
- **Validation**: Experimental testing of designed binders
- **Documentation**: Improved tutorials and examples

## Contact

For questions regarding this computational antivenom development pipeline or the scaffold-guided design methodology, please contact the research team.

---

**Note**: This represents a significant advancement in computational protein design for therapeutic applications, demonstrating the potential of AI-guided methods for addressing critical global health challenges like snakebite envenomation.
