# _Rhodnius prolixus_ Population Genetics Analysis Pipeline
## Evolutionary Genetics Assignment

### Overview
This pipeline extracts and analyzes DNA sequences for population structure and genetic differentiation analysis of *Rhodnius prolixus* populations across Colombia, Brazil, and Venezuela. The analysis follows methods from the reference paper (Nature Scientific Reports, 2025).

### Files in this Package
- `sequence_extraction_script.py` - Main script to extract sequences from NCBI GenBank
- `population_genetics_analyzer.py` - Basic population genetics analysis script  
- `requirements.txt` - Python package dependencies
- `Metadata.csv` - Sample metadata with accession numbers
- `README.md` - This instruction file

### Gene Loci Analyzed
Based on your metadata, the following 8 loci will be extracted:
1. **28S** - Nuclear ribosomal RNA gene
2. **CISP** - Nuclear protein-coding gene
3. **TRNA** - Transfer RNA gene
4. **CYTB** - Mitochondrial cytochrome b gene
5. **UPMETAL** - Nuclear protein-coding gene
6. **UPCA** - Nuclear protein-coding gene  
7. **PJH** - Nuclear protein-coding gene
8. **LSM** - Nuclear protein-coding gene

### Installation & Setup

#### 1. Install Required Software
```powershell
# Install Python packages
pip install -r requirements.txt
```

#### 2. Configure NCBI Access
Edit the `sequence_extraction_script.py` file and replace `"your.email@example.com"` with your actual email address (required by NCBI):

```python
extractor = SequenceExtractor(email="your.actual.email@domain.com")
```

### Running the Analysis

#### Step 1: Extract Sequences from NCBI GenBank
```powershell
python sequence_extraction_script.py
```

**What this does:**
- Reads accession numbers from `Metadata.csv`
- Downloads sequences from NCBI GenBank for all 8 gene loci
- Organizes sequences by gene type
- Creates mapping files and summary reports

**Expected Output:**
- `sequences/` directory with FASTA files for each gene
- `sample_sequence_mapping.csv` - Shows which sequences were retrieved per sample
- `extraction_summary.txt` - Summary statistics
- `sequence_extraction.log` - Detailed log file

#### Step 2: Basic Population Analysis
```powershell
python population_analysis_script.py
```

**What this does:**
- Calculates basic genetic diversity statistics
- Groups samples by country and region
- Generates summary visualizations
- Prepares data for advanced analyses

**Expected Output:**
- `population_analysis_summary.txt` - Comprehensive analysis summary
- `population_analysis_overview.png` - Visualization plots

### Key Results Expected

Based on the methods in your reference paper, you should obtain:

#### 1. Sequence Statistics (per gene)
| Metric | Description |
|--------|-------------|
| **Number of sequences** | Total sequences retrieved per locus |
| **Sequence length** | Base pairs per aligned locus |
| **Segregating sites (S)** | Variable nucleotide positions |
| **Nucleotide diversity (π)** | Average pairwise differences |
| **Haplotype diversity (h)** | Genetic variation measure |

#### 2. Population Structure Results
- **Sample distribution**: Colombia (~94 samples), Venezuela (~29 samples), Brazil (~6 samples)
- **Geographic clustering**: Samples grouped by country/state
- **Population differentiation**: FST values between populations
- **Gene flow patterns**: Connectivity between regions

#### 3. Files for Advanced Analysis
The extracted sequences can be used for:
- **MAFFT alignment** (as mentioned in your methods)
- **STRUCTURE/ADMIXTURE** clustering analysis
- **IQ-Tree/PhyML** phylogenetic reconstruction
- **DnaSP** neutrality tests (Tajima's D, Fu & Li's D)

### Next Steps for Complete Analysis

After running these scripts, you'll need specialized software for advanced analyses:

#### 1. Multiple Sequence Alignment
```bash
# Using MAFFT (online or local installation)
mafft --auto sequences/28S_sequences.fasta > 28S_aligned.fasta
```

#### 2. Population Structure (STRUCTURE software)
- Input: Concatenated nuclear loci alignment
- Parameters: K=1-20, 100,000 MCMC iterations
- Use STRUCTURE HARVESTER to determine optimal K

#### 3. Phylogenetic Analysis (IQ-Tree)
```bash
iqtree -s concatenated_alignment.fasta -bb 1000 -alrt 1000
```

#### 4. Population Genetics (DnaSP)
- Calculate FST, neutrality tests
- Estimate divergence times
- Analyze gene flow patterns

### Expected Key Results Summary

Your analysis should reveal:

1. **Population Structure**:
   - Primary clustering by country (Colombia, Venezuela, Brazil)
   - Secondary structure by geographic regions within countries
   - Possible isolation-by-distance patterns

2. **Genetic Differentiation**:
   - Highest FST between Brazil and other countries
   - Moderate differentiation between Colombian regions
   - Gene flow corridors along major river systems

3. **Diversity Patterns**:
   - Higher diversity in Colombian populations
   - Reduced diversity in peripheral populations
   - Different patterns between nuclear and mitochondrial loci

4. **Evolutionary Implications**:
   - Recent population expansion
   - Historical bottlenecks in some regions
   - Adaptation to different ecotopes (sylvatic/domestic)

### Troubleshooting

#### Common Issues:
1. **NCBI connection errors**: Add delays between requests, check internet connection
2. **Missing sequences**: Some accession numbers may be invalid or restricted
3. **Memory issues**: Process genes separately if dataset is too large

#### Error Logs:
Check `sequence_extraction.log` for detailed error messages and failed accessions.

### Citation & Methods Statement

**Methods Summary for Paper:**
"DNA sequences were retrieved from NCBI GenBank using accession numbers for 8 nuclear and mitochondrial loci. Sequences were aligned using MAFFT v7 and analyzed using custom Python scripts based on BioPython. Population structure was assessed using [specify methods used], genetic differentiation calculated via FST estimation, and phylogenetic relationships reconstructed using maximum likelihood in IQ-Tree."

### Contact & Support
For questions about this pipeline or evolutionary genetics analysis, consult:
- BioPython documentation: https://biopython.org/
- NCBI E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- Population genetics textbooks for theoretical background

---
*Analysis pipeline for Rhodnius prolixus population genetics - November 2025*



