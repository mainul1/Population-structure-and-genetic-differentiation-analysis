#!/usr/bin/env python3
"""
Population Genetics Analysis Script for Rhodnius prolixus
========================================================
Author: Mohammad Sarker
Date: November 2025
"""

import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import logging
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import subprocess
import os
import itertools
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PopulationGeneticsAnalyzer:
    """
    Population genetics analyzer for Rhodnius prolixus
    """
    
    def __init__(self, sequences_dir: str = "sequences", metadata_file: str = "Metadata.csv"):
        """Initialize analyzer"""
        self.sequences_dir = Path(sequences_dir)
        
        # Load metadata with encoding handling
        encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
        self.metadata = None
        
        for encoding in encodings:
            try:
                self.metadata = pd.read_csv(metadata_file, encoding=encoding)
                logger.info(f"Successfully loaded metadata with {encoding} encoding")
                break
            except UnicodeDecodeError:
                continue
        
        if self.metadata is None:
            raise ValueError("Could not read metadata CSV file with any encoding")
        
        # Initialize data structures
        self.sequences_by_gene = {}
        self.alignments = {}
        self.diversity_stats = {}
        self.population_groups = {}
        
        # Gene information
        self.nuclear_loci = {
            'TRNA': 597,
            'PJH': 601,
            'CISP': 522,
            'LSM': 569,
            'UPCA': 641,
            'UPMETAL': 506,
            '28S': 500
        }
        self.mitochondrial_loci = {
            'CYTB': 432
        }
        
        logger.info("Initialized analyzer")
    
    def load_sequences_by_gene(self) -> Dict[str, List]:
        """Load sequences by gene"""
        sequences_by_gene = {}
        
        for fasta_file in self.sequences_dir.glob("*_sequences.fasta"):
            gene = fasta_file.stem.replace("_sequences", "")
            sequences = list(SeqIO.parse(fasta_file, "fasta"))
            sequences_by_gene[gene] = sequences
            logger.info(f"Loaded {len(sequences)} sequences for {gene}")
        
        self.sequences_by_gene = sequences_by_gene
        return sequences_by_gene
    
    def perform_mafft_alignment(self, gene: str, sequences: List, output_dir: str = "alignments") -> Optional[str]:
        """
        Perform MAFFT alignment per locus
        """
        if len(sequences) < 2:
            logger.warning(f"Insufficient sequences for {gene} alignment")
            return None
        
        # Create alignment output directory
        align_dir = Path(output_dir)
        align_dir.mkdir(exist_ok=True)
        
        # Input and output files
        input_file = align_dir / f"{gene}_unaligned.fasta"
        output_file = align_dir / f"{gene}_aligned.fasta"
        
        # Write unaligned sequences
        SeqIO.write(sequences, input_file, "fasta")
        logger.info(f"Wrote {len(sequences)} unaligned {gene} sequences to {input_file}")
        
        # Try to run MAFFT (if installed locally)
        try:
            mafft_cmd = [
                "mafft",
                "--auto",
                "--inputorder",
                "--quiet",
                str(input_file)
            ]
            
            logger.info(f"Running MAFFT alignment for {gene}...")
            
            result = subprocess.run(
                mafft_cmd,
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if result.returncode == 0:
                with open(output_file, 'w') as f:
                    f.write(result.stdout)
                
                logger.info(f"MAFFT alignment completed for {gene} -> {output_file}")
                
                aligned_sequences = list(SeqIO.parse(output_file, "fasta"))
                self.alignments[gene] = aligned_sequences
                
                return str(output_file)
            else:
                logger.error(f"MAFFT failed for {gene}: {result.stderr}")
                return None
                
        except FileNotFoundError:
            logger.warning("MAFFT not found locally. Creating command for manual execution.")
            self._create_mafft_commands(gene, input_file, output_file)
            return None
        except subprocess.TimeoutExpired:
            logger.error(f"MAFFT alignment timed out for {gene}")
            return None
        except Exception as e:
            logger.error(f"MAFFT alignment failed for {gene}: {e}")
            return None
    
    def _create_mafft_commands(self, gene: str, input_file: Path, output_file: Path):
        """Create MAFFT commands for manual execution"""
        
        with open("mafft_commands.sh", 'a') as f:
            f.write(f"# MAFFT alignment for {gene}\n")
            f.write(f"mafft --auto --inputorder {input_file} > {output_file}\n")
            f.write(f"echo \"Completed alignment for {gene}\"\n\n")
        
        logger.info(f"MAFFT command for {gene} added to mafft_commands.sh")
    
    def align_all_sequences(self):
        """
        Perform MAFFT alignment for all gene loci
        """
        logger.info("Starting MAFFT alignment for all loci...")
        
        aligned_files = {}
        
        # Create master command file
        with open("mafft_commands.sh", 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# MAFFT Alignment Commands\n\n")
        
        for gene, sequences in self.sequences_by_gene.items():
            aligned_file = self.perform_mafft_alignment(gene, sequences)
            if aligned_file:
                aligned_files[gene] = aligned_file
        
        if not aligned_files and self.sequences_by_gene:
            logger.warning("MAFFT not available locally. Install MAFFT or run: bash mafft_commands.sh")
        
        return aligned_files
    
    def create_concatenated_alignments(self):
        """
        Create concatenated alignments:
        1. Nuclear alignment (7 loci, 3,936 bp)
        2. Full alignment (nuclear + mitochondrial, 4,368 bp)
        """
        logger.info("Creating concatenated alignments...")
        
        if not self.alignments:
            logger.error("No aligned sequences found. Run MAFFT alignment first.")
            return None, None
        
        nuclear_order = ['TRNA', 'PJH', 'CISP', 'LSM', 'UPCA', 'UPMETAL', '28S']
        mitochondrial_order = ['CYTB']
        
        nuclear_alignment = self._concatenate_loci(nuclear_order, "nuclear_alignment.fasta")
        full_order = nuclear_order + mitochondrial_order
        full_alignment = self._concatenate_loci(full_order, "full_alignment.fasta")
        
        return nuclear_alignment, full_alignment
    
    def _concatenate_loci(self, loci_order: List[str], output_file: str) -> Optional[str]:
        """Concatenate aligned loci in specified order"""
        
        available_loci = []
        for locus in loci_order:
            if locus in self.alignments:
                available_loci.append(locus)
            else:
                logger.warning(f"Aligned sequences for {locus} not found")
        
        if not available_loci:
            logger.error("No aligned loci available for concatenation")
            return None
        
        # Get sample IDs that are present in all loci
        sample_sets = []
        for locus in available_loci:
            locus_samples = {seq.id.split('_')[0] for seq in self.alignments[locus]}
            sample_sets.append(locus_samples)
        
        common_samples = set.intersection(*sample_sets) if sample_sets else set()
        
        if not common_samples:
            logger.error("No samples common to all loci for concatenation")
            return None
        
        logger.info(f"Concatenating {len(available_loci)} loci for {len(common_samples)} samples")
        
        # Create concatenated sequences
        concatenated_records = []
        
        for sample_id in sorted(common_samples):
            concatenated_seq = ""
            
            for locus in available_loci:
                locus_seq = None
                for seq in self.alignments[locus]:
                    if seq.id.startswith(sample_id):
                        locus_seq = str(seq.seq)
                        break
                
                if locus_seq:
                    concatenated_seq += locus_seq
                else:
                    expected_length = self.nuclear_loci.get(locus, self.mitochondrial_loci.get(locus, 500))
                    concatenated_seq += "N" * expected_length
            
            record = SeqRecord(
                Seq(concatenated_seq),
                id=sample_id,
                description=f"Concatenated sequence from {len(available_loci)} loci"
            )
            concatenated_records.append(record)
        
        # Save concatenated alignment
        output_path = Path("alignments") / output_file
        SeqIO.write(concatenated_records, output_path, "fasta")
        
        total_length = len(concatenated_records[0].seq) if concatenated_records else 0
        logger.info(f"Saved concatenated alignment: {output_path} ({total_length} bp, {len(concatenated_records)} samples)")
        
        return str(output_path)
    
    def create_population_groups(self) -> Dict[str, List[str]]:
        """
        Create population groups for analysis:
        - By Country (primary structure)
        - By Country-State (secondary structure)
        - By Geographic regions (tertiary structure)
        """
        groups = {}
        
        # Level 1: By Country
        for country in self.metadata['Country'].unique():
            if pd.notna(country):
                samples = self.metadata[self.metadata['Country'] == country]['ID'].tolist()
                groups[f"Country_{country}"] = samples
        
        # Level 2: By Country-State
        for _, row in self.metadata.iterrows():
            if pd.notna(row['Country']) and pd.notna(row['State']):
                group_name = f"{row['Country']}_{row['State']}"
                if group_name not in groups:
                    groups[group_name] = []
                groups[group_name].append(row['ID'])
        
        # Level 3: Geographic regions
        colombia_north = []  # Cesar, Magdalena
        colombia_central = []  # Casanare, Meta
        colombia_west = []   # Antioquia
        
        for _, row in self.metadata.iterrows():
            if row['Country'] == 'Colombia':
                if row['State'] in ['Cesar', 'Magdalena']:
                    colombia_north.append(row['ID'])
                elif row['State'] in ['Casanare', 'Meta']:
                    colombia_central.append(row['ID'])
                elif row['State'] in ['Antioquia']:
                    colombia_west.append(row['ID'])
        
        if colombia_north: groups['Colombia_North'] = colombia_north
        if colombia_central: groups['Colombia_Central'] = colombia_central  
        if colombia_west: groups['Colombia_West'] = colombia_west
        
        # Filter groups with at least 3 samples
        filtered_groups = {name: samples for name, samples in groups.items() 
                          if len(samples) >= 3}
        
        logger.info(f"Created {len(filtered_groups)} population groups")
        for name, samples in filtered_groups.items():
            logger.info(f"  {name}: {len(samples)} samples")
        
        self.population_groups = filtered_groups
        return filtered_groups
    
    def calculate_diversity_statistics(self, sequences: List[str], gene: str) -> Dict:
        """
        Calculate genetic diversity measures:
        - haplotype diversity (h)
        - number of segregating sites (S) 
        - population substitution rate (θ)
        - nucleotide diversity (π)
        - Tajima's D
        """
        n = len(sequences)
        if n < 2:
            return {'error': 'Insufficient sequences'}
        
        stats = {}
        stats['n_sequences'] = n
        stats['sequence_length'] = len(sequences[0]) if sequences else 0
        
        # 1. Segregating Sites (S)
        segregating_sites = 0
        polymorphic_sites = []
        
        for pos in range(len(sequences[0])):
            column = [seq[pos] for seq in sequences if pos < len(seq)]
            unique_bases = set(column)
            if len(unique_bases) > 1:
                segregating_sites += 1
                polymorphic_sites.append(pos)
        
        stats['segregating_sites'] = segregating_sites
        
        # 2. Population substitution rate (θ = S/a1)
        if n > 1:
            a1 = sum(1/i for i in range(1, n))
            theta = segregating_sites / a1 if a1 > 0 else 0
        else:
            theta = 0
        stats['theta'] = theta
        
        # 3. Nucleotide diversity (π)
        total_differences = 0
        comparisons = 0
        
        for i in range(n):
            for j in range(i + 1, n):
                differences = sum(1 for k in range(min(len(sequences[i]), len(sequences[j]))) 
                                if sequences[i][k] != sequences[j][k])
                valid_length = min(len(sequences[i]), len(sequences[j]))
                if valid_length > 0:
                    total_differences += differences / valid_length
                    comparisons += 1
        
        pi = total_differences / comparisons if comparisons > 0 else 0
        stats['nucleotide_diversity'] = pi
        
        # 4. Haplotype diversity (h)
        haplotype_counts = Counter(sequences)
        total = sum(haplotype_counts.values())
        h = 1 - sum((count/total)**2 for count in haplotype_counts.values())
        stats['haplotype_diversity'] = h
        stats['num_haplotypes'] = len(haplotype_counts)
        
        # 5. Tajima's D test for neutrality
        if segregating_sites > 0 and n > 1:
            a1 = sum(1/i for i in range(1, n))
            a2 = sum(1/(i**2) for i in range(1, n))
            
            b1 = (n + 1) / (3 * (n - 1))
            b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
            
            c1 = b1 - (1/a1)
            c2 = b2 - ((n + 2)/(a1 * n)) + (a2/(a1**2))
            
            e1 = c1 / a1
            e2 = c2 / (a1**2 + a2)
            
            var_D = e1 * segregating_sites + e2 * segregating_sites * (segregating_sites - 1)
            
            if var_D > 0:
                tajima_d = (pi - theta) / np.sqrt(var_D)
                stats['tajima_d'] = tajima_d
            else:
                stats['tajima_d'] = 0
        else:
            stats['tajima_d'] = 0
        
        logger.info(f"{gene} statistics: S={segregating_sites}, π={pi:.4f}, h={h:.4f}, D={stats['tajima_d']:.4f}")
        
        return stats
    
    def calculate_fst_hudson(self, pop1_sequences: List[str], pop2_sequences: List[str]) -> float:
        """
        Calculate FST using Hudson's method
        """
        n1, n2 = len(pop1_sequences), len(pop2_sequences)
        
        if n1 < 2 or n2 < 2:
            return 0.0
        
        # Calculate π within populations
        pi1 = self._calculate_pi_within(pop1_sequences)
        pi2 = self._calculate_pi_within(pop2_sequences)
        
        # Calculate π between populations
        pi_between = self._calculate_pi_between(pop1_sequences, pop2_sequences)
        
        # Average within-population π weighted by sample size
        pi_within = ((n1 - 1) * pi1 + (n2 - 1) * pi2) / (n1 + n2 - 2)
        
        # FST = (π_between - π_within) / π_between
        if pi_between > 0:
            fst = (pi_between - pi_within) / pi_between
            return max(0, fst)
        else:
            return 0.0
    
    def _calculate_pi_within(self, sequences: List[str]) -> float:
        """Calculate nucleotide diversity within population"""
        n = len(sequences)
        if n < 2:
            return 0.0
        
        total_diff = 0
        comparisons = 0
        
        for i in range(n):
            for j in range(i + 1, n):
                diff = sum(1 for k in range(min(len(sequences[i]), len(sequences[j]))) 
                          if sequences[i][k] != sequences[j][k])
                length = min(len(sequences[i]), len(sequences[j]))
                if length > 0:
                    total_diff += diff / length
                    comparisons += 1
        
        return total_diff / comparisons if comparisons > 0 else 0.0
    
    def _calculate_pi_between(self, pop1: List[str], pop2: List[str]) -> float:
        """Calculate nucleotide diversity between populations"""
        total_diff = 0
        comparisons = 0
        
        for seq1 in pop1:
            for seq2 in pop2:
                diff = sum(1 for k in range(min(len(seq1), len(seq2))) 
                          if seq1[k] != seq2[k])
                length = min(len(seq1), len(seq2))
                if length > 0:
                    total_diff += diff / length
                    comparisons += 1
        
        return total_diff / comparisons if comparisons > 0 else 0.0
    
    def calculate_pairwise_fst_matrix(self, sequences_by_pop: Dict[str, List[str]]) -> pd.DataFrame:
        """
        Calculate pairwise FST matrix using Hudson's method
        """
        populations = list(sequences_by_pop.keys())
        n_pops = len(populations)
        
        fst_matrix = np.zeros((n_pops, n_pops))
        
        for i in range(n_pops):
            for j in range(n_pops):
                if i == j:
                    fst_matrix[i, j] = 0.0
                else:
                    fst = self.calculate_fst_hudson(
                        sequences_by_pop[populations[i]], 
                        sequences_by_pop[populations[j]]
                    )
                    fst_matrix[i, j] = fst
        
        fst_df = pd.DataFrame(fst_matrix, index=populations, columns=populations)
        
        logger.info("Calculated pairwise FST matrix using Hudson's method")
        return fst_df
    
    def generate_results_table(self, output_file: str = "diversity_results.csv"):
        """
        Generate results table with diversity statistics
        """
        if not self.diversity_stats:
            logger.warning("No diversity statistics calculated yet. Run analysis first.")
            return
        
        results_data = []
        
        for gene, stats in self.diversity_stats.items():
            locus_type = "Nuclear" if gene in self.nuclear_loci else "Mitochondrial"
            expected_length = self.nuclear_loci.get(gene, self.mitochondrial_loci.get(gene, 0))
            
            result_row = {
                'Locus': gene,
                'Type': locus_type,
                'Length_bp': expected_length,
                'N_sequences': stats.get('n_sequences', 0),
                'N_haplotypes': stats.get('num_haplotypes', 0),
                'Segregating_sites': stats.get('segregating_sites', 0),
                'Theta': stats.get('theta', 0),
                'Pi': stats.get('nucleotide_diversity', 0),
                'Haplotype_diversity': stats.get('haplotype_diversity', 0),
                'Tajima_D': stats.get('tajima_d', 0)
            }
            results_data.append(result_row)
        
        # Create summary row for nuclear loci
        nuclear_stats = [stats for gene, stats in self.diversity_stats.items() 
                        if gene in self.nuclear_loci]
        
        if nuclear_stats:
            nuclear_summary = {
                'Locus': 'Nuclear_Combined',
                'Type': 'Nuclear',
                'Length_bp': sum(self.nuclear_loci.values()),
                'N_sequences': np.mean([s['n_sequences'] for s in nuclear_stats]),
                'N_haplotypes': np.sum([s['num_haplotypes'] for s in nuclear_stats]),
                'Segregating_sites': np.sum([s['segregating_sites'] for s in nuclear_stats]),
                'Theta': np.mean([s['theta'] for s in nuclear_stats]),
                'Pi': np.mean([s['nucleotide_diversity'] for s in nuclear_stats]),
                'Haplotype_diversity': np.mean([s['haplotype_diversity'] for s in nuclear_stats]),
                'Tajima_D': np.mean([s['tajima_d'] for s in nuclear_stats])
            }
            results_data.append(nuclear_summary)
        
        # Save results
        results_df = pd.DataFrame(results_data)
        results_df.to_csv(output_file, index=False)
        
        logger.info(f"Results table saved to {output_file}")
        
        # Print formatted table
        print("\n" + "="*80)
        print("GENETIC DIVERSITY RESULTS")
        print("="*80)
        print(f"{'Locus':<15} {'Type':<12} {'Length':<8} {'N':<5} {'H':<5} {'S':<5} {'θ':<8} {'π':<8} {'h':<8} {'D':<8}")
        print("-" * 80)
        
        for _, row in results_df.iterrows():
            print(f"{row['Locus']:<15} {row['Type']:<12} {row['Length_bp']:<8.0f} "
                  f"{row['N_sequences']:<5.0f} {row['N_haplotypes']:<5.0f} {row['Segregating_sites']:<5.0f} "
                  f"{row['Theta']:<8.4f} {row['Pi']:<8.4f} {row['Haplotype_diversity']:<8.4f} {row['Tajima_D']:<8.3f}")
        
        return results_df
    
    def run_complete_analysis(self):
        """
        Run complete population genetics analysis
        """
        logger.info("Starting complete analysis pipeline...")
        
        # Step 1: Load sequences
        self.load_sequences_by_gene()
        
        # Step 2: MAFFT Alignment
        logger.info("Step 2: Performing MAFFT alignment per locus...")
        aligned_files = self.align_all_sequences()
        
        # Step 3: Create concatenated alignments
        logger.info("Step 3: Creating concatenated alignments...")
        nuclear_alignment, full_alignment = self.create_concatenated_alignments()
        
        # Step 4: Create population groups
        self.create_population_groups()
        
        # Step 5: Calculate diversity statistics
        logger.info("Step 5: Calculating diversity statistics...")
        
        for gene, sequences in self.sequences_by_gene.items():
            if len(sequences) > 1:
                if gene in self.alignments:
                    seq_strings = [str(seq.seq) for seq in self.alignments[gene]]
                    logger.info(f"Using aligned sequences for {gene} statistics")
                else:
                    seq_strings = [str(seq.seq) for seq in sequences]
                    logger.info(f"Using raw sequences for {gene} statistics")
                
                self.diversity_stats[gene] = self.calculate_diversity_statistics(seq_strings, gene)
        
        # Step 6: Generate results table
        results_df = self.generate_results_table()
        
        logger.info("Complete analysis finished!")
        
        return {
            'diversity_stats': self.diversity_stats,
            'population_groups': self.population_groups,
            'results_table': results_df,
            'aligned_files': aligned_files,
            'nuclear_alignment': nuclear_alignment,
            'full_alignment': full_alignment
        }

def main():
    """
    Main function for population genetics analysis
    """
    print("RHODNIUS PROLIXUS POPULATION GENETICS ANALYSIS")
    print("=" * 60)
    
    try:
        # Check if sequences exist
        sequences_dir = Path("sequences")
        if not sequences_dir.exists():
            logger.error("Sequences directory not found. Please run sequence extraction first.")
            return
        
        # Initialize analyzer
        analyzer = PopulationGeneticsAnalyzer()
        
        # Run complete analysis
        results = analyzer.run_complete_analysis()
        
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE")
        print("="*60)
        
        print(f"\nAnalysis summary:")
        print(f"- Genes analyzed: {len(results['diversity_stats'])}")
        print(f"- Population groups: {len(results['population_groups'])}")
        print(f"- MAFFT alignments: {len(results.get('aligned_files', {}))}")
        
        if results.get('nuclear_alignment'):
            print(f"- Nuclear concatenated alignment: {results['nuclear_alignment']}")
        if results.get('full_alignment'):
            print(f"- Full concatenated alignment: {results['full_alignment']}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise

if __name__ == "__main__":
    main()
