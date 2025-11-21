#!/usr/bin/env python3
"""
Sequence Extraction Script for Rhodnius prolixus Population Genetics Study
=========================================================================

This script extracts DNA sequences from NCBI GenBank using accession numbers
from metadata CSV file for evolutionary genetics analysis.

Author: Evolutionary Genetics Analysis
Date: November 2025
Assignment: Population structure and genetic differentiation analysis
"""

import pandas as pd
import requests
import time
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Optional
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('sequence_extraction.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class SequenceExtractor:
    """
    Class for extracting sequences from NCBI GenBank using accession numbers
    """
    
    def __init__(self, email: str = "your.email@example.com"):
        """
        Initialize the SequenceExtractor
        
        Args:
            email: Your email address for NCBI Entrez queries
        """
        self.email = email
        Entrez.email = email
        self.gene_columns = {
            '28S': ['28s', '28s.1'],
            'CISP': ['CISP'],
            'TRNA': ['TRNA', 'TRNA.1'],
            'CYTB': ['CYTB'],
            'UPMETAL': ['UPMETAL', 'UPMETAL.1'],
            'UPCA': ['UPCA', 'UPCA.1'],
            'PJH': ['PJH', 'PJH.1'],
            'LSM': ['LSM', 'LSM.1']
        }
        self.sequences = {}
        self.failed_accessions = []
        
    def load_metadata(self, csv_file: str) -> pd.DataFrame:
        """
        Load metadata from CSV file
        
        Args:
            csv_file: Path to metadata CSV file
            
        Returns:
            DataFrame with metadata
        """
        try:
            # Try different encodings
            encodings = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
            df = None
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(csv_file, encoding=encoding)
                    logger.info(f"Successfully loaded metadata with {encoding} encoding")
                    break
                except UnicodeDecodeError:
                    continue
            
            if df is None:
                raise ValueError("Could not read CSV file with any encoding")
            
            logger.info(f"Loaded metadata with {len(df)} samples")
            logger.info(f"Columns: {list(df.columns)}")
            return df
        except Exception as e:
            logger.error(f"Error loading metadata: {e}")
            raise
    
    def extract_accessions(self, df: pd.DataFrame) -> Dict[str, List[str]]:
        """
        Extract all unique accession numbers by gene
        
        Args:
            df: Metadata DataFrame
            
        Returns:
            Dictionary with gene names as keys and lists of accessions as values
        """
        accessions_by_gene = {}
        
        for gene, columns in self.gene_columns.items():
            accessions = set()
            for col in columns:
                if col in df.columns:
                    # Remove NaN values and empty strings
                    valid_accessions = df[col].dropna().astype(str)
                    valid_accessions = valid_accessions[valid_accessions != '']
                    accessions.update(valid_accessions.tolist())
            
            accessions_by_gene[gene] = list(accessions)
            logger.info(f"{gene}: Found {len(accessions)} unique accessions")
        
        return accessions_by_gene
    
    def fetch_sequence_batch(self, accessions: List[str], retmax: int = 50) -> Dict[str, SeqRecord]:
        """
        Fetch sequences in batches from NCBI
        
        Args:
            accessions: List of accession numbers
            retmax: Maximum number of sequences per batch
            
        Returns:
            Dictionary with accession as key and SeqRecord as value
        """
        sequences = {}
        
        # Process in batches
        for i in range(0, len(accessions), retmax):
            batch = accessions[i:i + retmax]
            logger.info(f"Fetching batch {i//retmax + 1}: {len(batch)} sequences")
            
            try:
                # Search for sequences
                search_handle = Entrez.esearch(
                    db="nucleotide",
                    term=" OR ".join(batch),
                    retmax=retmax
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    # Fetch sequences
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=search_results["IdList"],
                        rettype="fasta",
                        retmode="text"
                    )
                    
                    # Parse sequences
                    for record in SeqIO.parse(fetch_handle, "fasta"):
                        # Extract accession from description
                        accession = self.extract_accession_from_id(record.id, record.description)
                        if accession:
                            sequences[accession] = record
                            logger.debug(f"Successfully fetched {accession}")
                    
                    fetch_handle.close()
                
                # Be nice to NCBI servers
                time.sleep(0.5)
                
            except Exception as e:
                logger.error(f"Error fetching batch: {e}")
                self.failed_accessions.extend(batch)
                time.sleep(2)  # Wait longer on error
        
        return sequences
    
    def extract_accession_from_id(self, seq_id: str, description: str) -> Optional[str]:
        """
        Extract accession number from sequence ID or description
        
        Args:
            seq_id: Sequence ID
            description: Sequence description
            
        Returns:
            Accession number if found, None otherwise
        """
        # Try to extract from ID first
        accession_match = re.search(r'([A-Z]{1,2}\d{6,8})', seq_id)
        if accession_match:
            return accession_match.group(1)
        
        # Try to extract from description
        accession_match = re.search(r'([A-Z]{1,2}\d{6,8})', description)
        if accession_match:
            return accession_match.group(1)
        
        return None
    
    def save_sequences_by_gene(self, sequences: Dict[str, SeqRecord], 
                               accessions_by_gene: Dict[str, List[str]], 
                               output_dir: str = "sequences"):
        """
        Save sequences organized by gene
        
        Args:
            sequences: Dictionary of all sequences
            accessions_by_gene: Accessions organized by gene
            output_dir: Output directory for sequence files
        """
        Path(output_dir).mkdir(exist_ok=True)
        
        for gene, accessions in accessions_by_gene.items():
            gene_sequences = []
            
            for acc in accessions:
                if acc in sequences:
                    seq_record = sequences[acc]
                    # Rename sequence with sample info
                    seq_record.id = f"{acc}_{gene}"
                    seq_record.description = f"{gene} gene sequence from {acc}"
                    gene_sequences.append(seq_record)
                else:
                    logger.warning(f"Sequence not found for {acc} ({gene})")
            
            if gene_sequences:
                output_file = Path(output_dir) / f"{gene}_sequences.fasta"
                SeqIO.write(gene_sequences, output_file, "fasta")
                logger.info(f"Saved {len(gene_sequences)} {gene} sequences to {output_file}")
            else:
                logger.warning(f"No sequences found for {gene}")
    
    def create_sample_sequence_map(self, df: pd.DataFrame, sequences: Dict[str, SeqRecord], 
                                   output_file: str = "sample_sequence_mapping.csv"):
        """
        Create a mapping file showing which sequences were successfully retrieved for each sample
        
        Args:
            df: Metadata DataFrame
            sequences: Dictionary of retrieved sequences
            output_file: Output CSV file name
        """
        mapping_data = []
        
        for _, row in df.iterrows():
            sample_id = row['ID']
            country = row['Country']
            species = row['Species name']
            
            sequence_status = {}
            
            for gene, columns in self.gene_columns.items():
                found_sequences = []
                for col in columns:
                    if col in row and pd.notna(row[col]) and row[col] != '':
                        accession = str(row[col])
                        if accession in sequences:
                            found_sequences.append(accession)
                
                sequence_status[f"{gene}_accessions"] = ";".join(found_sequences) if found_sequences else "Not found"
                sequence_status[f"{gene}_count"] = len(found_sequences)
            
            mapping_row = {
                'Sample_ID': sample_id,
                'Species': species,
                'Country': country,
                'State': row.get('State', ''),
                'Town': row.get('Town', ''),
                **sequence_status
            }
            mapping_data.append(mapping_row)
        
        mapping_df = pd.DataFrame(mapping_data)
        mapping_df.to_csv(output_file, index=False)
        logger.info(f"Sample-sequence mapping saved to {output_file}")
        return mapping_df
    
    def generate_summary_report(self, accessions_by_gene: Dict[str, List[str]], 
                                sequences: Dict[str, SeqRecord], 
                                output_file: str = "extraction_summary.txt"):
        """
        Generate a summary report of the extraction process
        
        Args:
            accessions_by_gene: Accessions organized by gene
            sequences: Dictionary of retrieved sequences
            output_file: Output text file name
        """
        with open(output_file, 'w') as f:
            f.write("SEQUENCE EXTRACTION SUMMARY REPORT\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("GENE-WISE EXTRACTION RESULTS:\n")
            f.write("-" * 30 + "\n")
            
            total_requested = 0
            total_retrieved = 0
            
            for gene, accessions in accessions_by_gene.items():
                requested = len(accessions)
                retrieved = sum(1 for acc in accessions if acc in sequences)
                success_rate = (retrieved / requested * 100) if requested > 0 else 0
                
                f.write(f"{gene}:\n")
                f.write(f"  Requested: {requested}\n")
                f.write(f"  Retrieved: {retrieved}\n")
                f.write(f"  Success rate: {success_rate:.1f}%\n\n")
                
                total_requested += requested
                total_retrieved += retrieved
            
            overall_success = (total_retrieved / total_requested * 100) if total_requested > 0 else 0
            
            f.write("OVERALL STATISTICS:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total accessions requested: {total_requested}\n")
            f.write(f"Total sequences retrieved: {total_retrieved}\n")
            f.write(f"Overall success rate: {overall_success:.1f}%\n\n")
            
            if self.failed_accessions:
                f.write("FAILED ACCESSIONS:\n")
                f.write("-" * 18 + "\n")
                for acc in self.failed_accessions:
                    f.write(f"{acc}\n")
        
        logger.info(f"Summary report saved to {output_file}")

def main():
    """
    Main function to run the sequence extraction pipeline
    """
    logger.info("Starting Rhodnius prolixus sequence extraction pipeline")
    
    # Initialize extractor
    extractor = SequenceExtractor(email="student.genetics@example.com")  # Using placeholder email
    
    try:
        # Load metadata
        df = extractor.load_metadata("Metadata.csv")
        
        # Extract accessions by gene
        logger.info("Extracting accession numbers by gene...")
        accessions_by_gene = extractor.extract_accessions(df)
        
        # Flatten all accessions for batch retrieval
        all_accessions = []
        for gene_accessions in accessions_by_gene.values():
            all_accessions.extend(gene_accessions)
        
        unique_accessions = list(set(all_accessions))
        logger.info(f"Total unique accessions to fetch: {len(unique_accessions)}")
        
        # Fetch sequences
        logger.info("Fetching sequences from NCBI GenBank...")
        sequences = extractor.fetch_sequence_batch(unique_accessions)
        
        logger.info(f"Successfully retrieved {len(sequences)} sequences")
        
        # Save sequences by gene
        logger.info("Organizing and saving sequences by gene...")
        extractor.save_sequences_by_gene(sequences, accessions_by_gene)
        
        # Create sample-sequence mapping
        logger.info("Creating sample-sequence mapping...")
        mapping_df = extractor.create_sample_sequence_map(df, sequences)
        
        # Generate summary report
        logger.info("Generating summary report...")
        extractor.generate_summary_report(accessions_by_gene, sequences)
        
        logger.info("Sequence extraction pipeline completed successfully!")
        
        # Print quick summary
        print("\n" + "="*50)
        print("EXTRACTION COMPLETE")
        print("="*50)
        print(f"Total samples in metadata: {len(df)}")
        print(f"Unique accessions found: {len(unique_accessions)}")
        print(f"Sequences successfully retrieved: {len(sequences)}")
        print(f"Success rate: {len(sequences)/len(unique_accessions)*100:.1f}%")
        print("\nOutput files created:")
        print("- sequences/[GENE]_sequences.fasta (one per gene)")
        print("- sample_sequence_mapping.csv")
        print("- extraction_summary.txt")
        print("- sequence_extraction.log")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        raise

if __name__ == "__main__":
    main()
