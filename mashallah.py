import os
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt, gc_fraction
import primer3
import nupack

# Function to read sequences from a FASTA file
def read_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences  # Returns a list of sequences (if multiple entries)

# Function to generate candidate primers
def generate_primers(sequence, primer_length=30):
    primers = []
    for i in range(len(sequence) - primer_length + 1):
        primer_seq = sequence[i:i + primer_length]
        primers.append(str(primer_seq))
    return primers

# Function to analyze primer secondary structure using NUPACK
def analyze_structure(primer):
    model = nupack.Model(material='dna', celsius=37)
    mfe_struct = nupack.mfe([primer], model=model)
    return mfe_struct[0].energy  # Gibbs Free Energy (∆G)

# Function to evaluate primers
def evaluate_primers(primers):
    evaluated = []
    for primer in primers:
        tm = mt.Tm_NN(primer)  # Nearest Neighbor Tm Calculation
        gc = gc_fraction(primer) * 100  # Convert gc_fraction to percentage
        hairpin = primer3.calc_hairpin(primer).tm  # Hairpin Melting Temperature
        homodimer = primer3.calc_homodimer(primer).tm  # Self-dimer Tm
        stability = analyze_structure(primer)  # Gibbs Free Energy (∆G)

        # Ensure RPA-compatible thermodynamic properties
        if 37 <= tm <= 42 and 30 <= gc <= 50 and hairpin < 35 and homodimer < 35:
            evaluated.append((primer, tm, gc, hairpin, homodimer, stability))

    return evaluated

# Function to process FASTA file and design RPA primers
def design_rpa_primers(fasta_file):
    fasta_sequences = read_fasta(fasta_file)
    best_primers = []

    for sequence in fasta_sequences:
        candidates = generate_primers(sequence)
        evaluated = evaluate_primers(candidates)

        # Rank by lowest Gibbs Free Energy (most stable)
        sorted_primers = sorted(evaluated, key=lambda x: x[5])
        best_primers.extend(sorted_primers[:5])  # Keep the top 5 per sequence

    return best_primers

# Main execution block
if __name__ == "__main__":
    fasta_file = "genome.fasta"  # Change this to your FASTA file
    if not os.path.exists(fasta_file):
        print(f"Error: {fasta_file} not found!")
    else:
        best_primers = design_rpa_primers(fasta_file)
        for primer in best_primers:
            print(f"Primer: {primer[0]}, Tm: {primer[1]:.2f}, GC: {primer[2]:.2f}%, ΔG: {primer[5]:.2f}")
