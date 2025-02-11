
### This script extracts the consensus sequences from a multiple sequence alignment.

from pathlib import Path
import argparse

import Bio
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def get_consensus_seq(filename: Path | str, threshold: float = 0.7) -> SeqRecord:
    """Generate a consensus sequence from a multiple sequence alignment, using the majority base when ambiguous.

    Args:
        filename (Path): Path to the alignment file.
        threshold (float): The minimum fraction of sequences required to call a consensus at each position.

    Returns:
        SeqRecord: The consensus sequence as a SeqRecord object.
    """
    # Read the alignment from the file
    alignment = AlignIO.read(filename, "fasta")

    # Create a list of the sequences from the alignment
    sequences = [str(record.seq) for record in alignment]

    # Generate the consensus sequence
    consensus_seq = []
    for i in range(len(sequences[0])):  # Iterate over the positions
        # Count the bases at position i
        bases = [seq[i] for seq in sequences]
        base_counts = {base: bases.count(base) for base in set(bases)}

        # Check if any base exceeds the threshold
        total_count = len(bases)
        consensus_base = None
        for base, count in base_counts.items():
            if count / total_count >= threshold:
                consensus_base = base
                break

        # If no base exceeds the threshold, select the majority base
        if consensus_base is None:
            # Select the majority base (the most frequent)
            consensus_base = max(base_counts, key=base_counts.get)

        consensus_seq.append(consensus_base)

    # Get the filename without extension for the header
    filename_without_extension = Path(filename).stem

    # Create a SeqRecord object with the consensus sequence
    consensus_record = SeqRecord(
        Seq("".join(consensus_seq)),
        id=filename_without_extension,  # Use the filename (without extension) as the ID
        description=f"Consensus sequence with threshold={threshold}"
    )

    return consensus_record

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a consensus sequence from an alignment file.")
    parser.add_argument("filename", type=str, help="Path to the alignment file in FASTA format.")
    parser.add_argument("--threshold", type=float, default=0.7, help="Threshold for consensus calling (default: 0.7).")

    args = parser.parse_args()

    consensus_seq = get_consensus_seq(args.filename, args.threshold)

    # Generate the output filename
    input_path = Path(args.filename)
    output_filename = input_path.with_name(f"{input_path.stem}_consensus.fasta")

    # Save the consensus sequence to a FASTA file
    with open(output_filename, "w") as output_file:
        output_file.write(f">{consensus_seq.id}\n{consensus_seq.seq}\n")

    print(f"Consensus sequence saved to {output_filename}")
