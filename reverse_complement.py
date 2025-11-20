def reverse_complement(dna_sequence):
    # Define the complement mapping
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Reverse the DNA sequence
    reversed_sequence = dna_sequence[::-1]
    # Get the complement of the reversed sequence
    reverse_complement_seq = ''.join(complement[base] for base in reversed_sequence)
    return reverse_complement_seq

# Example usage
original_sequence = 'GGGTGCTGGACATGCCGAT'
print("Original Sequence:", original_sequence)
print("Reverse Sequence:", original_sequence[::-1])
print("Reverse Complement:", reverse_complement(original_sequence))
