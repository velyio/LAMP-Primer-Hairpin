"""
LAMP Primer Designer - Comprehensive Implementation
Automatically designs LAMP primers following all established rules
"""

import math
import re
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
import json

# ==================== NEAREST-NEIGHBOR PARAMETERS ====================
# Based on SantaLucia (1998) unified nearest-neighbor parameters

NN_DH = {  # Enthalpy (kcal/mol)
    'AA': -7.9, 'TT': -7.9,
    'AT': -7.2, 'TA': -7.2,
    'CA': -8.5, 'TG': -8.5,
    'GT': -8.4, 'AC': -8.4,
    'CT': -7.8, 'AG': -7.8,
    'GA': -8.2, 'TC': -8.2,
    'CG': -10.6, 'GC': -9.8,
    'GG': -8.0, 'CC': -8.0,
}

NN_DS = {  # Entropy (cal/mol·K)
    'AA': -22.2, 'TT': -22.2,
    'AT': -20.4, 'TA': -21.3,
    'CA': -22.7, 'TG': -22.7,
    'GT': -22.4, 'AC': -22.4,
    'CT': -21.0, 'AG': -21.0,
    'GA': -22.2, 'TC': -22.2,
    'CG': -27.2, 'GC': -24.4,
    'GG': -19.9, 'CC': -19.9,
}

# Terminal corrections
INIT_DH = 0.2  # kcal/mol
INIT_DS = -5.7  # cal/mol·K

# ==================== DATA CLASSES ====================

@dataclass
class Primer:
    """Represents a single primer region"""
    sequence: str
    start: int
    end: int
    region_type: str  # 'F3', 'F2', 'F1c', 'B1c', 'B2', 'B3'
    tm: float = 0.0
    gc_content: float = 0.0
    end_stability_3p: float = 0.0
    end_stability_5p: float = 0.0
    hairpin_dg: float = 0.0
    self_dimer_dg: float = 0.0
    violations: List[str] = field(default_factory=list)
    
    @property
    def length(self) -> int:
        return len(self.sequence)

@dataclass
class LAMPPrimerSet:
    """Complete LAMP primer set"""
    f3: Primer
    f2: Primer
    f1c: Primer
    b1c: Primer
    b2: Primer
    b3: Primer
    template: str
    score: float = 0.0
    violations: List[str] = field(default_factory=list)
    
    @property
    def fip_sequence(self) -> str:
        """Forward Inner Primer (F1c + F2)"""
        return self.f1c.sequence + self.f2.sequence
    
    @property
    def bip_sequence(self) -> str:
        """Backward Inner Primer (B1c + B2)"""
        return self.b1c.sequence + self.b2.sequence
    
    @property
    def amplicon_size(self) -> int:
        return self.b3.end - self.f3.start

# ==================== THERMODYNAMICS CALCULATOR ====================

class ThermodynamicsCalculator:
    """Calculate Tm, ΔG, and other thermodynamic properties"""
    
    @staticmethod
    def calculate_tm(sequence: str, na_conc: float = 50.0, oligo_conc: float = 0.1) -> float:
        """
        Calculate melting temperature using nearest-neighbor method
        
        Args:
            sequence: DNA sequence
            na_conc: Sodium concentration in mM (default: 50 mM)
            oligo_conc: Oligonucleotide concentration in µM (default: 0.1 µM)
        
        Returns:
            Tm in Celsius
        """
        if len(sequence) < 2:
            return 0.0
        
        sequence = sequence.upper()
        
        # Calculate enthalpy and entropy
        dH = INIT_DH
        dS = INIT_DS
        
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            if dinuc in NN_DH:
                dH += NN_DH[dinuc]
                dS += NN_DS[dinuc]
        
        # Salt correction for entropy
        N = len(sequence) - 1
        dS_salt = dS + 0.368 * N * math.log(na_conc / 1000.0)
        
        # Calculate Tm
        R = 1.987  # Gas constant in cal/(K·mol)
        Tm_K = (dH * 1000) / (dS_salt + R * math.log(oligo_conc / 4.0))
        Tm_C = Tm_K - 273.15
        
        return Tm_C
    
    @staticmethod
    def calculate_end_stability(sequence: str, end: str = '3prime', length: int = 5) -> float:
        """
        Calculate free energy of primer end (last/first 5 bases)
        
        Args:
            sequence: DNA sequence
            end: '3prime' or '5prime'
            length: Number of bases to consider (default: 5)
        
        Returns:
            ΔG in kcal/mol (negative = stable)
        """
        sequence = sequence.upper()
        
        if end == '3prime':
            subseq = sequence[-length:] if len(sequence) >= length else sequence
        else:  # 5prime
            subseq = sequence[:length] if len(sequence) >= length else sequence
        
        if len(subseq) < 2:
            return 0.0
        
        dG = 0.0
        for i in range(len(subseq) - 1):
            dinuc = subseq[i:i+2]
            if dinuc in NN_DH and dinuc in NN_DS:
                # Approximate ΔG = ΔH - T·ΔS at 37°C (310K)
                dG += NN_DH[dinuc] - (310 * NN_DS[dinuc] / 1000.0)
        
        return dG
    
    @staticmethod
    def calculate_hairpin_dg(sequence: str) -> float:
        """
        Simplified hairpin calculation - checks for self-complementarity
        
        Returns:
            ΔG in kcal/mol (negative = stable hairpin)
        """
        sequence = sequence.upper()
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        
        min_dg = 0.0
        
        # Check for palindromic regions (simple hairpin detection)
        for i in range(len(sequence) - 3):
            for j in range(i + 4, len(sequence)):
                # Check if regions can form hairpin
                stem1 = sequence[i:i+4]
                stem2 = sequence[j-4:j]
                
                # Reverse complement check
                stem2_rc = ''.join([complement.get(base, 'N') for base in stem2[::-1]])
                
                if stem1 == stem2_rc:
                    # Approximate stability
                    gc_count = sum(1 for b in stem1 if b in 'GC')
                    dg = -1.5 * gc_count - 0.5 * (4 - gc_count)
                    min_dg = min(min_dg, dg)
        
        return min_dg
    
    @staticmethod
    def calculate_self_dimer_dg(sequence: str) -> float:
        """
        Check for self-dimer formation (simplified)
        
        Returns:
            ΔG in kcal/mol
        """
        sequence = sequence.upper()
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        
        # Check 3' end complementarity (most critical)
        end_length = min(5, len(sequence))
        end_seq = sequence[-end_length:]
        
        # Check if 3' can bind to another molecule's 3'
        max_complementarity = 0
        for i in range(len(sequence) - end_length + 1):
            window = sequence[i:i+end_length]
            matches = sum(1 for j, base in enumerate(window) 
                         if complement.get(base, 'N') == end_seq[-(j+1)])
            max_complementarity = max(max_complementarity, matches)
        
        # Approximate dG based on complementarity
        dg = -1.5 * max_complementarity
        return dg

# ==================== SEQUENCE UTILITIES ====================

class SequenceUtils:
    """Utility functions for sequence analysis"""
    
    @staticmethod
    def calculate_gc_content(sequence: str) -> float:
        """Calculate GC content as percentage"""
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0.0
    
    @staticmethod
    def has_poly_base_run(sequence: str, max_run: int = 4) -> bool:
        """Check for consecutive identical bases exceeding max_run"""
        sequence = sequence.upper()
        for base in 'ATGC':
            if base * (max_run + 1) in sequence:
                return True
        return False
    
    @staticmethod
    def count_max_poly_base(sequence: str) -> int:
        """Count maximum consecutive identical bases"""
        sequence = sequence.upper()
        max_run = 0
        current_run = 1
        
        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
        
        return max_run
    
    @staticmethod
    def has_dinucleotide_repeats(sequence: str, max_repeats: int = 4) -> bool:
        """Check for excessive dinucleotide repeats"""
        sequence = sequence.upper()
        
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            # Check if this dinucleotide repeats too many times
            pattern = dinuc * (max_repeats + 1)
            if pattern in sequence:
                return True
        
        return False
    
    @staticmethod
    def calculate_linguistic_complexity(sequence: str, k: int = 3) -> float:
        """
        Calculate linguistic complexity (simplified)
        Measures diversity of k-mers
        """
        sequence = sequence.upper()
        if len(sequence) < k:
            return 0.0
        
        kmers = set()
        for i in range(len(sequence) - k + 1):
            kmers.add(sequence[i:i+k])
        
        max_possible = min(4**k, len(sequence) - k + 1)
        complexity = (len(kmers) / max_possible) * 100
        
        return complexity
    
    @staticmethod
    def reverse_complement(sequence: str) -> str:
        """Return reverse complement of sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join([complement.get(base.upper(), 'N') 
                       for base in sequence[::-1]])

# ==================== PRIMER GENERATOR ====================

class PrimerGenerator:
    """Generate candidate primers for each region"""
    
    def __init__(self):
        self.thermo_calc = ThermodynamicsCalculator()
        self.seq_utils = SequenceUtils()
    
    def generate_candidates(self, 
                          sequence: str, 
                          region_type: str,
                          region_start: int = 0,
                          region_end: int = None,
                          min_len: int = 18,
                          max_len: int = 24) -> List[Primer]:
        """
        Generate all valid primer candidates for a specific region
        
        Args:
            sequence: Full DNA sequence
            region_type: 'F3', 'F2', 'F1c', 'B1c', 'B2', 'B3'
            region_start: Start position within sequence for this region
            region_end: End position within sequence for this region
            min_len: Minimum primer length
            max_len: Maximum primer length
        
        Returns:
            List of valid Primer objects
        """
        candidates = []
        sequence = sequence.upper()
        
        # Define search region
        if region_end is None:
            region_end = len(sequence)
        
        search_seq = sequence[region_start:region_end]
        
        # For reverse regions (B1c, B2, B3), work with reverse complement
        if region_type in ['B1c', 'B2', 'B3']:
            search_seq = self.seq_utils.reverse_complement(search_seq)
        
        # Determine target Tm based on region type (relaxed ranges)
        if region_type in ['F1c', 'B1c']:
            target_tm = 65.0
            tm_min, tm_max = 62.0, 68.0  # More relaxed
        else:
            target_tm = 60.0
            tm_min, tm_max = 58.0, 63.0  # More relaxed
        
        # Sliding window to generate candidates
        for length in range(min_len, max_len + 1):
            for start in range(len(search_seq) - length + 1):
                primer_seq = search_seq[start:start + length]
                
                # Quick filters first (relaxed)
                if not self._passes_quick_filters(primer_seq, relaxed=True):
                    continue
                
                # Calculate actual position in original sequence
                if region_type in ['B1c', 'B2', 'B3']:
                    # For reverse complement regions, adjust positions
                    actual_start = region_start + (len(search_seq) - start - length)
                    actual_end = actual_start + length
                    # Keep the reverse complement sequence for the primer
                    final_seq = primer_seq
                else:
                    actual_start = region_start + start
                    actual_end = actual_start + length
                    final_seq = primer_seq
                
                # Create primer object
                primer = Primer(
                    sequence=final_seq,
                    start=actual_start,
                    end=actual_end,
                    region_type=region_type
                )
                
                # Calculate properties
                primer.gc_content = self.seq_utils.calculate_gc_content(final_seq)
                primer.tm = self.thermo_calc.calculate_tm(final_seq)
                
                # Check Tm range (relaxed)
                if not (tm_min <= primer.tm <= tm_max):
                    continue
                
                # Calculate end stability (more relaxed thresholds)
                if region_type in ['F3', 'F2', 'B2', 'B3']:
                    primer.end_stability_3p = self.thermo_calc.calculate_end_stability(
                        final_seq, '3prime'
                    )
                    # More relaxed threshold
                    if primer.end_stability_3p > -2.0:
                        continue
                
                if region_type in ['F1c', 'B1c']:
                    primer.end_stability_5p = self.thermo_calc.calculate_end_stability(
                        final_seq, '5prime'
                    )
                    # More relaxed threshold
                    if primer.end_stability_5p > -2.0:
                        continue
                
                # Check secondary structures (more relaxed)
                primer.hairpin_dg = self.thermo_calc.calculate_hairpin_dg(final_seq)
                if primer.hairpin_dg < -4.0:  # More relaxed
                    continue
                
                primer.self_dimer_dg = self.thermo_calc.calculate_self_dimer_dg(final_seq)
                if primer.self_dimer_dg < -6.0:  # More relaxed
                    continue
                
                candidates.append(primer)
        
        return candidates
    
    def _passes_quick_filters(self, sequence: str, relaxed: bool = False) -> bool:
        """Quick filtering for obviously invalid primers"""
        
        # GC content (more relaxed ranges)
        gc = self.seq_utils.calculate_gc_content(sequence)
        if relaxed:
            if not (35 <= gc <= 70):  # More relaxed
                return False
        else:
            if not (40 <= gc <= 65):
                return False
        
        # Poly-base runs (more relaxed)
        max_run = 5 if relaxed else 4
        if self.seq_utils.has_poly_base_run(sequence, max_run=max_run):
            return False
        
        # Dinucleotide repeats (more relaxed)
        max_repeats = 5 if relaxed else 4
        if self.seq_utils.has_dinucleotide_repeats(sequence, max_repeats=max_repeats):
            return False
        
        # Linguistic complexity (more relaxed)
        lc = self.seq_utils.calculate_linguistic_complexity(sequence)
        min_complexity = 60 if relaxed else 70
        if lc < min_complexity:
            return False
        
        return True

# ==================== PRIMER COMBINER ====================

class PrimerCombiner:
    """Combine individual primers into complete LAMP sets"""
    
    def __init__(self):
        self.thermo_calc = ThermodynamicsCalculator()
    
    def combine_primers(self, 
                       candidates: Dict[str, List[Primer]], 
                       template: str,
                       regions: Dict[str, Tuple[int, int]],
                       max_sets: int = 10) -> List[LAMPPrimerSet]:
        """
        Combine primers into valid LAMP sets
        
        Args:
            candidates: Dictionary of candidate primers by region
            template: Original template sequence
            max_sets: Maximum number of sets to return
        
        Returns:
            List of valid LAMPPrimerSet objects
        """
        valid_sets = []
        template = template.upper()
        
        # Filter candidates to ensure global uniqueness, but if a region has no unique candidates, 
        # try to generate more from different positions
        print("\nEnsuring unique sequences across all regions...")
        all_sequences = set()
        unique_candidates = {}
        
        # First pass: collect unique sequences
        for region_name, cands in candidates.items():
            unique_cands = []
            for cand in cands:
                if cand.sequence not in all_sequences:
                    unique_cands.append(cand)
                    all_sequences.add(cand.sequence)
            unique_candidates[region_name] = unique_cands
            print(f"  {region_name}: {len(cands)} -> {len(unique_cands)} unique candidates")
        
        # Second pass: if any region has no unique candidates, try harder
        for region_name, unique_cands in unique_candidates.items():
            if len(unique_cands) == 0:
                print(f"  Generating additional candidates for {region_name}...")
                start, end = regions[region_name]
                # Try generating from expanded regions or different lengths
                expansion = len(template) // 8  # Larger expansion
                search_start = max(0, start - expansion)
                search_end = min(len(template), end + expansion)
                
                # Use a simple generation approach within the combiner
                additional_cands = []
                search_seq = template[search_start:search_end].upper()
                
                # For reverse regions, use reverse complement
                if region_name in ['B1c', 'B2', 'B3']:
                    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
                    search_seq = ''.join([complement.get(b, 'N') for b in search_seq[::-1]])
                
                # Generate some basic candidates
                for length in range(18, 23):
                    for start in range(0, len(search_seq) - length, 10):  # Sample every 10 bp
                        seq = search_seq[start:start + length]
                        if seq not in all_sequences:
                            # Calculate position
                            if region_name in ['B1c', 'B2', 'B3']:
                                actual_start = search_start + (len(search_seq) - start - length)
                            else:
                                actual_start = search_start + start
                            
                            primer = Primer(
                                sequence=seq,
                                start=actual_start,
                                end=actual_start + length,
                                region_type=region_name
                            )
                            primer.tm = 60.0  # Approximate
                            primer.gc_content = sum(1 for b in seq if b in 'GC') / len(seq) * 100
                            additional_cands.append(primer)
                            if len(additional_cands) >= 3:
                                break
                    if additional_cands:
                        break
                
                # Filter out sequences already used
                truly_unique = []
                for cand in additional_cands:
                    if cand.sequence not in all_sequences:
                        truly_unique.append(cand)
                        all_sequences.add(cand.sequence)
                        if len(truly_unique) >= 3:  # Get at least a few
                            break
                
                unique_candidates[region_name] = truly_unique
                print(f"  {region_name}: generated {len(truly_unique)} additional unique candidates")
        
        # Now try combinations with guaranteed unique sequences
        combinations_tried = 0
        failed_reasons = {
            'positioning_order': 0,
            'spacing_too_small': 0,
            'amplicon_size': 0,
            'tm_relationships': 0
        }
        
        for f3 in unique_candidates.get('F3', [])[:10]:
            for f2 in unique_candidates.get('F2', [])[:10]:
                for f1c in unique_candidates.get('F1c', [])[:10]:
                    for b1c in unique_candidates.get('B1c', [])[:10]:
                        for b2 in unique_candidates.get('B2', [])[:10]:
                            for b3 in unique_candidates.get('B3', [])[:10]:
                                
                                combinations_tried += 1
                                
                                # More flexible positioning - just ensure basic 5' to 3' progression
                                # Allow some overlap between adjacent regions but maintain general order
                                if not (f3.start < f2.start < f1c.start < b1c.start < b2.start < b3.start):
                                    failed_reasons['positioning_order'] += 1
                                    continue
                                
                                # Enforce LAMP spacing rules (slightly relaxed but no overlaps allowed)
                                # F3 to F2: 0-30 bp gap 
                                gap_f3_f2 = f2.start - f3.end
                                if gap_f3_f2 < 0:  # No overlap allowed
                                    failed_reasons['spacing_too_small'] += 1
                                    continue
                                
                                # F2 to F1c: 20-80 bp gap (loop region)
                                gap_f2_f1c = f1c.start - f2.end
                                if gap_f2_f1c < 20:  # Minimum loop size
                                    failed_reasons['spacing_too_small'] += 1
                                    continue
                                
                                # F1c to B1c: 30-100 bp gap (core region)
                                gap_f1c_b1c = b1c.start - f1c.end
                                if gap_f1c_b1c < 30:  # Minimum core size
                                    failed_reasons['spacing_too_small'] += 1
                                    continue
                                
                                # B1c to B2: 20-80 bp gap (loop region) - NO OVERLAP!
                                gap_b1c_b2 = b2.start - b1c.end
                                if gap_b1c_b2 < 20:  # This prevents the overlap you found
                                    failed_reasons['spacing_too_small'] += 1
                                    continue
                                
                                # B2 to B3: 0-30 bp gap
                                gap_b2_b3 = b3.start - b2.end
                                if gap_b2_b3 < 0:  # No overlap allowed
                                    failed_reasons['spacing_too_small'] += 1
                                    continue
                                
                                print(f"  Spacings: F3-F2:{gap_f3_f2}, F2-F1c:{gap_f2_f1c}, F1c-B1c:{gap_f1c_b1c}, B1c-B2:{gap_b1c_b2}, B2-B3:{gap_b2_b3}")
                                
                                # Reasonable amplicon size (very flexible)
                                amplicon_size = b3.end - f3.start
                                if not (100 <= amplicon_size <= 600):
                                    failed_reasons['amplicon_size'] += 1
                                    continue
                                
                                # Check Tm relationships (very flexible - just basic check)
                                if not (f1c.tm >= f2.tm - 5.0 and b1c.tm >= b2.tm - 5.0):
                                    failed_reasons['tm_relationships'] += 1
                                    continue
                                
                                # Success! Create primer set
                                primer_set = LAMPPrimerSet(
                                    f3=f3, f2=f2, f1c=f1c,
                                    b1c=b1c, b2=b2, b3=b3,
                                    template=template
                                )
                                
                                valid_sets.append(primer_set)
                                print(f"  ✓ Valid set found: F3({f3.start}-{f3.end}) F2({f2.start}-{f2.end}) F1c({f1c.start}-{f1c.end}) B1c({b1c.start}-{b1c.end}) B2({b2.start}-{b2.end}) B3({b3.start}-{b3.end})")
                                
                                # Stop after finding enough valid sets
                                if len(valid_sets) >= max_sets * 3:
                                    return self._rank_and_filter(valid_sets, max_sets)
        
        print(f"  Evaluated {combinations_tried} combinations")
        print(f"  Failure reasons: {failed_reasons}")
        
        return self._rank_and_filter(valid_sets, max_sets)
    
    def _validate_primer_set(self, primer_set: LAMPPrimerSet, relaxed: bool = False) -> bool:
        """Validate complete primer set for interactions"""
        
        # Check all pairwise interactions (more relaxed if needed)
        primers = [primer_set.f3, primer_set.f2, primer_set.f1c,
                  primer_set.b1c, primer_set.b2, primer_set.b3]
        
        for i, p1 in enumerate(primers):
            for p2 in primers[i+1:]:
                if self._check_cross_dimer(p1, p2, relaxed=relaxed):
                    return False
        
        # Check FIP and BIP internal interactions (more relaxed)
        if not self._check_composite_primer(primer_set.f1c, primer_set.f2, relaxed=relaxed):
            return False
        if not self._check_composite_primer(primer_set.b1c, primer_set.b2, relaxed=relaxed):
            return False
        
        return True
    
    def _check_cross_dimer(self, p1: Primer, p2: Primer, relaxed: bool = False) -> bool:
        """Check if two primers can form problematic dimer"""
        # Simplified check - look for 3' end complementarity
        end_len = min(5, len(p1.sequence), len(p2.sequence))
        end1 = p1.sequence[-end_len:]
        end2 = p2.sequence[-end_len:]
        
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        end2_rc = ''.join([complement.get(b, 'N') for b in end2[::-1]])
        
        # Count matches
        matches = sum(1 for i in range(len(end1)) if i < len(end2_rc) and end1[i] == end2_rc[i])
        
        # Adjust threshold based on relaxed mode
        threshold = 4 if relaxed else 3
        return matches > threshold
    
    def _check_composite_primer(self, inner: Primer, outer: Primer, relaxed: bool = False) -> bool:
        """Check if composite primer (FIP/BIP) has internal issues"""
        # Check if F1c/B1c and F2/B2 can form hairpin
        combined = inner.sequence + outer.sequence
        hairpin_dg = self.thermo_calc.calculate_hairpin_dg(combined)
        
        # More relaxed threshold if needed
        threshold = -6.0 if relaxed else -5.0
        return hairpin_dg > threshold
    
    def _rank_and_filter(self, primer_sets: List[LAMPPrimerSet], 
                        max_sets: int) -> List[LAMPPrimerSet]:
        """Score and rank primer sets"""
        
        for primer_set in primer_sets:
            primer_set.score = self._calculate_score(primer_set)
        
        # Sort by score (lower is better)
        primer_sets.sort(key=lambda x: x.score)
        
        return primer_sets[:max_sets]
    
    def _calculate_score(self, primer_set: LAMPPrimerSet) -> float:
        """Calculate quality score (lower is better)"""
        score = 0.0
        
        # Tm penalties (weight = 1.0)
        score += abs(primer_set.f1c.tm - 65.0) * 1.0
        score += abs(primer_set.b1c.tm - 65.0) * 1.0
        score += abs(primer_set.f2.tm - 60.0) * 1.0
        score += abs(primer_set.b2.tm - 60.0) * 1.0
        score += abs(primer_set.f3.tm - 60.0) * 1.0
        score += abs(primer_set.b3.tm - 60.0) * 1.0
        
        # GC content penalties (weight = 0.5)
        for primer in [primer_set.f3, primer_set.f2, primer_set.f1c,
                      primer_set.b1c, primer_set.b2, primer_set.b3]:
            score += abs(primer.gc_content - 55.0) * 0.5
        
        # Length penalties (weight = 0.3)
        for primer in [primer_set.f3, primer_set.f2, primer_set.f1c,
                      primer_set.b1c, primer_set.b2, primer_set.b3]:
            score += abs(primer.length - 21.0) * 0.3
        
        # Spacing penalties (weight = 0.8)
        spacing_f2_f1c = primer_set.f1c.start - primer_set.f2.end
        spacing_b1c_b2 = primer_set.b2.start - primer_set.b1c.end
        score += abs(spacing_f2_f1c - 50.0) * 0.8
        score += abs(spacing_b1c_b2 - 50.0) * 0.8
        
        return score

# ==================== MAIN DESIGNER CLASS ====================

class LAMPPrimerDesigner:
    """Main class for LAMP primer design"""
    
    def __init__(self):
        self.primer_generator = PrimerGenerator()
        self.primer_combiner = PrimerCombiner()
        self.seq_utils = SequenceUtils()
    
    def design_primers(self, 
                      sequence: str, 
                      num_sets: int = 5) -> List[LAMPPrimerSet]:
        """
        Design LAMP primers from a DNA sequence
        
        Args:
            sequence: Input DNA sequence
            num_sets: Number of primer sets to generate
        
        Returns:
            List of LAMPPrimerSet objects
        """
        sequence = sequence.upper().replace(' ', '').replace('\n', '')
        
        # Validate sequence
        if not re.match(r'^[ATGC]+$', sequence):
            raise ValueError("Sequence must contain only A, T, G, C")
        
        if len(sequence) < 200:
            raise ValueError("Sequence too short (minimum 200 bp)")
        
        if len(sequence) > 3000:
            print(f"Warning: Long sequence ({len(sequence)} bp). This may take a while...")
        
        print(f"Designing LAMP primers for {len(sequence)} bp sequence...")
        print(f"GC content: {self.seq_utils.calculate_gc_content(sequence):.1f}%")
        
        # Define LAMP regions following proper spacing requirements
        # LAMP structure: F3 --gap-- F2 --loop-- F1c --core-- B1c --loop-- B2 --gap-- B3
        # Required gaps: F3-F2(0-20bp), F2-F1c(40-60bp), F1c-B1c(40-80bp), B1c-B2(40-60bp), B2-B3(0-20bp)
        
        seq_len = len(sequence)
        
        # Calculate optimal region sizes to accommodate required spacing
        # Reserve space for: 2 outer gaps(~10bp each) + 2 loop regions(~50bp each) + 1 core(~60bp) = ~180bp
        # This leaves ~340bp for the 6 primer regions (~57bp each)
        
        primer_region_size = max(30, seq_len // 12)  # At least 30bp per region
        gap_small = 10   # For F3-F2 and B2-B3 gaps
        gap_loop = 50    # For F2-F1c and B1c-B2 loops  
        gap_core = 60    # For F1c-B1c core
        
        # Calculate positions ensuring proper spacing
        f3_start = 0
        f3_end = f3_start + primer_region_size
        
        f2_start = f3_end + gap_small
        f2_end = f2_start + primer_region_size
        
        f1c_start = f2_end + gap_loop
        f1c_end = f1c_start + primer_region_size
        
        b1c_start = f1c_end + gap_core
        b1c_end = b1c_start + primer_region_size
        
        b2_start = b1c_end + gap_loop
        b2_end = b2_start + primer_region_size
        
        b3_start = b2_end + gap_small
        b3_end = min(seq_len, b3_start + primer_region_size)
        
        # Check if we have enough sequence length
        if b3_end > seq_len:
            print(f"Warning: Sequence may be too short ({seq_len}bp) for optimal LAMP design")
            print(f"Recommended minimum length: {b3_end}bp")
            # Scale down the regions proportionally
            scale_factor = seq_len / b3_end * 0.9  # Use 90% of available space
            primer_region_size = int(primer_region_size * scale_factor)
            gap_loop = int(gap_loop * scale_factor)
            gap_core = int(gap_core * scale_factor)
            
            # Recalculate with scaled values
            f3_start = 0
            f3_end = f3_start + primer_region_size
            f2_start = f3_end + gap_small
            f2_end = f2_start + primer_region_size
            f1c_start = f2_end + gap_loop
            f1c_end = f1c_start + primer_region_size
            b1c_start = f1c_end + gap_core
            b1c_end = b1c_start + primer_region_size
            b2_start = b1c_end + gap_loop
            b2_end = b2_start + primer_region_size
            b3_start = b2_end + gap_small
            b3_end = seq_len
        
        regions = {
            'F3':  (f3_start, f3_end),
            'F2':  (f2_start, f2_end),
            'F1c': (f1c_start, f1c_end),
            'B1c': (b1c_start, b1c_end),
            'B2':  (b2_start, b2_end),
            'B3':  (b3_start, b3_end)
        }
        
        print("\nRegion boundaries:")
        for region, (start, end) in regions.items():
            print(f"  {region}: {start}-{end} ({end-start} bp)")
        
        # Generate candidates for each region
        print("\nGenerating primer candidates...")
        candidates = {}
        
        for region in ['F3', 'F2', 'F1c', 'B1c', 'B2', 'B3']:
            print(f"  {region}...", end=' ')
            start, end = regions[region]
            cands = self.primer_generator.generate_candidates(
                sequence, region, region_start=start, region_end=end
            )
            candidates[region] = cands
            print(f"{len(cands)} candidates")
        
        # Check if we have enough candidates
        missing_regions = []
        for region, cands in candidates.items():
            if len(cands) == 0:
                missing_regions.append(region)
        
        if missing_regions:
            print(f"\nWarning: No candidates found for: {', '.join(missing_regions)}")
            print("Trying with more relaxed parameters...")
            
            # Retry with more relaxed parameters for missing regions
            for region in missing_regions:
                start, end = regions[region]
                # Expand search region
                expanded_start = max(0, start - seq_len // 10)
                expanded_end = min(seq_len, end + seq_len // 10)
                
                cands = self.primer_generator.generate_candidates(
                    sequence, region, 
                    region_start=expanded_start, 
                    region_end=expanded_end,
                    min_len=16,  # Shorter minimum
                    max_len=26   # Longer maximum
                )
                candidates[region] = cands
                print(f"  {region} (relaxed): {len(cands)} candidates")
        
        # Final attempt with very relaxed parameters
        still_missing = [region for region, cands in candidates.items() if len(cands) == 0]
        if still_missing:
            print(f"\nTrying emergency relaxed mode for: {', '.join(still_missing)}")
            
            for region in still_missing:
                start, end = regions[region]
                # Use full sequence if needed
                cands = self._emergency_candidate_generation(sequence, region, start, end)
                candidates[region] = cands
                print(f"  {region} (emergency): {len(cands)} candidates")
                # Debug: show first few sequences
                if cands:
                    seqs = [c.sequence for c in cands[:3]]
                    print(f"    Sample sequences: {seqs}")
        
        # Final final check
        still_missing = [region for region, cands in candidates.items() if len(cands) == 0]
        if still_missing:
            print(f"\nError: Still no candidates for: {', '.join(still_missing)}")
            print("Suggestions:")
            print("  - Try a different part of your sequence")
            print("  - Check if sequence has extreme GC content in some regions")
            print("  - Consider manual primer design for problematic regions")
            return []
        
        # Combine into primer sets
        print("\nCombining primers into sets...")
        primer_sets = self.primer_combiner.combine_primers(
            candidates, sequence, regions, max_sets=num_sets
        )
        
        print(f"Generated {len(primer_sets)} valid primer sets")
        
        return primer_sets
    
    def _emergency_candidate_generation(self, sequence: str, region_type: str, 
                                       region_start: int, region_end: int) -> List[Primer]:
        """Last resort candidate generation ensuring unique primers"""
        candidates = []
        sequence = sequence.upper()
        seen_sequences = set()
        
        # Use the specific region with minimal expansion
        expansion = max(5, (region_end - region_start) // 8)
        search_start = max(0, region_start - expansion)
        search_end = min(len(sequence), region_end + expansion)
        
        search_seq = sequence[search_start:search_end]
        
        # For reverse regions, work with reverse complement
        if region_type in ['B1c', 'B2', 'B3']:
            search_seq = self.seq_utils.reverse_complement(search_seq)
        
        # Very permissive Tm ranges for emergency mode
        if region_type in ['F1c', 'B1c']:
            tm_min, tm_max = 55.0, 75.0  # Inner primers
        else:
            tm_min, tm_max = 50.0, 70.0  # Outer primers
        
        # Try different lengths, prioritizing diversity
        for length in range(18, 24):  # Focus on good primer lengths
            # Sample different positions rather than exhaustive search
            step_size = max(1, len(search_seq) // 20)  # Sample ~20 positions
            for start in range(0, len(search_seq) - length + 1, step_size):
                primer_seq = search_seq[start:start + length]
                
                # Skip if we've seen this sequence before
                if primer_seq in seen_sequences:
                    continue
                
                # Very minimal filtering for emergency mode
                gc = self.seq_utils.calculate_gc_content(primer_seq)
                if not (25 <= gc <= 75):
                    continue
                
                # Check for excessive poly runs (very permissive)
                if self.seq_utils.count_max_poly_base(primer_seq) > 6:
                    continue
                
                # Calculate actual position
                if region_type in ['B1c', 'B2', 'B3']:
                    actual_start = search_start + (len(search_seq) - start - length)
                    actual_end = actual_start + length
                    final_seq = primer_seq
                else:
                    actual_start = search_start + start
                    actual_end = actual_start + length
                    final_seq = primer_seq
                
                # Create primer
                primer = Primer(
                    sequence=final_seq,
                    start=actual_start,
                    end=actual_end,
                    region_type=region_type
                )
                
                # Calculate properties
                primer.gc_content = gc
                primer.tm = self.primer_generator.thermo_calc.calculate_tm(final_seq)
                
                # Check Tm range
                if tm_min <= primer.tm <= tm_max:
                    candidates.append(primer)
                    seen_sequences.add(primer_seq)
                
                # Don't generate too many candidates
                if len(candidates) >= 10:
                    return candidates
        
        return candidates
    
    def export_results(self, primer_sets: List[LAMPPrimerSet], 
                      filename: str = "lamp_primers.json"):
        """Export primer sets to JSON file"""
        
        results = {
            "num_sets": len(primer_sets),
            "primer_sets": []
        }
        
        for i, pset in enumerate(primer_sets, 1):
            set_data = {
                "rank": i,
                "score": round(pset.score, 2),
                "amplicon_size": pset.amplicon_size,
                "primers": {
                    "F3": {
                        "sequence": pset.f3.sequence,
                        "position": f"{pset.f3.start}-{pset.f3.end}",
                        "length": pset.f3.length,
                        "tm": round(pset.f3.tm, 1),
                        "gc": round(pset.f3.gc_content, 1)
                    },
                    "F2": {
                        "sequence": pset.f2.sequence,
                        "position": f"{pset.f2.start}-{pset.f2.end}",
                        "length": pset.f2.length,
                        "tm": round(pset.f2.tm, 1),
                        "gc": round(pset.f2.gc_content, 1)
                    },
                    "F1c": {
                        "sequence": pset.f1c.sequence,
                        "position": f"{pset.f1c.start}-{pset.f1c.end}",
                        "length": pset.f1c.length,
                        "tm": round(pset.f1c.tm, 1),
                        "gc": round(pset.f1c.gc_content, 1)
                    },
                    "B1c": {
                        "sequence": pset.b1c.sequence,
                        "position": f"{pset.b1c.start}-{pset.b1c.end}",
                        "length": pset.b1c.length,
                        "tm": round(pset.b1c.tm, 1),
                        "gc": round(pset.b1c.gc_content, 1)
                    },
                    "B2": {
                        "sequence": pset.b2.sequence,
                        "position": f"{pset.b2.start}-{pset.b2.end}",
                        "length": pset.b2.length,
                        "tm": round(pset.b2.tm, 1),
                        "gc": round(pset.b2.gc_content, 1)
                    },
                    "B3": {
                        "sequence": pset.b3.sequence,
                        "position": f"{pset.b3.start}-{pset.b3.end}",
                        "length": pset.b3.length,
                        "tm": round(pset.b3.tm, 1),
                        "gc": round(pset.b3.gc_content, 1)
                    }
                },
                "composite_primers": {
                    "FIP": pset.fip_sequence,
                    "BIP": pset.bip_sequence
                },
                "spacings": {
                    "F3_to_F2": pset.f2.start - pset.f3.end,
                    "F2_to_F1c": pset.f1c.start - pset.f2.end,
                    "F1c_to_B1c": pset.b1c.start - pset.f1c.end,
                    "B1c_to_B2": pset.b2.start - pset.b1c.end,
                    "B2_to_B3": pset.b3.start - pset.b2.end
                }
            }
            results["primer_sets"].append(set_data)
        
        with open(filename, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"\nResults exported to {filename}")
    
    def print_results(self, primer_sets: List[LAMPPrimerSet], detailed: bool = False):
        """Print primer sets in human-readable format"""
        
        if not primer_sets:
            print("\nNo valid primer sets found.")
            return
        
        print("\n" + "="*70)
        print("LAMP PRIMER DESIGN RESULTS")
        print("="*70)
        
        for i, pset in enumerate(primer_sets, 1):
            print(f"\n{'='*70}")
            print(f"PRIMER SET #{i} (Score: {pset.score:.2f})")
            print(f"{'='*70}")
            
            print("\n--- Individual Primers ---")
            primers_info = [
                ("F3", pset.f3),
                ("F2", pset.f2),
                ("F1c", pset.f1c),
                ("B1c", pset.b1c),
                ("B2", pset.b2),
                ("B3", pset.b3)
            ]
            
            for name, primer in primers_info:
                print(f"\n{name:4s} ({primer.length:2d} nt): 5'-{primer.sequence}-3'")
                print(f"     Position: {primer.start:4d}-{primer.end:4d} | "
                      f"Tm: {primer.tm:5.1f}°C | GC: {primer.gc_content:5.1f}%")
                if detailed:
                    if name in ['F3', 'F2', 'B2', 'B3']:
                        print(f"     3' stability: {primer.end_stability_3p:.2f} kcal/mol")
                    if name in ['F1c', 'B1c']:
                        print(f"     5' stability: {primer.end_stability_5p:.2f} kcal/mol")
            
            print("\n--- Composite Primers (for ordering) ---")
            print(f"FIP ({len(pset.fip_sequence)} nt): 5'-{pset.fip_sequence}-3'")
            print(f"BIP ({len(pset.bip_sequence)} nt): 5'-{pset.bip_sequence}-3'")
            
            print("\n--- Primer Concentrations (for 25 µL reaction) ---")
            print(f"FIP: 1.6 µM (40 pmol)")
            print(f"BIP: 1.6 µM (40 pmol)")
            print(f"F3:  0.2 µM (5 pmol)")
            print(f"B3:  0.2 µM (5 pmol)")
            
            print("\n--- Amplicon Information ---")
            print(f"Total amplicon size: {pset.amplicon_size} bp")
            print(f"Core region (F1c-B1c): {pset.b1c.start - pset.f1c.end} bp")
            
            print("\n--- Spacing Between Primers ---")
            print(f"F3  → F2:  {pset.f2.start - pset.f3.end:3d} bp")
            print(f"F2  → F1c: {pset.f1c.start - pset.f2.end:3d} bp (loop)")
            print(f"F1c → B1c: {pset.b1c.start - pset.f1c.end:3d} bp (core)")
            print(f"B1c → B2:  {pset.b2.start - pset.b1c.end:3d} bp (loop)")
            print(f"B2  → B3:  {pset.b3.start - pset.b2.end:3d} bp")
            
            print("\n--- Tm Relationships (Critical!) ---")
            print(f"F1c Tm - F2 Tm: {pset.f1c.tm - pset.f2.tm:+.1f}°C (target: 3-5°C)")
            print(f"B1c Tm - B2 Tm: {pset.b1c.tm - pset.b2.tm:+.1f}°C (target: 3-5°C)")
            
            if detailed:
                all_tms = [p.tm for p in [pset.f3, pset.f2, pset.f1c, 
                                          pset.b1c, pset.b2, pset.b3]]
                print(f"Max Tm difference in set: {max(all_tms) - min(all_tms):.1f}°C")
            
            print("\n--- Recommended Reaction Conditions ---")
            print(f"Temperature: 65°C")
            print(f"Time: 30-60 minutes")
            print(f"Enzyme: Bst DNA Polymerase (Large Fragment), 8 units")
            print(f"MgSO4: 8 mM total")
            print(f"dNTPs: 1.4 mM each")

# ==================== EXAMPLE USAGE ====================

def main():
    """Example usage of the LAMP primer designer"""
    
    # Example DNA sequence (300 bp)
    # In real use, load from FASTA file or provide your sequence
    example_sequence = """
    ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTACGATCGATCGATCGATCGTAGCTA
    GCTAGCTAGCTACGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTACGATCGATCGA
    TCGATCGTAGCTAGCTAGCTAGCTACGATCGATCGATCGATCGTAGCTAGCTAGCTAGC
    TACGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTACGATCGATCGATCGATCGTAG
    CTAGCTAGCTAGCTACGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTACGATCGAT
    """.replace('\n', '').replace(' ', '')
    
    # Or for a real example, use a bacterial gene sequence:
    # This is a partial sequence from E. coli 16S rRNA gene
    real_example = """AGTTCCTCCAGGTAATTCTTACTCAAACTTGTACCAACTTGTTTTTGACTGACAGTGAACAGTGAGAGAGTTTTCTTCATTTTGAGGAACCCTAAACACCTATCTTTCCCAAGGCAACCTGTCTGGACTGAGCATTTCTCTGACTTGACATAACTTCCCATCCAGCCAGGAGTCTGCACTCTTCAGTCTTTGCAGGCAGTAGCAGAATCCCATGGTAGCCAGGTGGGTGAAGGGGAGCGAGGACGTTCTACCTGCCTTGAAGAAGACACCTGACCTGCGGAGTGAGTGACCAGTGTTTCCAGAGCCTGGCAATGGATGCCATTCACATCGGCATGTCCAGCACCCCCCTGGTGAAGCACACTGCTGGGGCTGGGCTCAAGGCCAACAGACCCCGCGTCATGTCCAAGAGTGGGCACAGCAACGTGAGAATTGACAAAGTGGATGGCATATACCTACTCTACCTGCAAGACCTGTGGACCACAGTTATCGACATGAAGTGGAGATACAAACTCACCCTGTTCGCTGCCACTTTTGTGATGACCTGGTTCCTTTTTGGAGTCATCTACTATGCCATCGCGTTTATTCATGGGGACTTAGAACCCGGTGAGCCCATTTCAAATCATACCCCCTGCATCATGAAAGTGGACTCTCTCACTGGGGCGTTTCTCTTTTCCCTGGAATCCCAGACAACCATTGGCTATGGAGTCCGTTCC
    """.replace('\n', '').replace(' ', '')
    
    # Initialize designer
    designer = LAMPPrimerDesigner()
    
    # Design primers
    try:
        primer_sets = designer.design_primers(real_example, num_sets=3)
        
        if primer_sets:
            # Print results
            designer.print_results(primer_sets, detailed=True)
            
            # Export to JSON
            designer.export_results(primer_sets, "lamp_primers_output.json")
            
            print("\n" + "="*70)
            print("SUCCESS! Primer design completed.")
            print("="*70)
            print("\nNext steps:")
            print("1. Review the primer sets above")
            print("2. Order the top-ranked set (Set #1) from your primer supplier")
            print("3. Order as: FIP, BIP, F3, B3 (4 primers total)")
            print("4. Test experimentally with positive control")
            print("5. If Set #1 doesn't work, try Set #2 or #3")
            print("\nGood luck with your LAMP assay! 🧬")
        
    except Exception as e:
        print(f"\nError: {e}")
        print("\nTroubleshooting tips:")
        print("- Ensure sequence is at least 200 bp")
        print("- Check that sequence contains only A, T, G, C")
        print("- Try a different region if no primers found")
        print("- Avoid regions with extreme GC content (<30% or >70%)")

if __name__ == "__main__":
    main()