"""
Mutation processing utilities.

This module contains functions for processing mutation strings,
extracting genomic positions, and validating mutation formats.
"""

import re
from typing import List
from interface import MutationType
from typing import Optional


def get_symbols_for_mutation_type(mutation_type: MutationType) -> List[str]:
    """Returns the list of symbols (amino acids or nucleotides) for the given mutation type.
    
    Args:
        mutation_type (MutationType): The type of mutation (AMINO_ACID or NUCLEOTIDE)
        
    Returns:
        List[str]: List of valid symbols for the mutation type
        
    Raises:
        ValueError: If the mutation type is unknown
    """
    if mutation_type == MutationType.AMINO_ACID:
        return ["A", "C", "D", "E", "F", "G", "H", "I", "K", 
                "L", "M", "N", "P", "Q", "R", "S", "T", 
                "V", "W", "Y"]
    elif mutation_type == MutationType.NUCLEOTIDE:
        return ["A", "T", "C", "G"]
    else:
        raise ValueError(f"Unknown mutation type: {mutation_type}")


def extract_position(mutation_str: str) -> int:
    """Extract the genomic position number from a mutation string.
    
    This function works for both nucleotide and amino acid mutations:
    - Nucleotide mutations: "C345T", "456-", "748G" -> returns position number
    - Amino acid mutations: "ORF1a:T103L", "S:N126K" -> returns position number
    
    Args:
        mutation_str (str): The mutation string to parse
        
    Returns:
        int: The genomic position number, or 0 if parsing fails
        
    Examples:
        >>> extract_position("C345T")
        345
        >>> extract_position("ORF1a:T103L")
        103
        >>> extract_position("456-")
        456
        >>> extract_position("S:N126K")
        126
    """
    mutation_str = mutation_str.strip()
    
    # Check if it's an amino acid mutation (contains ':')
    if ':' in mutation_str:
        # Amino acid mutation format: "ORF1a:T103L", "S:N126K"
        try:
            # Split by colon and get the mutation part
            gene_part, mutation_part = mutation_str.split(':', 1)
            
            # Extract position from the mutation part
            # Pattern: [REF][POSITION][ALT] where REF and ALT are amino acids
            amino_acids = get_symbols_for_mutation_type(MutationType.AMINO_ACID)
            amino_acid_pattern = '|'.join(amino_acids)
            
            # Match patterns like T103L, N126K, etc.
            match = re.match(rf"^({amino_acid_pattern})?(\d+)({amino_acid_pattern}|-)?$", mutation_part.upper())
            if match:
                return int(match.group(2))
                
        except (ValueError, IndexError):
            pass
    else:
        # Nucleotide mutation format: "C345T", "456-", "748G"
        try:
            nucleotides = get_symbols_for_mutation_type(MutationType.NUCLEOTIDE)
            nucleotide_pattern = '|'.join(nucleotides)
            
            # Match patterns like C345T, 456-, 748G
            match = re.match(rf"^({nucleotide_pattern})?(\d+)({nucleotide_pattern}|-)?$", mutation_str.upper())
            if match:
                return int(match.group(2))
                
        except ValueError:
            pass
    
    # Fallback: try to extract any number from the string
    try:
        numbers = re.findall(r'\d+', mutation_str)
        if numbers:
            return int(numbers[0])
    except (ValueError, IndexError):
        pass
    
    return 0  # Return 0 if position extraction fails

def lapis_mutation_to_pos(mutation: str) -> Optional[str]:
    """
    Convert a LAPIS nucleotide mutation string to tallymut pos format.

    LAPIS returns mutations as "C241T" (ref + position + alt).
    Tallymut format used by LolliPop is "241T" (position + alt, no ref).

    Args:
        mutation: LAPIS mutation string e.g. "C241T", "T508-"

    Returns:
        str: tallymut pos format e.g. "241T", "508-"
        None: if the mutation string is malformed or an insertion
    """
    if not mutation or len(mutation) < 3:
        return None
    match = re.match(r'^[A-Z](\d+)([A-Z\-])$', mutation)
    if not match:
        return None
    position = match.group(1)
    alt = match.group(2)
    if alt == '-':
        return f"{position}-"
    return f"{position}{alt}"

def sort_mutations_by_position(mutations: List[str]) -> List[str]:
    """Sort a list of mutations by their genomic position in ascending order.
    
    Args:
        mutations (List[str]): List of mutation strings to sort
        
    Returns:
        List[str]: Sorted list of mutations by genomic position
        
    Examples:
        >>> mutations = ["C500T", "A100G", "T300C"]
        >>> sort_mutations_by_position(mutations)
        ["A100G", "T300C", "C500T"]
    """
    return sorted(mutations, key=extract_position)


def validate_mutation(mutation_str: str, mutation_type: MutationType) -> bool:
    """Validate if a mutation string conforms to the expected format for the given mutation type.
    
    Args:
        mutation_str (str): The mutation string to validate
        mutation_type (MutationType): The type of mutation (AMINO_ACID or NUCLEOTIDE)
        
    Returns:
        bool: True if the mutation string is valid, False otherwise
        
    Examples:
        >>> validate_mutation("C345T", MutationType.NUCLEOTIDE)
        True
        >>> validate_mutation("ORF1a:T103L", MutationType.AMINO_ACID)
        True
        >>> validate_mutation("X123Y", MutationType.AMINO_ACID)
        False
    """
    mutation_str = mutation_str.strip()
    
    if mutation_type == MutationType.AMINO_ACID:
        # Amino acid mutation format: "ORF1a:T103L", "S:N126K"
        if ':' not in mutation_str:
            return False
        try:
            gene_part, mutation_part = mutation_str.split(':', 1)
            amino_acids = get_symbols_for_mutation_type(MutationType.AMINO_ACID)
            amino_acid_pattern = '|'.join(amino_acids)
            match = re.match(rf"^({amino_acid_pattern})?(\d+)({amino_acid_pattern}|-)?$", mutation_part.upper())
            return match is not None
        except (ValueError, IndexError):
            return False
            
    elif mutation_type == MutationType.NUCLEOTIDE:
        # Nucleotide mutation format: "C345T", "456-", "748G"
        try:
            nucleotides = get_symbols_for_mutation_type(MutationType.NUCLEOTIDE)
            nucleotide_pattern = '|'.join(nucleotides)
            match = re.match(rf"^({nucleotide_pattern})?(\d+)({nucleotide_pattern}|-)?$", mutation_str.upper())
            return match is not None
        except ValueError:
            return False
            
    return False  # Unknown mutation type


def possible_mutations_at_position(position: int, mutation_type: MutationType, gene: Optional[str] = None, include_reference: bool = True) -> List[str]:
    """Generate all possible mutations at a given genomic position for the specified mutation type.
    
    Args:
        position (int): The genomic position number
        mutation_type (MutationType): The type of mutation (AMINO_ACID or NUCLEOTIDE)
        gene (str, optional): The gene name for amino acid mutations (e.g., "ORF1ab", "S")
        include_reference (bool): Whether to include reference base in mutation strings.
                                If False, generates simpler mutations like "100T" instead of "A100T"
        
    Returns:
        List[str]: List of possible mutation strings at the given position
        
    Examples:
        >>> possible_mutations_at_position(100, MutationType.NUCLEOTIDE, include_reference=True)
        ["A100T", "A100C", "A100G", "T100A", "T100C", "T100G", ...]
        
        >>> possible_mutations_at_position(100, MutationType.NUCLEOTIDE, include_reference=False)
        ["100A", "100T", "100C", "100G", "100-"]
        
        >>> possible_mutations_at_position(50, MutationType.AMINO_ACID, "ORF1ab", include_reference=False)
        ["ORF1ab:50A", "ORF1ab:50C", "ORF1ab:50D", ..., "ORF1ab:50Y"]
    """
    symbols = get_symbols_for_mutation_type(mutation_type)
    mutations = []
    
    if include_reference:
        # Original behavior: generate all ref->alt combinations
        for ref in symbols:
            for alt in symbols:
                if ref != alt:
                    if mutation_type == MutationType.AMINO_ACID and gene:
                        mutations.append(f"{gene}:{ref}{position}{alt}")
                    else:
                        mutations.append(f"{ref}{position}{alt}")
    else:
        # Simplified behavior: generate position->alt mutations (no reference base)
        for alt in symbols:
            if mutation_type == MutationType.AMINO_ACID and gene:
                mutations.append(f"{gene}:{position}{alt}")
            else:
                mutations.append(f"{position}{alt}")
        
        # Add deletion mutation
        if mutation_type == MutationType.AMINO_ACID and gene:
            mutations.append(f"{gene}:{position}-")
        else:
            mutations.append(f"{position}-")
    
    return mutations