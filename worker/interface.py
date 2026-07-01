"""Shared inferface types for the app."""
from enum import Enum

# Define the mutation types as an Enum for better type safety
class MutationType(Enum):
    AMINO_ACID = "aminoAcid"
    NUCLEOTIDE = "nucleotide"
