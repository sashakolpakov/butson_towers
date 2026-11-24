"""Sage-free KKPR-style tools for root-of-unity Hadamard matrices."""

from .coefficients.cyclotomic import CyclotomicElement, CyclotomicField
from .hadamard.family import FamilyVerification, MonomialHadamardFamily
from .hadamard.group_developed import FiniteGroupTable, GroupDevelopedTorus, cyclic_group
from .hadamard.model import HadamardTorus
from .laurent.ideal import LaurentIdeal
from .laurent.polynomial import LaurentPolynomial, LaurentRing

__all__ = [
    "CyclotomicElement",
    "CyclotomicField",
    "FamilyVerification",
    "FiniteGroupTable",
    "GroupDevelopedTorus",
    "HadamardTorus",
    "MonomialHadamardFamily",
    "LaurentIdeal",
    "LaurentPolynomial",
    "LaurentRing",
    "cyclic_group",
]

__version__ = "0.2.0"
