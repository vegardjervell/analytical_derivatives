"""
The best docs. The very best.
"""
from __future__ import annotations
import collections.abc
import typing
__all__ = ['Analytical', 'PolyExp', 'PolyFrac', 'Polynomial', 'get_partitions', 'partition_multiplicity']
class Analytical:
    def derivative(self, arg0: typing.SupportsFloat, arg1: typing.SupportsInt) -> float:
        ...
    def numerical_derivative(self, arg0: typing.SupportsFloat, arg1: typing.SupportsInt, arg2: typing.SupportsFloat) -> float:
        ...

class PolyExp(Analytical):
    def __init__(self, arg0: Polynomial, arg1: Polynomial) -> None:
        ...

class PolyFrac(Analytical):
    def __init__(self, arg0: Polynomial, arg1: Polynomial) -> None:
        ...

class Polynomial(Analytical):
    def __init__(self, arg0: typing.SupportsInt, arg1: typing.SupportsInt, arg2: collections.abc.Sequence[typing.SupportsFloat], arg3: typing.SupportsInt) -> None:
        ...

def get_partitions(arg0: typing.SupportsInt) -> list[list[int]]:
    ...

def partition_multiplicity(arg0: collections.abc.Sequence[typing.SupportsInt]) -> float:
    ...
