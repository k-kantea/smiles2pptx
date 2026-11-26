"""2D coordinates of the atoms (the third stage): layout(mol) -> one point per atom.

The bond length is 1 and y points up. Ring systems are placed in rings.py (bridged ones in
springs.py) and everything else in place.py; untangle.py clears overlaps, and lengthen.py
lengthens bonds that would look too short between their labels.
"""
from .geometry import Point
from .place import layout

__all__ = ["Point", "layout"]
