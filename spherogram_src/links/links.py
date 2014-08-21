from .links_base import Crossing, Strand, CrossingStrand, CrossingEntryPoint
try:
    from .invariants import Link   # Only works in Sage
except ImportError:
    from .links_base import Link   # Lacks some methods for computing invariants

