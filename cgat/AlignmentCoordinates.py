"""AlignmentCoordinates.py - Pure-Python coordinate mapping for alignments
==========================================================================

A lightweight replacement for the unmaintained ``alignlib_lite`` package.
Only the coordinate-mapping surface used by cgat-apps is implemented:
building maps from pairs/diagonals/blocks, querying mapped positions, and
comparing maps (identity, overlap, shortest distance).

No dynamic-programming sequence alignment is performed.
"""

from __future__ import annotations

import sys


class Alignment:
    """Map row coordinates to column coordinates.

    Internally stores ``row -> col``. Adding the same row with a different
    column marks that row as ambiguous; :meth:`mapRowToCol` then raises
    :class:`RuntimeError` if any ambiguity exists in the map (matching
    alignlib's uniqueness check used by ``diff_chains.validateChain``).
    """

    def __init__(self):
        self._map = {}
        self._ambiguous = set()
        self._row_from = None
        self._row_to = None
        self._col_from = None
        self._col_to = None

    def clear(self):
        self._map.clear()
        self._ambiguous.clear()
        self._row_from = None
        self._row_to = None
        self._col_from = None
        self._col_to = None

    def getNew(self):
        """Return a new empty alignment of the same type."""
        return Alignment()

    def copy(self, other):
        """Copy *other* into this alignment (replacing current contents)."""
        self.clear()
        self._map = dict(other._map)
        self._ambiguous = set(other._ambiguous)
        self._row_from = other._row_from
        self._row_to = other._row_to
        self._col_from = other._col_from
        self._col_to = other._col_to
        return self

    def _update_bounds(self, row, col):
        if self._row_from is None or row < self._row_from:
            self._row_from = row
        if self._row_to is None or row > self._row_to:
            self._row_to = row
        if self._col_from is None or col < self._col_from:
            self._col_from = col
        if self._col_to is None or col > self._col_to:
            self._col_to = col

    def addPair(self, row, col):
        """Add a mapped pair ``(row, col)``."""
        if row in self._ambiguous:
            return
        if row in self._map:
            if self._map[row] != col:
                del self._map[row]
                self._ambiguous.add(row)
            return
        self._map[row] = col
        self._update_bounds(row, col)

    def addPairExplicit(self, row, col, score=0):
        """Add a mapped pair; *score* is accepted for API compatibility."""
        self.addPair(row, col)

    def addDiagonal(self, row_from, row_to, offset):
        """Add pairs ``(r, r + offset)`` for ``r`` in ``[row_from, row_to)``."""
        for row in range(row_from, row_to):
            self.addPair(row, row + offset)

    def addBlocks(self, query_starts, sbjct_starts, sizes):
        """Add contiguous blocks given parallel start/size lists."""
        for q, s, size in zip(query_starts, sbjct_starts, sizes):
            self.addDiagonal(q, q + size, s - q)

    def mapRowToCol(self, row):
        """Return the column for *row*, or ``-1`` if unmapped.

        Raises
        ------
        RuntimeError
            If any row in this alignment maps to conflicting columns.
        """
        if self._ambiguous:
            raise RuntimeError(
                "alignment has ambiguous row mappings: %s" %
                sorted(self._ambiguous)[:5])
        return self._map.get(row, -1)

    def getNumAligned(self):
        """Return the number of unambiguously aligned residue pairs."""
        return len(self._map)

    def getRowFrom(self):
        """Return the first (minimum) aligned row, or 0 if empty."""
        return 0 if self._row_from is None else self._row_from

    def getRowTo(self):
        """Return one past the last aligned row (exclusive end), or 0 if empty.

        Matches alignlib's ``getRowTo()`` convention used by gff2psl/Blat.
        """
        return 0 if self._row_to is None else self._row_to + 1

    def getColFrom(self):
        """Return the first (minimum) aligned column, or 0 if empty."""
        return 0 if self._col_from is None else self._col_from

    def getColTo(self):
        """Return one past the last aligned column (exclusive end), or 0 if empty."""
        return 0 if self._col_to is None else self._col_to + 1

    def getBlocks(self):
        """Return contiguous diagonal blocks as ``(row_start, col_start, size)``.

        Blocks are sorted by increasing row start.
        """
        if not self._map:
            return []

        items = sorted(self._map.items())
        blocks = []
        row0, col0 = items[0]
        size = 1
        prev_row, prev_col = row0, col0

        for row, col in items[1:]:
            if row == prev_row + 1 and col == prev_col + 1:
                size += 1
            else:
                blocks.append((row0, col0, size))
                row0, col0 = row, col
                size = 1
            prev_row, prev_col = row, col

        blocks.append((row0, col0, size))
        return blocks

    def alignedRows(self):
        """Return a frozenset of unambiguously aligned rows."""
        return frozenset(self._map)

    def __len__(self):
        return len(self._map)

    def __bool__(self):
        return bool(self._map)


def makeAlignmentBlocks():
    """Return a new empty :class:`Alignment` (alignlib-compatible name)."""
    return Alignment()


def makeAlignmentVector():
    """Return a new empty :class:`Alignment` (alignlib-compatible name)."""
    return Alignment()


def addDiagonal2Alignment(aln, row_from, row_to, offset):
    """Add a diagonal to *aln* (module-level alias)."""
    aln.addDiagonal(row_from, row_to, offset)


def getAlignmentIdentity(a, b):
    """Count rows aligned in both maps to the same column (row-row)."""
    n = 0
    # Iterate the smaller map for speed.
    if len(a._map) > len(b._map):
        a, b = b, a
    for row, col in a._map.items():
        if b._map.get(row, sys.maxsize) == col:
            n += 1
    return n


def getAlignmentOverlap(a, b):
    """Count rows that are aligned in both maps (row-row), ignoring columns."""
    return len(a._map.keys() & b._map.keys())


def getAlignmentShortestDistance(a, b):
    """Minimum gap between the aligned row ranges of *a* and *b*.

    Returns 0 if the row ranges overlap or either alignment is empty.
    """
    if not a._map or not b._map:
        return 0

    a_from, a_to = a.getRowFrom(), a.getRowTo()
    b_from, b_to = b.getRowFrom(), b.getRowTo()

    # getRowTo() is exclusive, so adjacent ranges [0,10) and [10,15) touch.
    if a_to <= b_from:
        return b_from - a_to
    if b_to <= a_from:
        return a_from - b_to
    return 0


def formatEmissions(aln):
    """Return a simple emissions-style string for an alignment.

    Format: alternating run lengths of aligned pairs as ``+N+N`` segments,
    sufficient for ``diff_fasta --correct-shift`` diagnostic output.
    """
    blocks = aln.getBlocks()
    if not blocks:
        return ""

    parts = []
    prev_row_end = None
    prev_col_end = None
    for row, col, size in blocks:
        if prev_row_end is not None:
            row_gap = row - prev_row_end
            col_gap = col - prev_col_end
            if row_gap:
                parts.append("-%i" % row_gap)
            if col_gap:
                parts.append("+%i" % col_gap)
        parts.append("+%i+%i" % (size, size))
        prev_row_end = row + size
        prev_col_end = col + size
    return "".join(parts)
