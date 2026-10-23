"""Turn catalogue date strings into (year, precision) pairs.

Rules from PLAN.md: ranges become their midpoint, and anything vaguer than 25 years
is rejected. Precision is the half-width of the interval the string allows.
"""
import re

MAX_PRECISION = 25

_IE = re.compile(r"i\.\s?e\.?,?\s*\[?\s*(1[4-9]\d\d)")
_RANGE = re.compile(r"(1[4-9]\d\d)\s*(?:-|\u2013|\u2014|to|/|bis|tot|and|&)\s*(1[4-9]\d\d|\d\d)\b")  # hyphen, en and em dash
_DECADE = re.compile(r"\b(1[4-9]\d)[-?_u]\b|\b(1[4-9]\d)[-?_u](?=\W|$)")
_CENTURY = re.compile(r"\b(1[4-9])(?:--|__|uu|\?\?)")
_YEAR = re.compile(r"(?<!\d)(1[4-9]\d\d)(?!\d)")


def parse_year(text):
    """Return (year, precision) or None."""
    if not text:
        return None
    s = str(text)
    m = _IE.search(s)
    if m:
        return float(m.group(1)), 0.0
    m = _RANGE.search(s)
    if m:
        a = int(m.group(1))
        b = m.group(2)
        b = int(b) if len(b) == 4 else int(str(a)[:2] + b)
        if b >= a:
            half = (b - a) / 2
            return (a + half, half) if half <= MAX_PRECISION else None
    m = _YEAR.search(s)
    if m:
        return float(m.group(1)), 0.0
    m = _DECADE.search(s)
    if m:
        d = int(m.group(1) or m.group(2))
        return d * 10 + 5.0, 5.0
    if _CENTURY.search(s):
        return None  # century only: too vague
    return None
