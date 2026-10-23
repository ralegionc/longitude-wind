import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from dates import parse_year  # noqa: E402


@pytest.mark.parametrize("text, expected", [
    ("1862", (1862, 0)),
    ("ca. 1700", (1700, 0)),
    ("[1700?]", (1700, 0)),
    ("c1890", (1890, 0)),
    ("[1816-1823]", (1819.5, 3.5)),
    ("1816-23", (1819.5, 3.5)),
    ("1782 [i.e. 1790]", (1790, 0)),
    ("185-?", (1855, 5)),
    ("[185-]", (1855, 5)),
    ("Amsterdam, 1650.", (1650, 0)),
])
def test_parses(text, expected):
    y, p = parse_year(text)
    assert y == pytest.approx(expected[0]) and p == pytest.approx(expected[1])


@pytest.mark.parametrize("text", ["", None, "[17--]", "18??", "1700-1800", "undated", "2004"])
def test_rejects_vague_or_out_of_range(text):
    assert parse_year(text) is None
