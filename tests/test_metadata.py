import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from analyze import great_circle_km, region  # noqa: E402
from metadata import FACSIMILE, geocode, load_gazetteer, norm_maker  # noqa: E402

GAZ = load_gazetteer()


@pytest.mark.parametrize("place, lat, lon", [
    ("Amstelodami", 52.37, 4.89), ("Lugduni Batavorum", 52.16, 4.49), ("Paris", 48.85, 2.35),
    ("Venice", 45.44, 12.33), ("'s-Gravenhage", 52.08, 4.30), ("Edinburgh,London; Dublin", 55.95, -3.19),
    ("Hildburghausen", 50.43, 10.73), ("New York", 40.71, -74.01),
])
def test_geocode_known_imprints(place, lat, lon):
    got = geocode(place, GAZ)
    assert got is not None
    assert great_circle_km(got[0], got[1], lat, lon) < 30


def test_geocode_rejects_sine_loco():
    assert geocode("[S.l.]", GAZ) is None


def test_maker_normalisation_strips_dates():
    assert norm_maker("Blaeu, Willem Janszoon, 1571-1638") == "blaeu willem janszoon"
    assert norm_maker(["Mitchell, Samuel Augustus, 1792-1868", "Other"]) == "mitchell samuel augustus"


def test_facsimile_flag():
    assert FACSIMILE.search("Facsimile-atlas to the early history of cartography")
    assert not FACSIMILE.search("Nouvel atlas de la France")


@pytest.mark.parametrize("lon, lat, name", [
    (2, 48, "Europe"), (-75, 40, "North America"), (-60, -20, "Latin America"),
    (20, 0, "Africa & Middle East"), (100, 30, "Asia"), (140, -25, "Oceania & Pacific"),
])
def test_regions(lon, lat, name):
    assert region(lon, lat) == name
