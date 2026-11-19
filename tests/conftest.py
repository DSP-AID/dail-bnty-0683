import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from seestar_eff import pipeline


@pytest.fixture(scope="session", autouse=True)
def _data():
    """Ship data as tarball only; extract members before any test reads them."""
    pipeline.ensure_data(Path(__file__).resolve().parents[1] / "data")
