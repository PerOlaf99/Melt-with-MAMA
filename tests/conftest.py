"""Shared pytest setup for the varmelt test suite."""
import os
import sys

import pytest

# Import varmelt from this checkout, wherever it lives.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from varmelt import inpcr  # noqa: E402


@pytest.fixture(autouse=True)
def _fresh_inpcr_cache():
    """The in-silico PCR scan cache is keyed on (assembly, chrom, ...), so
    tests that stub different synthetic chromosomes under the same name would
    otherwise see the first test's sequence."""
    inpcr.clear_cache()
    yield
    inpcr.clear_cache()
