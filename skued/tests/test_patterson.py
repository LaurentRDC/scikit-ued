# -*- coding: utf-8 -*-


import numpy as np

from crystals import Crystal
from skued import patterson, powdersim


def test_patterson_output_shape():
    """Test that the output shape is as expected."""
    # Simulate a powder pattern first
    crystal = Crystal.from_database("vo2-m1")
    q = np.linspace(0.2, 10, 1024)
    I = powdersim(crystal=crystal, q=q)

    radii = np.arange(0.1, 5, 1 / 50)
    pairdist = patterson(q=q, I=I, crystal=crystal, radii=radii)

    assert radii.shape == pairdist.shape
