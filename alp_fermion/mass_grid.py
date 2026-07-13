"""Shared BC10 ALP mass grid with the unsupported pole points retained.

The sensitivity output needs marker rows at the eta and eta-prime poles so
plotting cannot bridge them. Generators inspect ``model`` (or the standalone
window constants in the ROOT-only template generator) and skip those points.
"""

ALP_MASS_GRID = sorted({round(x, 2) for x in (
    [0.22 + 0.02 * i for i in range(40)]
    + [1.00 + 0.05 * i for i in range(41)]
    + [3.00 + 0.10 * i for i in range(17)]
    + [4.65, 4.70, 4.75]
)})

