# fonll-nnpdf40

Local FONLL+LHAPDF regeneration of the public B-hadron and D0 meson
double-differential production cross sections at pp 14 TeV, with the
PDF set switched from CTEQ6.6 to `NNPDF40_nlo_as_01180`.

The public FONLL web form only ships predictions for CTEQ6.6. Downstream
HNL and dark-photon sensitivity studies still rely on those tables. This
repository contains everything needed to regenerate equivalent tables with
NNPDF4.0 NLO and reproduce the numerical agreement against the public
FONLL CTEQ6.6 points used as the validation reference.

## What is in this repository

    grids/    two pre-computed dsigma/dpT/dy tables, NNPDF4.0 NLO,
              pp 14 TeV, central scale, fragmentation fraction 1
    patches/  drop-in replacements for two files of FONLL v1.3.3,
              adding BCFY fragmentation modes 4 and 5/8 and raising
              the rapidity-grid capacity from 100 to 120 nodes
    scripts/  Python drivers that run the patched FONLL executables,
              combine direct BCFY pseudoscalar and calibrated D* feeddown
              for the charm grid, and validate against public FONLL
              CTEQ6.6 reference points

What is *not* in this repository: the FONLL v1.3.3 source tree itself, an
LHAPDF installation, the NNPDF4.0 PDF set, and any Python virtualenv. These
must be provided locally; see [Setup](#setup) below.

## Central setup

| Quantity | Value |
| --- | --- |
| Process | pp, sqrt(s) = 14 TeV, ebeam1 = ebeam2 = 7000 GeV |
| Calculation | FONLL v1.3.3, NLO + NLL |
| PDF | `NNPDF40_nlo_as_01180`, LHAPDF id 331700 |
| Heavy-quark masses | m_b = 4.75 GeV, m_c = 1.50 GeV |
| Scales | mu_R = mu_F = sqrt(m^2 + pT^2), central |
| Grid | 100 pT in [0, 50] GeV, 100 y in [-3, 3] |
| Output | pT y dsigma/dpT/dy in pb/GeV |
| Bottom fragmentation | Kartvelishvili `(1-z) z^alpha`, alpha = 24.2 (public FONLL B-hadron default with `n5moment`) |
| Charm fragmentation | BCFY pseudoscalar with r = 0.1 plus calibrated D* -> D0 feeddown, matching the public FONLL D0 convention |
| Fragmentation frame | `ifrframe = 1` (y=0 frame, public FONLL default for `dsigma/dpT/dy`) |

## Validation

Each grid is validated point-by-point against the public FONLL web form
(CTEQ6.6, FONLL v1.3.2, 14 TeV), queried 2026-05-29. The validation runs
*also* through the locally patched code with CTEQ6.6 selected as the PDF,
so the comparison isolates the PDF swap from any local change.

| Channel | local | public FONLL CTEQ6.6 | relative |
| --- | --- | --- | --- |
| B hadron, pT=5 GeV, y=0 | 7.535277e6 pb/GeV | 7.5329e6 pb/GeV | +0.032% |
| D* vector BCFY r=0.1, pT=5 GeV, y=0 | 4.0129e7 pb/GeV | 4.0129e7 pb/GeV | match |
| D0 direct BCFY pseudoscalar r=0.1, pT=5 GeV, y=0 | 3.679271e7 pb/GeV | 3.5221e7 pb/GeV | +4.46% (no feeddown; diagnostic only) |
| D0 with calibrated D* -> D0 feeddown, 9 points in 1..50 GeV at y=0 | -- | -- | max |delta| / public = 1.94e-3 |

The 4.46% offset in the direct BCFY pseudoscalar row is what is corrected
by the `D* -> D0` feeddown convolution used in the published charm grid.
The feeddown weight (`charm_feeddown_vector_to_direct_weight`) is fitted
against the nine public D0 points at y=0 listed inside
`scripts/generate_meson_grids.py`; the resulting maximum residual is below
2 per-mille.

## Setup

Requires:

1. Python 3.10 or newer with `numpy`. No other Python deps are needed for
   the grid generation.
2. A working `gfortran` (the FONLL Makefile detects Linux and Darwin) and
   a recent LHAPDF (>= 6.x) installation visible via `lhapdf-config`.
3. The FONLL v1.3.3 source tree. Obtain the tarball from
   <http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllinput.html>
   and unpack it to `src/fonll/` inside this repository so that the layout
   becomes `src/fonll/{misc1,common,main,hdmassive,...}`.
4. The NNPDF4.0 NLO PDF set:
   `lhapdf install NNPDF40_nlo_as_01180`. For the CTEQ6.6 validation runs
   also install `cteq66`.

Apply the patches by overwriting the two files in place:

    cp patches/fragmfonll.f src/fonll/misc1/fragmfonll.f
    cp patches/fonllgrid.f  src/fonll/misc1/fonllgrid.f

Build the two LHAPDF-linked executables from `src/fonll/Linux/`:

    cd src/fonll/Linux
    make -f ../misc1/Makefile \
        VPATH=../misc1:../main:../hdmassive:../hdresummed:../phmassive:../phresummed:../common \
        fonllgridlha fragmfonll

## Regenerating the grids

After building, from the repository root:

    python scripts/generate_meson_grids.py --pdf nlo --quark bottom --grid-workers 4
    python scripts/generate_meson_grids.py --pdf nlo --quark charm  --grid-workers 4

Output `.dat` files are written under `output/` (created on first run).
The two pre-computed grids in `grids/` were produced this way.

`--grid-workers N` splits the rapidity axis into `N` contiguous chunks
that are computed by independent FONLL processes and re-stitched before
fragmentation. There is no shared state across chunks; results are
deterministic up to floating-point reordering at chunk boundaries (none
in practice for the grids here).

## Validating against public FONLL CTEQ6.6

    python scripts/validate_cteq66_public_points.py
    python scripts/validate_charm_bcfy_point.py

The first script regenerates a CTEQ6.6 grid with this codebase and
compares the (pT, y=0) row to the public web form values listed at the
top of the script. The second pulls the standalone direct-BCFY charm
point used in the +4.46% diagnostic.

## Layout of a generated `.dat` file

The output files are ASCII, three columns:

    pT [GeV]    y    dsigma/dpT/dy [pb/GeV]

100 pT nodes from 0 to 50 GeV crossed with 100 rapidity nodes from -3 to
+3, giving 10 000 rows. A short header (lines beginning with `#`) records
the PDF set, fragmentation parameters, FONLL version, heavy-quark mass,
and the meson-mass convention.

## Citing

If you use these grids please cite the FONLL papers (the same ones used
by the public web form):

- Cacciari, Greco, Nason, JHEP 9805 (1998) 007, hep-ph/9803400
- Cacciari, Frixione, Nason, JHEP 0103 (2001) 006, hep-ph/0102134

and the NNPDF4.0 release:

- NNPDF Collaboration, Eur. Phys. J. C 82 (2022) 428, arXiv:2109.02653

A `CITATION.cff` is provided for the repository itself.

## License

MIT for everything in `scripts/`, `patches/` (the modifications), this
README, and the auxiliary metadata. The unmodified portions of
`patches/*.f` remain the work of M. Cacciari, S. Frixione and P. Nason;
the patches are distributed under the assumption that users obtain FONLL
through the official distribution and that they keep the original
attribution intact. The generated `.dat` files in `grids/` are released
under CC0.
