# Patches against FONLL v1.3.3

These two files are drop-in replacements for the corresponding files
inside an unpacked FONLL v1.3.3 source tree:

    patches/fragmfonll.f  ->  src/fonll/misc1/fragmfonll.f
    patches/fonllgrid.f   ->  src/fonll/misc1/fonllgrid.f

The official FONLL distribution is available from
<http://www.lpthe.jussieu.fr/~cacciari/fonll/fonllinput.html>. Unpack the
tarball so that `misc1/`, `common/`, `main/`, etc. live at
`src/fonll/<name>/` relative to the root of this repository, then copy
the two files above over the originals.

## Summary of changes

### `fonllgrid.f`

One line modified:

    line 54:  parameter (nymax=120,nptmax=100)

(was `nymax=100`). This sets the size of the static rapidity arrays used
by the grid driver. A 100-node rapidity grid plus its terminating
sentinel needs `nymax >= 101`; the bump to 120 leaves room for slightly
denser grids without further changes.

No physics is altered. The arrays are simple storage and the loop bounds
elsewhere in the file are driven by the actual number of rapidity nodes
read from the input card, not by `nymax`.

### `fragmfonll.f`

Two kinds of changes.

(1) The rapidity array dimension `nymx` is raised from 100 to 120 at
each of the four places where it is declared in the file:

    line 246:  parameter(nymx=120,nptmx=250,maxfiles=10)
    line 672:  parameter(nymx=120,nptmx=250,maxfiles=10)
    line 741:  parameter (nymx=120,nptmx=250)
    line 1072: parameter(nymx=120,nptmx=250,maxfiles=10)

The four parameter declarations must agree for the COMMON-block layout
to be consistent across subroutines; that is why the change appears in
four places. As for `fonllgrid.f`, no physics formula is affected.

(2) Two new fragmentation options are added in
`fragmfonll.f::fragfun`, between the existing Kartvelishvili branch and
the trailing error stop. They implement the BCFY S-wave heavy-light
fragmentation functions of Braaten, Cheung, Fleming, Yuan,
Phys. Rev. D51 (1995) 4819 [hep-ph/9408231]:

    ifrag = 4  ->  BCFY vector (D*-like)
    ifrag = 5  ->  BCFY pseudoscalar (D0-like)
    ifrag = 8  ->  BCFY pseudoscalar (alias of 5, used by some downstream
                                       drivers that distinguish pseudoscalar
                                       species labels)

In all three cases the BCFY parameter is passed in through the existing
`ep` common block (renamed locally to `r` inside `fragfun`), which is
how the existing Kartvelishvili `alpha` is also passed. The two prompts
in `fragmfonll.f` that list the available fragmentation options were
updated to mention the new entries.

The added block is:

      elseif(ifrag.eq.4) then
    c BCFY vector S-wave heavy-light fragmentation model.
         r=ep
         den=1.d0-(1.d0-r)*z
         poly=2.d0-2.d0*(3.d0-2.d0*r)*z
        #        +3.d0*(3.d0-2.d0*r+4.d0*r**2)*z**2
        #        -2.d0*(1.d0-r)*(4.d0-r+2.d0*r**2)*z**3
        #        +(1.d0-r)**2*(3.d0-2.d0*r+2.d0*r**2)*z**4
         fragfun=xnorm*3.d0*r*z*(1.d0-z)**2*poly/den**6
      elseif(ifrag.eq.5.or.ifrag.eq.8) then
    c BCFY pseudoscalar S-wave heavy-light fragmentation model.
         r=ep
         den=1.d0-(1.d0-r)*z
         poly=6.d0-18.d0*(1.d0-2.d0*r)*z
        #        +(21.d0-74.d0*r+68.d0*r**2)*z**2
        #        -2.d0*(1.d0-r)*(6.d0-19.d0*r+18.d0*r**2)
        #        *z**3
        #        +3.d0*(1.d0-r)**2*(1.d0-2.d0*r+2.d0*r**2)
        #        *z**4
         fragfun=xnorm*r*z*(1.d0-z)**2*poly/den**6

(see `patches/fragmfonll.f` around line 595 for the full version in
context).

## License of the patched files

The unmodified portions of `fragmfonll.f` and `fonllgrid.f` are the work
of M. Cacciari, S. Frixione and P. Nason and ship under FONLL's
distribution terms. Only the modifications described above are released
under this repository's MIT license. The copies in this directory retain
the original file content verbatim apart from the changes listed here;
users are expected to obtain FONLL through the official channel and to
keep its original attribution and citation requirements intact.
