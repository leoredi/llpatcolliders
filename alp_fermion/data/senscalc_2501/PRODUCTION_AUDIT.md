# BC10 production audit against arXiv:2501.04525

_Audit date: 2026-07-13._ This is an implementation audit for the GRENDEL
geometry, not a claim that every production mode discussed in the paper must
be included in every experiment. GRENDEL is a transverse detector covering
approximately `|eta| < 0.5`; both the total production probability and the
accepted angular/energy distribution matter.

## Pinned sources

- SensCalc `v.1.3.3`, commit
  `0bca050633aae16e148d47f21840fa07ff4b8724`.
- Official arXiv:2501.04525 companion repository
  <https://github.com/maksymovchynnikov/ALPs-phenomenology>, commit
  `e1443ac1a07ed8a702173e7ead97aa231eedb563` (2026-02-03).
- GKOZ arXiv:2310.03524 for the LHC channel hierarchy and the exclusive
  `B -> K^(i) a` calculation already used by GRENDEL.

The relevant SensCalc source hashes are:

| File | SHA-256 |
|---|---|
| `3. ALP-fermion sensitivity.nb` | `48883fb7ecce956bd73965709ecf44d0cae286e637067c672c1566d736072ef3` |
| `codes/generic.nb` | `95cc528268e41b273182377275bf9eeb00d8cc1cb92d149ddeaa9e72967f0d4a` |
| `codes/LLP distribution/prod-pheno-ALP-fermion.nb` | `e1a7b91a4e82421081d165561a1b7f23e90aefe451efca32a5a3459b1ac0d0f7` |
| `codes/Acceptances/ALP-fermion.nb` | `8d205b65da456fb1a8fac8d1803eaa73ae839128a60d0624eeea17cacee5d4b4` |
| `codes/experiments.nb` | `3e60809105f49d8c5274de9a352cd0063bac9cb8feba31188e5e49ca6e468bdf` |
| `Coefficients-ALP-fermion-Lambda=1.-TeV.mx` | `e6885569ee8d33fb5310723ba2f020bf8b4e3506c2b5ad258b5230f9aa981f43` |
| `ProductionProbability-Fragmentation-ALP-fermion.mx` | `4804205cfb14c0d40322fcabb7ec65843023ee44e91f6918603bca0e8835e9eb` |
| `DoubleDistr_ALP-fermion_Fragmentation_LHC.m` | `43116651d526358c446af54a16d43ca9b2233fb42ff7ab4bae6946fe21477a79` |
| `BrRatios-Msquared-LightMesonDecays-ALP-fermion-Lambda=1.-TeV.mx` | `c7468dcf8863a8ac182c2720d34d12b2dd12693700aec032460cd5d678985f71` |
| `ProdProb_ALP-fermion_Bremsstrahlung-AP_LHC.m` | `435022c25a2bb908ffd5c4b12a5645a0d5aa7f6914c49aab5e33fe0b8b60e3d2` |
| `sigmaDrellYan_ALP-fermion_LHC.txt` | `12e013f6df9a65c07514e75d8f1c47d53b53ebc72b21624ab13eb70afaed185d` |

## Channel decisions

| Production mode | GRENDEL status | Decision and evidence |
|---|---|---|
| Exclusive `B` decays | Included | The current code samples `B+` and `B0` from the 14 TeV FONLL spectrum and sums the nine-state kaon tower. This is the GKOZ channel set retained by the 2501 implementation and is more differential than SensCalc's merged `Bcharged -> PiCharged + a` event topology. Keep it. |
| Drell-Yan / gluon fusion | Not included | GKOZ's LHC production plot places this channel roughly five to six orders below the `B` tower over the current GRENDEL island. The tabulation starts at 1.5 GeV, while the present island closes near 2.7 GeV and `B` production remains open to about 4.8 GeV. It cannot affect the published island at current precision; document the omission. |
| Old flux-times-mixing | Deliberately excluded | ArXiv:2501.04525 identifies this approximation as chiral-rotation dependent and kinematically ambiguous. SensCalc still exposes `Old-Mixing-Pi0/Eta/EtaPr`; they must not be enabled as a shortcut. |
| Proton bremsstrahlung | Not included | The 2501 calculation is a forward quasi-real approximation with a large theory uncertainty. For a transverse detector it requires a dedicated angular acceptance calculation. It is not justified to add a total-rate reweighting to the FONLL `B` sample. |
| Quark fragmentation | Audit pending | The generalized mixing construction is the correct 2501 treatment below about 2 GeV. SensCalc's revised-production path exposes fragmentation and contains both probability and differential-spectrum assets. Decode them and evaluate the transverse acceptance before the final scan. |
| Light-meson decays | Audit pending at low mass | The 2501 modes are kinematically confined to the low-mass part of the scan (at most the parent-meson mass). Their very large parent flux means they cannot be dismissed from total rates alone. Decode the branching/matrix-element asset and test with LHC meson kinematics before finalizing masses below roughly 0.8 GeV. |
| `B_s -> phi a` | Not added | This is absent from GKOZ, ALPINIST, SensCalc's channel set, and the cited Boiarska kaon-tower form factors. It is an optional new calculation, not a missing implementation of arXiv:2501.04525. |

## Transverse parent-spectrum check

The five readable SensCalc LHC parent spectra are shape tables normalized by
the sampler; their physical multiplicities per collision are separate entries in
`codes/experiments.nb`. The source-pinned sampler in
`alp_fermion/production_spectra.py` gives the following parent-level fractions
inside `|eta| < 0.5`:

| Parent | Yield / collision | Parent fraction in `|eta| < 0.5` |
|---|---:|---:|
| eta | 3.64 | 0.0596 |
| eta-prime | 0.46 | 0.0653 |
| omega | 4.42 | 0.0555 |
| charged rho | 8.68181 | 0.0546 |
| K-short | 3.1 | 0.0606 |

These are not ALP acceptances: the daughter decay kinematics and the
mass-dependent branching ratios still have to be applied. The source tables
also use a 6.5 TeV beam energy (13 TeV collisions), whereas the GRENDEL FONLL
baseline is 14 TeV. They are suitable for deciding whether a dedicated 14 TeV
calculation is necessary, but not for silently mixing final event weights.

Reproduce the source check and parent-level summary with:

```bash
python alp_fermion/tools/export_senscalc_2501_production.py \
  --check-only /path/to/SensCalc
python -m alp_fermion.production_spectra /path/to/SensCalc
```

After Wolfram Engine is activated, omit `--check-only` to decode the binary
probability, matrix-element, bremsstrahlung, and fragmentation-grid assets into
the gitignored `alp_fermion/tmp/senscalc_2501_production/` directory.

SensCalc applies the raw meson-decay, fragmentation, and effective-production
coefficients as `(Fpi*g_Y/(2*v_h))^2 * coefficient_raw`, with
`Fpi = 0.093 GeV`. On the GRENDEL BNT axis, `1/f_BNT = g_Y/v_h`; consequently
the exported probability or branching coefficient is
`coefficient_bnt = Fpi^2 * coefficient_raw / 4`, and the physical probability
is `(1/f_BNT)^2 * coefficient_bnt`. The exporter writes only explicitly named
BNT coefficient tables; raw Drell-Yan cross-section coefficients remain marked
as raw and carry their separate normalization formula in metadata.

The pinned binary sources record `Windows-x86-64` in their MX headers. Wolfram
MX is system-dependent, so the source system ID is reported during preflight
and the decoder `$SystemID` is retained in the export manifest. A successful
cross-platform import must be demonstrated before any generated table is
adopted.

## SensCalc integration boundary

ArXiv:2305.13383 is the SensCalc methodology paper. SensCalc is beneficial
here as a pinned source of decay tables, branching ratios, matrix elements,
and production cross-checks. It should not replace the GRENDEL detector and
reconstruction simulation.

The `v.1.3.3` ALP-fermion sensitivity notebook exposes meson decays,
bremsstrahlung, Drell-Yan, and fragmentation when the revised `2501.04525`
production description is selected. Legacy neutral-meson mixing is retained
only as the alternate old description and must not be combined with the
revised channels. GRENDEL imports verified physics inputs and applies its own
geometry, decay sampling, stable-particle handling, and four-hit track
selection rather than substituting the SensCalc detector acceptance.

## Gate for the final 600k scan

Run the expensive production/template/sensitivity campaign once, after:

1. the 2501 decay widths, branching ratios, and squared matrix elements have
   been exported and validated together;
2. the low-mass fragmentation and light-meson production assets have been
   decoded and their daughter-level transverse accepted yields compared with
   the `B` tower;
3. any retained production mode has a real four-vector sample rather than a
   total-rate-only correction.
