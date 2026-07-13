# SensCalc arXiv:2501.04525 decay inputs

This directory is the reviewable export of the universally fermion-coupled
ALP decay data used by GRENDEL BC10. The source is the official SensCalc
release `v.1.3.3`, pinned at commit
`0bca050633aae16e148d47f21840fa07ff4b8724`.

Upstream repository: <https://github.com/maksymovchynnikov/SensCalc>

SensCalc methodology: arXiv:2305.13383.

Decay phenomenology: arXiv:2501.04525. This description is invariant under
the chiral-rotation parameter choice and includes mixing with heavy
pseudoscalar excitations. SensCalc `v.1.3.3` is used because its release notes
specifically record a fix to the decay-mode description of universally
fermion-coupled ALPs.

## Source files

The binary Wolfram MX dumps are not copied into this repository. The exporter
checks these exact SHA-256 values before decoding them:

| Input | SHA-256 |
|---|---|
| `Widths-model-ALP-fermion-scale-1000.-GeV-2501.04525.m` | `d12fb78d28edff0aa081c9fb66d829b42b4ec71202684019d7b9047ecb40b869` |
| `Br-ratios-SensCalc-model-ALP-fermion-scale-1000.-GeV-2501.04525.m` | `34a09eed87d081bfffe79b464094741454022f79478c5b28bc236bc361049286` |
| `Matrix-elements-squared-model-ALP-fermion-scale-1000.-GeV-2501.04525.m` | `f959c2257fa349e5af6966795db8cbf0da2e3ff8d097b5fccc3316d270ff17ed` |

All three inputs are required. Widths determine the lifetime, branching ratios
determine the visible-channel mixture, and the squared matrix elements
determine the charged-daughter kinematics used by the reconstruction
acceptance.

## Regeneration

Obtain and activate the free Wolfram Engine, then run:

```bash
git clone --branch v.1.3.3 --depth 1 \
  https://github.com/maksymovchynnikov/SensCalc /tmp/SensCalc-v1.3.3
wolframscript -activate
python alp_fermion/tools/export_senscalc_2501.py /tmp/SensCalc-v1.3.3
```

On the Homebrew macOS installation, the wrapper automatically selects the
nested kernel at:

```text
/Applications/Wolfram Engine.app/Contents/Resources/Wolfram Player.app/Contents/MacOS/WolframKernel
```

Use `--check-only` to validate the SensCalc commit and source hashes without a
Wolfram license.

## Coupling conversion

SensCalc applies its raw width coefficient as

```text
Gamma = (g_Y / (2 v_h))^2 * coefficient_raw.
```

The GRENDEL BNT-axis convention is `1/f_BNT = g_Y/v_h`, so the exported BNT
coefficient is `coefficient_raw/4`. The branching ratios are convention
independent. Before the tables become the default model input, the conversion
must pass the physical anchors documented in `model.py`: the approximately
9% dimuon branching ratio at 1 GeV and the charm onset near 2.58 GeV.
