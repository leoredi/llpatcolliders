"""
production/madgraph/_mg5_common.py

Shared plumbing for MG5-driven production pipelines (W/Z, prompt-tau).

Centralises the macOS LHAPDF / Mach-O dyld fix worked out in the W/Z driver
session so both drivers carry it. Three things this module provides:

1. MG5 executable resolution with a four-step fallback chain.
2. LHAPDF wiring: the resolver for lhapdf-config + the runtime libdir, an
   env-builder that injects DYLD_LIBRARY_PATH / LD_LIBRARY_PATH for any MG5
   subprocess, and a per-job patcher for me5_configuration.txt.
3. Mach-O rpath workaround: symlink libLHAPDF.dylib into every
   SubProcesses/P*_*/ directory so the compiled madevent/gensym binaries
   resolve @rpath/libLHAPDF.dylib via @loader_path.

The forced bin/madevent compile pre-build (also from the W/Z session) is
exposed as force_compile_subprocesses, so the rpath patch can run on real
binaries rather than empty source trees.
"""

import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
PYTHON_EXE = Path(sys.executable)


def _resolve_mg5_exe():
    """Resolve the mg5_aMC executable.

    Resolution order:
      1. $HNL_MG5_EXE (explicit path)
      2. hnl/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC (drop-in vendoring)
      3. projects-root shared vendored/MG5_aMC_v3_6_6/bin/mg5_aMC
         (siblings llpatcolliders_FONLL / llpatcolliders_MATT / ...
         share this install)
      4. sibling llpatcolliders_FONLL/vendored/MG5_aMC_v3_6_6/bin/mg5_aMC

    The 148 MB MG5 install is intentionally not committed; see
    hnl/vendored/PROVENANCE.md.
    """
    env_path = os.environ.get("HNL_MG5_EXE")
    if env_path:
        return Path(env_path)

    vendored = PROJECT_ROOT / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    if vendored.exists():
        return vendored

    shared = PROJECT_ROOT.parent.parent / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    if shared.exists():
        return shared

    rel = Path("llpatcolliders_FONLL") / "vendored" / "MG5_aMC_v3_6_6" / "bin" / "mg5_aMC"
    for base in (PROJECT_ROOT.parent, PROJECT_ROOT.parent.parent):
        candidate = base / rel
        if candidate.exists():
            return candidate
    return PROJECT_ROOT.parent.parent / rel


MG5_EXE = _resolve_mg5_exe()


# LHAPDF install used for NNPDF4.0 NLO. MG5 run_cards use `pdlabel = lhapdf`,
# so MadEvent needs lhapdf-config on its path; we point it at the shared
# install under projects-root NNPDF40/fonll-local/env/, which already carries
# NNPDF40_nlo_as_01180 (the same set the FONLL meson tables were generated
# with). Overridable via $HNL_LHAPDF_CONFIG.
LHAPDF_CONFIG = Path(os.environ.get(
    "HNL_LHAPDF_CONFIG",
    PROJECT_ROOT.parent.parent / "NNPDF40" / "fonll-local" / "env" / "bin" / "lhapdf-config",
))
LHAPDF_LIBDIR = LHAPDF_CONFIG.parent.parent / "lib"


def mg5_subprocess_env():
    """Process env for MG5 subprocesses: inject LHAPDF lib path on the loader.

    On macOS DYLD_LIBRARY_PATH is stripped by SIP across some shells but
    survives a direct subprocess.run(env=...) hand-off, which is how we
    call MG5. We set both DYLD_* and LD_* so the same wiring works on
    Linux without branching.
    """
    env = os.environ.copy()
    for key in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
        existing = env.get(key, "")
        env[key] = f"{LHAPDF_LIBDIR}:{existing}" if existing else str(LHAPDF_LIBDIR)
    return env


def patch_me5_configuration(cards_dir):
    """Per-job MG5 config: disable browser auto-open and wire LHAPDF.

    Writing `lhapdf_py3 = <path-to-lhapdf-config>` makes MadEvent discover
    the shared LHAPDF install (and therefore the NNPDF40 PDF grids) without
    touching the global mg5_configuration.txt under the shared MG5 install.
    """
    cfg = cards_dir / "me5_configuration.txt"
    if not cfg.exists():
        return
    text = cfg.read_text()
    new_text = text.replace("# automatic_html_opening = True", "automatic_html_opening = False")
    new_text = new_text.replace("# automatic_html_opening = False", "automatic_html_opening = False")
    new_text = new_text.replace("automatic_html_opening = True", "automatic_html_opening = False")
    new_text = new_text.replace("# web_browser = None", "web_browser = None")
    if "automatic_html_opening" not in new_text:
        new_text += "\nautomatic_html_opening = False\n"

    lhapdf_line = f"lhapdf_py3 = {LHAPDF_CONFIG}"
    if "# lhapdf_py3 = lhapdf-config" in new_text:
        new_text = new_text.replace("# lhapdf_py3 = lhapdf-config", lhapdf_line)
    elif "lhapdf_py3 =" not in new_text:
        new_text += f"\n{lhapdf_line}\n"

    if new_text != text:
        cfg.write_text(new_text)


def force_compile_subprocesses(work_subdir, log_path, timeout=900):
    """Force Fortran compilation of every SubProcess up-front.

    MG5's `output` step emits source only; the gensym/madevent binaries
    are compiled lazily on the first `generate_events`. We must drive a
    pre-compile so patch_rpath_for_lhapdf has real binaries to symlink
    next to. The interactive `madevent` shell accepts a `compile`
    command; we drive it via stdin and route output to log_path.

    Returns True if compilation succeeded (process exited 0 *and* at
    least one P*_*/madevent binary exists), False otherwise. Callers
    should refuse to continue on False because the dyld-rpath symlink
    step has nothing to patch and the subsequent generate_events will
    fail with a useless 'compile directory' error.
    """
    with open(log_path, "w") as log:
        result = subprocess.run(
            [str(PYTHON_EXE), "bin/madevent"],
            input="compile\nquit\n", text=True,
            cwd=work_subdir, env=mg5_subprocess_env(),
            stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
        )
    if result.returncode != 0:
        return False
    # Sanity check: at least one madevent binary must exist.
    for sub in (work_subdir / "SubProcesses").glob("P*_*"):
        if (sub / "madevent").exists():
            return True
    return False


def patch_rpath_for_lhapdf(work_subdir):
    """macOS-only: symlink libLHAPDF.dylib next to every compiled MG5 binary.

    Background: on macOS DYLD_LIBRARY_PATH is not reliably forwarded to
    MG5's grandchild Fortran subprocesses (MadEvent launches them through
    its own multiprocess machinery; SIP also rewrites the env for any
    binary loaded from a protected path). The compiled `madevent` /
    `gensym` binaries reference `@rpath/libLHAPDF.dylib` and their
    embedded rpath set is gcc's homebrew install dir.

    `install_name_tool -add_rpath` cannot extend the rpath in place
    because MG5's link step did not reserve header padding (the tool
    errors with "larger updated load commands do not fit"). Re-linking
    each subprocess would require modifying MG5's makefiles.

    The portable workaround: dyld searches `@loader_path` (the binary's
    own directory) when resolving `@rpath/...`. We drop a symlink to
    libLHAPDF.dylib into every SubProcesses/P*_*/ directory next to its
    madevent/gensym binary, and dyld picks it up.

    Idempotent; safe to re-run.
    """
    if sys.platform != "darwin":
        return
    src = LHAPDF_LIBDIR / "libLHAPDF.dylib"
    if not src.exists():
        return
    for sub in (work_subdir / "SubProcesses").glob("P*_*"):
        if not sub.is_dir():
            continue
        link = sub / "libLHAPDF.dylib"
        if link.exists() or link.is_symlink():
            continue
        link.symlink_to(src)
