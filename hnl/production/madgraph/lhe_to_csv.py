#!/usr/bin/env python3
"""
production/madgraph/lhe_to_csv.py

Convert MadGraph LHE files to CSV for HNL or tau 4-vectors.

Two modes:
  - HNL mode: Extract HNL (PDG 9900012) 4-vectors from W/Z → ℓ N events.
    Output: weight,E,px,py,pz (headerless)
  - Tau mode: Extract tau (PDG ±15) 4-vectors from W/Z → τ ν or τ+τ- events.
    Output: weight,E,px,py,pz (headerless)

Adapted from the upstream llpatcolliders_FONLL lhe_to_csv.py.
"""

import gzip
import numpy as np
from pathlib import Path


class LHEParser:
    """Parse LHE files and extract particle 4-vectors."""

    PDG_HNL = 9900012
    PDG_TAU_PLUS = -15
    PDG_TAU_MINUS = 15

    def __init__(self, lhe_path):
        self.lhe_path = Path(lhe_path)
        if not self.lhe_path.exists():
            raise FileNotFoundError(f"LHE file not found: {lhe_path}")

    def _open(self):
        if self.lhe_path.suffix == '.gz':
            return gzip.open(self.lhe_path, 'rt', encoding='utf-8')
        return open(self.lhe_path, 'r', encoding='utf-8')

    def extract_particles(self, target_pdg_ids, split_event_weight=False):
        """
        Parse LHE and extract 4-vectors for target PDG IDs.

        Parameters
        ----------
        target_pdg_ids : set of int
            PDG IDs to extract (e.g., {9900012} for HNL, {15, -15} for tau).
        split_event_weight : bool
            If True, split each event weight evenly across all extracted particles
            in that event.

        Yields
        ------
        dict with 'weight', 'E', 'px', 'py', 'pz', 'pdgid', 'origin' for each
        match. ``origin`` is the PDG of the target particle's mother (LHE
        MOTHUP1, 1-based index resolved against the same event), or 0 if the
        mother slot is empty / the index is out of range.
        """
        in_event = False
        header_parsed = False
        event_weight = 1.0
        all_event_pdgs = []     # PDG of every particle in the current event, in order
        candidate_records = []  # tuples (target_dict, mother_idx_1based)

        with self._open() as f:
            for line in f:
                stripped = line.strip()

                if stripped.startswith('<event>'):
                    in_event = True
                    header_parsed = False
                    event_weight = 1.0
                    all_event_pdgs = []
                    candidate_records = []
                    continue

                if stripped.startswith('</event>'):
                    if candidate_records:
                        if split_event_weight:
                            particle_weight = event_weight / len(candidate_records)
                        else:
                            particle_weight = event_weight
                        for record, mother_idx in candidate_records:
                            record['weight'] = particle_weight
                            # MOTHUP is 1-based; 0 means no mother. Resolve.
                            origin = 0
                            if 1 <= mother_idx <= len(all_event_pdgs):
                                origin = all_event_pdgs[mother_idx - 1]
                            record['origin'] = origin
                            yield record
                    in_event = False
                    continue

                if in_event:
                    if stripped.startswith('<') or stripped.startswith('#'):
                        continue

                    if not header_parsed:
                        parts = stripped.split()
                        if len(parts) >= 3:
                            event_weight = float(parts[2])
                        header_parsed = True
                        continue

                    parts = stripped.split()
                    if len(parts) >= 11:
                        try:
                            pdgid = int(parts[0])
                        except (ValueError, IndexError):
                            pdgid = 0  # malformed; keep ordering via placeholder
                        all_event_pdgs.append(pdgid)
                        if pdgid in target_pdg_ids:
                            # LHE columns (0-indexed): 0=IDUP 1=ISTUP
                            # 2=MOTHUP1 3=MOTHUP2 4=COL1 5=COL2
                            # 6=PUP1(px) 7=PUP2(py) 8=PUP3(pz) 9=PUP4(E)
                            # 10=PUP5(m) ...
                            try:
                                mother_idx = int(parts[2])
                                candidate_records.append((
                                    {
                                        'weight': event_weight,
                                        'pdgid': pdgid,
                                        'E': float(parts[9]),
                                        'px': float(parts[6]),
                                        'py': float(parts[7]),
                                        'pz': float(parts[8]),
                                    },
                                    mother_idx,
                                ))
                            except (ValueError, IndexError):
                                pass

    def write_hnl_csv(self, output_path):
        """Extract HNL 4-vectors and write headerless CSV."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        rows = []
        for p in self.extract_particles({self.PDG_HNL}):
            rows.append([p['weight'], p['E'], p['px'], p['py'], p['pz']])

        if rows:
            data = np.array(rows)
            np.savetxt(output_path, data, delimiter=",", fmt="%.8e")
        else:
            output_path.write_text("")

        return len(rows)

    def write_tau_csv(self, output_path):
        """Extract tau 4-vectors and write headerless CSV.

        Schema: ``weight, E, px, py, pz, origin`` (6 columns).
        ``origin`` is the mother PDG (24=W+, -24=W-, 23=Z, ...) so a
        downstream consumer can apply per-process polarisation when
        decaying tau -> N + X.

        Each tau row carries the *full* event weight. For HNL production every
        tau in the event is an independent potential parent, so a Z/gamma* ->
        tau+ tau- event contributes one expected N decay per tau (weight x BR
        per row), not half. Splitting the event weight across the two taus
        would underweight the Drell-Yan contribution by 2x. The W -> tau nu
        branch is unaffected (one tau per event, so split-vs-not is identical).
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        rows = []
        for p in self.extract_particles(
            {self.PDG_TAU_MINUS, self.PDG_TAU_PLUS},
            split_event_weight=False,
        ):
            rows.append([p['weight'], p['E'], p['px'], p['py'], p['pz'], p['origin']])

        if rows:
            data = np.array(rows)
            np.savetxt(output_path, data, delimiter=",", fmt="%.8e")
        else:
            output_path.write_text("")

        return len(rows)


def main():
    """Command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(description="Convert LHE to CSV")
    parser.add_argument("lhe_file", help="Input LHE file (.lhe or .lhe.gz)")
    parser.add_argument("csv_file", help="Output CSV file")
    parser.add_argument("--mode", choices=["hnl", "tau"], default="hnl",
                        help="Extraction mode (default: hnl)")
    args = parser.parse_args()

    lhe = LHEParser(args.lhe_file)
    if args.mode == "hnl":
        n = lhe.write_hnl_csv(args.csv_file)
    else:
        n = lhe.write_tau_csv(args.csv_file)

    print(f"Extracted {n} particles → {args.csv_file}")


if __name__ == "__main__":
    main()
