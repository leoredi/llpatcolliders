#!/usr/bin/env python3
"""
LHE to CSV converter for HNL events

Parses LesHouches Event (LHE) files from MadGraph and extracts:
- HNL particles (PDG 9900012 = N1)
- Parent W/Z bosons (PDG 24, -24, 23)
- 4-momenta and event weights

Output CSV format matches analysis pipeline requirements.

Usage:
    python lhe_to_csv.py input.lhe output.csv --mass 15.0 --flavour muon
    python lhe_to_csv.py input.lhe.gz output.csv --mass 15.0 --flavour muon
"""

import sys
import os
import gzip
import argparse
import re
from pathlib import Path


class LHEParser:
    """Parse LesHouches Event (LHE) files and extract HNL events"""

    # PDG codes
    PDG_HNL_N1 = 9900012
    PDG_WPLUS = 24
    PDG_WMINUS = -24
    PDG_Z = 23

    def __init__(self, lhe_path, mass_gev, flavour):
        """
        Initialize LHE parser

        Args:
            lhe_path: Path to LHE file (can be .lhe or .lhe.gz)
            mass_gev: HNL mass in GeV
            flavour: Lepton flavour (electron, muon, tau)
        """
        self.lhe_path = Path(lhe_path)
        self.mass_gev = float(mass_gev)
        self.flavour = flavour

        if not self.lhe_path.exists():
            raise FileNotFoundError(f"LHE file not found: {lhe_path}")

    def _open_lhe(self):
        """Open LHE file (handles both .lhe and .lhe.gz)"""
        if self.lhe_path.suffix == '.gz':
            return gzip.open(self.lhe_path, 'rt', encoding='utf-8')
        else:
            return open(self.lhe_path, 'r', encoding='utf-8')

    def parse_events(self):
        """
        Parse LHE file and yield event dictionaries

        Uses MATHUSLA approach: Extract HNL 4-vectors only.
        Parent W/Z may not appear in LHE if off-shell (controlled by bw_cut).

        Yields:
            dict with keys:
                - event_id: int
                - hnl_pdgid: int (N1)
                - mass_hnl_GeV: float
                - weight: float
                - hnl_E_GeV, hnl_px_GeV, hnl_py_GeV, hnl_pz_GeV: float
                - hnl_pt_GeV, hnl_eta, hnl_phi: float (derived)
        """
        event_id = 0
        in_event = False
        event_weight = 1.0
        particles = []
        header_parsed = False

        with self._open_lhe() as f:
            for line in f:
                line = line.strip()

                # Start of event block
                if line.startswith('<event>'):
                    in_event = True
                    particles = []
                    event_weight = 1.0
                    header_parsed = False
                    continue

                # End of event block
                if line.startswith('</event>'):
                    in_event = False
                    event_id += 1

                    # Extract HNL 4-vector (MATHUSLA approach)
                    event_data = self._extract_hnl(
                        particles, event_id, event_weight
                    )

                    if event_data is not None:
                        yield event_data

                    continue

                # Parse event content
                if in_event:
                    # Skip XML tags and comments
                    if line.startswith('<') or line.startswith('#'):
                        continue

                    # First non-comment line is event header
                    if not header_parsed:
                        # Event header format: nup idprup xwgtup scalup aqedup aqcdup
                        parts = line.split()
                        if len(parts) >= 3:
                            event_weight = float(parts[2])  # xwgtup
                        header_parsed = True
                        continue

                    # Subsequent lines are particles
                    if line:
                        parts = line.split()
                        if len(parts) >= 11:
                            try:
                                particle = {
                                    'pdgid': int(parts[0]),
                                    'status': int(parts[1]),
                                    'mother1': int(parts[2]),
                                    'mother2': int(parts[3]),
                                    'px': float(parts[6]),
                                    'py': float(parts[7]),
                                    'pz': float(parts[8]),
                                    'E': float(parts[9]),
                                    'mass': float(parts[10]),
                                }
                                particles.append(particle)
                            except (ValueError, IndexError):
                                # Skip malformed lines
                                pass

        print(f"Parsed {event_id} events from {self.lhe_path.name}")

    def _find_parent_boson(self, particles, hnl):
        """
        Traverse mother chain to find parent W/Z boson.

        Args:
            particles: List of particle dicts from LHE event
            hnl: HNL particle dict

        Returns:
            int: PDG ID of parent W/Z, or 0 if not found
        """
        visited = set()
        current_idx = hnl['mother1']

        while 1 <= current_idx <= len(particles) and current_idx not in visited:
            visited.add(current_idx)
            parent = particles[current_idx - 1]  # LHE indices are 1-based
            if parent['pdgid'] in [self.PDG_WPLUS, self.PDG_WMINUS, self.PDG_Z]:
                return parent['pdgid']
            current_idx = parent['mother1']

        return 0  # Parent boson not found in chain

    def _extract_hnl(self, particles, event_id, weight):
        """
        Extract HNL 4-vector and parent info from particle list

        Attempts to extract parent W/Z if present in LHE. If not found
        (can happen for off-shell bosons), uses default parent_pdg=0.

        Output format matches Pythia CSV for analysis pipeline compatibility.

        Args:
            particles: List of particle dicts
            event_id: Event number
            weight: Event weight

        Returns:
            dict with event data, or None if no HNL found
        """
        import math

        # Find HNL (N1) - should be only one per event
        hnl = None
        for p in particles:
            if p['pdgid'] == self.PDG_HNL_N1:
                hnl = p
                break

        if hnl is None:
            return None

        # Try to find parent W/Z (may not exist if off-shell)
        parent_pdg = 0  # Default if not found
        mother1_idx = hnl['mother1']

        if 1 <= mother1_idx <= len(particles):
            parent_candidate = particles[mother1_idx - 1]
            if parent_candidate['pdgid'] in [self.PDG_WPLUS, self.PDG_WMINUS, self.PDG_Z]:
                parent_pdg = parent_candidate['pdgid']

        # Fallback: traverse mother chain to find parent W/Z
        if parent_pdg == 0:
            parent_pdg = self._find_parent_boson(particles, hnl)

        # Extract 4-momentum
        px = hnl['px']
        py = hnl['py']
        pz = hnl['pz']
        E = hnl['E']

        # Compute derived quantities
        pt = math.sqrt(px**2 + py**2)
        p = math.sqrt(px**2 + py**2 + pz**2)

        # Eta (pseudorapidity): η = -ln(tan(θ/2))
        # Use theta-based calculation for numerical stability (avoids log of negative numbers)
        if pt < 1e-10:
            eta = 999.0 if pz > 0 else -999.0
        else:
            theta = math.atan2(pt, pz)
            if theta < 1e-10:
                eta = 999.0
            elif theta > math.pi - 1e-10:
                eta = -999.0
            else:
                eta = -math.log(math.tan(theta / 2.0))

        # Phi (azimuthal angle)
        phi = math.atan2(py, px)

        # Boost factor: β γ = p / m (NOT the Lorentz factor γ = E / m)
        beta_gamma = p / self.mass_gev if self.mass_gev > 0 else 0.0

        # Production vertex: MadGraph produces at IP, so (0,0,0) in mm
        prod_x_mm = 0.0
        prod_y_mm = 0.0
        prod_z_mm = 0.0

        # Build event data dictionary (EXACT Pythia CSV format)
        return {
            'event': event_id,
            'weight': weight,
            'hnl_id': self.PDG_HNL_N1,
            'parent_pdg': parent_pdg,
            'pt': pt,
            'eta': eta,
            'phi': phi,
            'p': p,
            'E': E,
            'mass': self.mass_gev,
            'prod_x_mm': prod_x_mm,
            'prod_y_mm': prod_y_mm,
            'prod_z_mm': prod_z_mm,
            'beta_gamma': beta_gamma,
        }

    def write_csv(self, output_path):
        """
        Parse LHE and write CSV file

        Output format EXACTLY matches Pythia CSV from meson production:
        event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma

        Args:
            output_path: Path to output CSV file

        Returns:
            int: Number of events written
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # CSV header (EXACT Pythia format for analysis pipeline compatibility)
        header = "event,weight,hnl_id,parent_pdg,pt,eta,phi,p,E,mass,prod_x_mm,prod_y_mm,prod_z_mm,beta_gamma"

        n_events = 0
        n_no_parent = 0  # Count events where parent W/Z not found

        with open(output_path, 'w') as f:
            f.write(header + '\n')

            for event in self.parse_events():
                # Track missing parents
                if event['parent_pdg'] == 0:
                    n_no_parent += 1

                # Write CSV row (EXACT format, values NOT in scientific notation for compatibility)
                row = (
                    f"{event['event']},"
                    f"{event['weight']},"
                    f"{event['hnl_id']},"
                    f"{event['parent_pdg']},"
                    f"{event['pt']:.6g},"
                    f"{event['eta']:.6g},"
                    f"{event['phi']:.6g},"
                    f"{event['p']:.6g},"
                    f"{event['E']:.6g},"
                    f"{event['mass']:.6g},"
                    f"{event['prod_x_mm']:.6g},"
                    f"{event['prod_y_mm']:.6g},"
                    f"{event['prod_z_mm']:.6g},"
                    f"{event['beta_gamma']:.6g}"
                )
                f.write(row + '\n')
                n_events += 1

        print(f"Wrote {n_events} HNL events to {output_path}")
        if n_no_parent > 0:
            print(f"  Warning: {n_no_parent}/{n_events} events had parent_pdg=0 (W/Z not in LHE)")
            print(f"           This is expected for off-shell bosons (controlled by bw_cut)")

        return n_events


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description='Convert LHE file to CSV for HNL analysis'
    )
    parser.add_argument(
        'lhe_file',
        help='Input LHE file (.lhe or .lhe.gz)'
    )
    parser.add_argument(
        'csv_file',
        help='Output CSV file'
    )
    parser.add_argument(
        '--mass',
        type=float,
        required=True,
        help='HNL mass in GeV'
    )
    parser.add_argument(
        '--flavour',
        choices=['electron', 'muon', 'tau'],
        required=True,
        help='Lepton flavour'
    )

    args = parser.parse_args()

    # Parse and convert
    lhe_parser = LHEParser(args.lhe_file, args.mass, args.flavour)
    n_events = lhe_parser.write_csv(args.csv_file)

    print(f"\nConversion complete:")
    print(f"  Input:  {args.lhe_file}")
    print(f"  Output: {args.csv_file}")
    print(f"  Events: {n_events}")
    print(f"  Mass:   {args.mass} GeV")
    print(f"  Flavour: {args.flavour}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
