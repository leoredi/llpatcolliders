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

        Yields:
            dict with keys:
                - event_id: int
                - parent_pdgid: int (W+/W-/Z)
                - hnl_pdgid: int (N1)
                - mass_hnl_GeV: float
                - weight: float
                - parent_E_GeV, parent_px_GeV, parent_py_GeV, parent_pz_GeV: float
                - hnl_E_GeV, hnl_px_GeV, hnl_py_GeV, hnl_pz_GeV: float
        """
        event_id = 0
        in_event = False
        event_weight = 1.0
        particles = []

        with self._open_lhe() as f:
            for line in f:
                line = line.strip()

                # Start of event block
                if line.startswith('<event>'):
                    in_event = True
                    particles = []
                    event_weight = 1.0
                    continue

                # End of event block
                if line.startswith('</event>'):
                    in_event = False
                    event_id += 1

                    # Extract HNL and parent information
                    event_data = self._extract_hnl_and_parent(
                        particles, event_id, event_weight
                    )

                    if event_data is not None:
                        yield event_data

                    continue

                # Parse event content
                if in_event:
                    # First line after <event> is event header
                    if len(particles) == 0 and not line.startswith('#'):
                        # Event header format: nup idprup xwgtup scalup aqedup aqcdup
                        parts = line.split()
                        if len(parts) >= 3:
                            event_weight = float(parts[2])  # xwgtup

                    # Particle lines
                    elif not line.startswith('#') and line:
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

    def _extract_hnl_and_parent(self, particles, event_id, weight):
        """
        Extract HNL and its parent W/Z from particle list

        Args:
            particles: List of particle dicts
            event_id: Event number
            weight: Event weight

        Returns:
            dict with event data, or None if no HNL found
        """
        # Find HNL (N1)
        hnl = None
        for i, p in enumerate(particles):
            if p['pdgid'] == self.PDG_HNL_N1:
                hnl = p
                hnl_index = i + 1  # LHE uses 1-based indexing
                break

        if hnl is None:
            return None

        # Find parent W/Z using mother indices
        parent = None
        mother1_idx = hnl['mother1']

        if 1 <= mother1_idx <= len(particles):
            parent_candidate = particles[mother1_idx - 1]
            if parent_candidate['pdgid'] in [self.PDG_WPLUS, self.PDG_WMINUS, self.PDG_Z]:
                parent = parent_candidate

        # If direct mother is not W/Z, search all particles
        # (sometimes there are intermediate resonances)
        if parent is None:
            for p in particles:
                if p['pdgid'] in [self.PDG_WPLUS, self.PDG_WMINUS, self.PDG_Z]:
                    # Check if this W/Z is outgoing (status 1) or intermediate (status 2)
                    parent = p
                    break

        if parent is None:
            # No W/Z found - this shouldn't happen for our processes
            # but handle gracefully
            print(f"Warning: Event {event_id} has HNL but no W/Z parent found")
            return None

        # Build event data dictionary
        return {
            'event_id': event_id,
            'parent_pdgid': parent['pdgid'],
            'hnl_pdgid': hnl['pdgid'],
            'mass_hnl_GeV': self.mass_gev,
            'weight': weight,
            'parent_E_GeV': parent['E'],
            'parent_px_GeV': parent['px'],
            'parent_py_GeV': parent['py'],
            'parent_pz_GeV': parent['pz'],
            'hnl_E_GeV': hnl['E'],
            'hnl_px_GeV': hnl['px'],
            'hnl_py_GeV': hnl['py'],
            'hnl_pz_GeV': hnl['pz'],
        }

    def write_csv(self, output_path):
        """
        Parse LHE and write CSV file

        Args:
            output_path: Path to output CSV file

        Returns:
            int: Number of events written
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # CSV header (EXACT format required by analysis pipeline)
        header = "event_id,parent_pdgid,hnl_pdgid,mass_hnl_GeV,weight,parent_E_GeV,parent_px_GeV,parent_py_GeV,parent_pz_GeV,hnl_E_GeV,hnl_px_GeV,hnl_py_GeV,hnl_pz_GeV"

        n_events = 0
        with open(output_path, 'w') as f:
            f.write(header + '\n')

            for event in self.parse_events():
                # Write CSV row
                row = (
                    f"{event['event_id']},"
                    f"{event['parent_pdgid']},"
                    f"{event['hnl_pdgid']},"
                    f"{event['mass_hnl_GeV']:.1f},"
                    f"{event['weight']:.6e},"
                    f"{event['parent_E_GeV']:.6e},"
                    f"{event['parent_px_GeV']:.6e},"
                    f"{event['parent_py_GeV']:.6e},"
                    f"{event['parent_pz_GeV']:.6e},"
                    f"{event['hnl_E_GeV']:.6e},"
                    f"{event['hnl_px_GeV']:.6e},"
                    f"{event['hnl_py_GeV']:.6e},"
                    f"{event['hnl_pz_GeV']:.6e}"
                )
                f.write(row + '\n')
                n_events += 1

        print(f"Wrote {n_events} HNL events to {output_path}")
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
