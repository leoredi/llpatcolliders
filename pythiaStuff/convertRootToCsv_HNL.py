#!/usr/bin/env python3
"""
Convert ROOT file from PYTHIA simulation to CSV format for HNL analysis.
Extracts HNL particles (PDG ID 9900014 and -9900014) and their kinematic properties.
"""

import ROOT as r
import csv
import sys
import numpy as np

# Don't show ROOT graphics
r.gROOT.SetBatch()

def calculate_pt_eta_phi(px, py, pz):
    """
    Calculate transverse momentum (pt), pseudorapidity (eta), and azimuthal angle (phi)
    from momentum components.

    Args:
        px, py, pz: Momentum components in GeV/c

    Returns:
        pt, eta, phi
    """
    pt = np.sqrt(px**2 + py**2)
    p = np.sqrt(px**2 + py**2 + pz**2)

    # Calculate eta
    if p > 0:
        eta = 0.5 * np.log((p + pz) / (p - pz)) if abs(p - pz) > 1e-10 else 0.0
    else:
        eta = 0.0

    # Calculate phi
    phi = np.arctan2(py, px)

    return pt, eta, phi

# Input/output files
input_filename = "main144.root"
output_filename = "../LLP.csv"

if len(sys.argv) > 1:
    input_filename = sys.argv[1]
if len(sys.argv) > 2:
    output_filename = sys.argv[2]

print(f"Reading ROOT file: {input_filename}")
print(f"Output CSV file: {output_filename}")

# Open ROOT file
try:
    inputFile = r.TFile(input_filename)
    if not inputFile or inputFile.IsZombie():
        print(f"Error: Could not open file {input_filename}")
        sys.exit(1)
except Exception as e:
    print(f"Error opening ROOT file: {e}")
    sys.exit(1)

# Get tree
inputTree = inputFile.Get("t")
if not inputTree:
    print("Error: Could not find tree 't' in ROOT file")
    sys.exit(1)

print(f"Found tree with {inputTree.GetEntries()} entries")

# PDG ID for HNL (muon-type HNL and anti-HNL)
HNL_PDG = 9900014
ANTI_HNL_PDG = -9900014

# Open CSV file for writing
with open(output_filename, 'w', newline='') as csvfile:
    fieldnames = ['event', 'id', 'pt', 'eta', 'phi', 'momentum', 'mass']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    event_count = 0
    hnl_count = 0

    # Process each event
    for event_idx, event in enumerate(inputTree):
        # Find HNL particles in this event
        hnls_in_event = []

        # Access particles in the event
        # The structure depends on how main144 stores particles
        # Typical structure: event.particles is a collection

        try:
            if hasattr(event, 'particles'):
                particles = event.particles
            elif hasattr(event, 'particle'):
                particles = event.particle
            else:
                # If particles are stored differently, try direct access
                # This may need adjustment based on actual ROOT file structure
                print("Warning: Could not determine particle storage structure")
                continue

            for particle in particles:
                # Check if this is an HNL
                pid = particle.pid if hasattr(particle, 'pid') else particle.id

                if abs(pid) == HNL_PDG:
                    # Get momentum components
                    px = particle.px if hasattr(particle, 'px') else 0.0
                    py = particle.py if hasattr(particle, 'py') else 0.0
                    pz = particle.pz if hasattr(particle, 'pz') else 0.0

                    # Calculate pt, eta, phi
                    pt, eta, phi = calculate_pt_eta_phi(px, py, pz)

                    # Calculate total momentum
                    momentum = np.sqrt(px**2 + py**2 + pz**2)

                    # Get mass
                    if hasattr(particle, 'mass'):
                        mass = particle.mass
                    elif hasattr(particle, 'm'):
                        mass = particle.m
                    else:
                        # Calculate mass from energy and momentum if needed
                        if hasattr(particle, 'energy') or hasattr(particle, 'e'):
                            energy = particle.energy if hasattr(particle, 'energy') else particle.e
                            mass = np.sqrt(max(0, energy**2 - momentum**2))
                        else:
                            mass = 5.0  # Default HNL mass

                    hnls_in_event.append({
                        'event': event_count,
                        'id': pid,
                        'pt': pt,
                        'eta': eta,
                        'phi': phi,
                        'momentum': momentum,
                        'mass': mass
                    })
                    hnl_count += 1

            # Write HNLs from this event to CSV
            if hnls_in_event:
                for hnl in hnls_in_event:
                    writer.writerow(hnl)
                event_count += 1

        except Exception as e:
            print(f"Warning: Error processing event {event_idx}: {e}")
            continue

        # Progress update
        if (event_idx + 1) % 1000 == 0:
            print(f"Processed {event_idx + 1} events, found {hnl_count} HNLs in {event_count} events with HNLs")

inputFile.Close()

print(f"\nConversion complete!")
print(f"Total events processed: {inputTree.GetEntries()}")
print(f"Events with HNLs: {event_count}")
print(f"Total HNLs found: {hnl_count}")
print(f"Output written to: {output_filename}")
