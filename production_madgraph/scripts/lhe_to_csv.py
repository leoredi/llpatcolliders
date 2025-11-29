import os
import sys
import gzip
import csv
import argparse
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_lhe_event(lines, event_id, mass_hnl):
    """
    Parses a single event block from LHE.
    Returns a dictionary of event data or None if HNL/Parent not found.
    """
    # First line: Npart ID weight scale alpha_em alpha_s
    header = lines[0].strip().split()
    try:
        weight = float(header[2])
    except (IndexError, ValueError):
        logger.warning(f"Event {event_id}: Could not parse weight from '{lines[0]}'")
        return None

    particles = []
    # Parse particles
    # LHE particle line: pdg status mother1 mother2 color1 color2 px py pz E mass spin (lifetime)
    for i, line in enumerate(lines[1:]):
        parts = line.strip().split()
        if len(parts) < 10:
            continue

        try:
            pdg = int(parts[0])
            status = int(parts[1])
            mother1 = int(parts[2])
            mother2 = int(parts[3])
            px = float(parts[6])
            py = float(parts[7])
            pz = float(parts[8])
            e = float(parts[9])

            # 1-based index to 0-based
            particles.append({
                'index': i + 1,
                'pdg': pdg,
                'status': status,
                'mother1': mother1,
                'px': px,
                'py': py,
                'pz': pz,
                'e': e
            })
        except ValueError:
            continue

    # Find HNL (9900012)
    hnl = None
    for p in particles:
        if abs(p['pdg']) == 9900012:
            hnl = p
            break

    if not hnl:
        # It's possible we missed it or something went wrong
        # logger.debug(f"Event {event_id}: No HNL found.")
        return None

    # Find Parent (W: 24, Z: 23)
    parent = None
    m_idx = hnl['mother1']

    # Strategy 1: Explicit parent in event record
    if m_idx > 0 and m_idx <= len(particles):
        potential_parent = particles[m_idx - 1]
        if abs(potential_parent['pdg']) in [23, 24]:
            parent = potential_parent

    # Strategy 2: Reconstruct parent from final state siblings
    # If MG didn't write the W/Z, the HNL's mother is likely the beam particles (idx 1 or 2)
    # The sibling lepton + HNL = Parent W/Z.
    if not parent:
        # Find the sibling lepton.
        # It should have the same mother as HNL (or come from same vertex)
        # and be a lepton (11, 12, 13, 14, 15, 16)
        siblings = []
        for p in particles:
            if p == hnl:
                continue
            # Check if same mother (simple check)
            if p['mother1'] == hnl['mother1']:
                 # Check if lepton
                 if abs(p['pdg']) in [11, 12, 13, 14, 15, 16]:
                     siblings.append(p)

        if len(siblings) == 1:
            sib = siblings[0]
            # Reconstruct parent
            parent_px = hnl['px'] + sib['px']
            parent_py = hnl['py'] + sib['py']
            parent_pz = hnl['pz'] + sib['pz']
            parent_e  = hnl['e'] + sib['e']

            # Determine PDG
            # Charge conservation: Q_parent = Q_hnl + Q_sib
            # HNL is neutral (N1).
            # If sib is e- (11), Q=-1 -> Parent W- (-24)
            # If sib is e+ (-11), Q=+1 -> Parent W+ (24)
            # If sib is nu (12,14,16), Q=0 -> Parent Z (23)
            # If sib is anti-nu (-12,-14,-16), Q=0 -> Parent Z (23)

            sib_pdg = sib['pdg']
            parent_pdg = 0
            if abs(sib_pdg) in [11, 13, 15]: # Charged lepton
                if sib_pdg > 0: # e- (11) has pdg>0 in Pythia? No.
                # PDG: e- is 11. Charge -1.
                # e+ is -11. Charge +1.
                # W+ is 24. W- is -24.
                    parent_pdg = -24 # W- -> l- N
                else:
                    parent_pdg = 24  # W+ -> l+ N
            elif abs(sib_pdg) in [12, 14, 16]: # Neutrino
                parent_pdg = 23 # Z -> nu N

            if parent_pdg != 0:
                parent = {
                    'pdg': parent_pdg,
                    'px': parent_px,
                    'py': parent_py,
                    'pz': parent_pz,
                    'e': parent_e
                }

    if not parent:
        # Could not identify or reconstruct parent
        return None

    return {
        'event_id': event_id,
        'parent_pdgid': parent['pdg'],
        'hnl_pdgid': hnl['pdg'],
        'mass_hnl_GeV': mass_hnl,
        'weight': weight,
        'parent_E_GeV': parent['e'],
        'parent_px_GeV': parent['px'],
        'parent_py_GeV': parent['py'],
        'parent_pz_GeV': parent['pz'],
        'hnl_E_GeV': hnl['e'],
        'hnl_px_GeV': hnl['px'],
        'hnl_py_GeV': hnl['py'],
        'hnl_pz_GeV': hnl['pz']
    }

def process_lhe_file(lhe_path, output_csv_path, mass_hnl):
    logger.info(f"Processing {lhe_path} -> {output_csv_path}")

    if lhe_path.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    events_data = []

    try:
        with open_func(lhe_path, mode) as f:
            in_event = False
            event_lines = []
            event_count = 0

            for line in f:
                if '<event>' in line:
                    in_event = True
                    event_lines = []
                    continue
                if '</event>' in line:
                    in_event = False
                    event_count += 1
                    data = parse_lhe_event(event_lines, event_count, mass_hnl)
                    if data:
                        events_data.append(data)
                    continue

                if in_event:
                    event_lines.append(line)

    except FileNotFoundError:
        logger.error(f"File not found: {lhe_path}")
        return False
    except Exception as e:
        logger.error(f"Error processing LHE: {e}")
        return False

    # Write to CSV
    # Headers: event_id,parent_pdgid,hnl_pdgid,mass_hnl_GeV,weight,parent_E_GeV,parent_px_GeV,parent_py_GeV,parent_pz_GeV,hnl_E_GeV,hnl_px_GeV,hnl_py_GeV,hnl_pz_GeV
    headers = [
        'event_id', 'parent_pdgid', 'hnl_pdgid', 'mass_hnl_GeV', 'weight',
        'parent_E_GeV', 'parent_px_GeV', 'parent_py_GeV', 'parent_pz_GeV',
        'hnl_E_GeV', 'hnl_px_GeV', 'hnl_py_GeV', 'hnl_pz_GeV'
    ]

    if os.path.dirname(output_csv_path):
        os.makedirs(os.path.dirname(output_csv_path), exist_ok=True)

    with open(output_csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(events_data)

    logger.info(f"Wrote {len(events_data)} events to {output_csv_path}")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert LHE to CSV for HNL analysis.')
    parser.add_argument('lhe_file', help='Path to input LHE file')
    parser.add_argument('csv_file', help='Path to output CSV file')
    parser.add_argument('mass', type=float, help='HNL Mass in GeV')

    args = parser.parse_args()

    success = process_lhe_file(args.lhe_file, args.csv_file, args.mass)
    if not success:
        sys.exit(1)
