"""
Per-Parent Geometric Efficiency Calculator

Computes ε_geom(parent, mass, cτ) for the drainage gallery detector.

This module:
1. Reads Pythia CSV with parent_id tracking
2. Groups events by parent species (K, D, B, W, Z, ...)
3. Computes decay probability in detector volume for each parent
4. Outputs per-parent efficiency maps: eps_geom[parent_pdg][ctau]

Replaces the old single-efficiency approach with proper parent separation
for PBC-grade analysis.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pickle
from dataclasses import dataclass


@dataclass
class DetectorGeometry:
    """
    Drainage gallery detector geometry.

    Coordinates relative to CMS interaction point.
    """
    # Detector position (m)
    x_min: float = -5.0
    x_max: float = 5.0
    y_min: float = 18.0  # ~20m above IP
    y_max: float = 22.0
    z_min: float = -5.0
    z_max: float = 5.0

    # Shielding (for background studies, not used in acceptance)
    concrete_thickness_m: float = 1.0
    earth_thickness_m: float = 18.0

    def volume_m3(self) -> float:
        """Active detector volume in m³."""
        return ((self.x_max - self.x_min) *
                (self.y_max - self.y_min) *
                (self.z_max - self.z_min))

    def contains_point(self, x: float, y: float, z: float) -> bool:
        """Check if a point is inside the detector volume."""
        return (self.x_min <= x <= self.x_max and
                self.y_min <= y <= self.y_max and
                self.z_min <= z <= self.z_min)


class PerParentEfficiency:
    """
    Calculate geometric efficiency per parent species.

    For each (mass, flavour) mass point:
    - Reads CSV from Pythia with parent_id column
    - Groups events by parent PDG ID
    - For each parent and cτ hypothesis:
      - Computes P(decay in drainage gallery | kinematics, cτ)
    - Saves efficiency map: eps_geom[parent_pdg][ctau]
    """

    def __init__(self, detector: Optional[DetectorGeometry] = None):
        """
        Initialize efficiency calculator.

        Parameters:
        -----------
        detector : DetectorGeometry, optional
            Detector geometry (default: drainage gallery)
        """
        self.detector = detector if detector else DetectorGeometry()
        self.ctau_scan_points = np.logspace(-2, 3, 100)  # 0.01 m to 1 km

    def load_csv(self, csv_path: Path) -> pd.DataFrame:
        """
        Load Pythia CSV with production vertices and parent IDs.

        Expected columns:
        - event, weight, id, parent_id
        - pt, eta, phi, momentum, energy, mass
        - prod_x_m, prod_y_m, prod_z_m
        """
        df = pd.read_csv(csv_path)

        # Verify required columns
        required = ['event', 'weight', 'id', 'parent_id',
                   'pt', 'eta', 'phi', 'momentum', 'energy', 'mass',
                   'prod_x_m', 'prod_y_m', 'prod_z_m']

        missing = set(required) - set(df.columns)
        if missing:
            raise ValueError(f"CSV missing columns: {missing}")

        return df

    def decay_probability_in_detector(self, df: pd.DataFrame,
                                      ctau_m: float) -> np.ndarray:
        """
        Calculate P(decay in detector) for each event at given cτ.

        Uses proper Lorentz boost: decay length in lab frame is γβcτ.

        Parameters:
        -----------
        df : pd.DataFrame
            Event data with kinematics and production vertices
        ctau_m : float
            Proper lifetime cτ in meters

        Returns:
        --------
        np.ndarray
            Decay probability for each event (0 to 1)
        """
        # Production vertex
        x0 = df['prod_x_m'].values
        y0 = df['prod_y_m'].values
        z0 = df['prod_z_m'].values

        # HNL momentum and energy
        p = df['momentum'].values
        E = df['energy'].values
        m = df['mass'].values

        # Boost factor: γβ = p/m
        gamma_beta = p / m

        # Mean decay length in lab frame
        L_decay_lab = gamma_beta * ctau_m

        # Direction unit vector
        pt = df['pt'].values
        eta = df['eta'].values
        phi = df['phi'].values

        # Momentum components (unit vector)
        px_hat = pt * np.cos(phi) / p
        py_hat = pt * np.sin(phi) / p
        pz_hat = np.sqrt(1 - (pt/p)**2) * np.sign(np.sinh(eta))

        # Ray-tracing: find intersection with detector box
        # For each event, compute distance to detector entry and exit

        # Simplified: compute distance to detector center plane (y ≈ 20 m)
        # More sophisticated: full box intersection

        y_detector_center = (self.detector.y_max + self.detector.y_min) / 2

        # Distance to y = y_detector_center plane
        # y0 + t * py_hat = y_detector_center
        # t = (y_detector_center - y0) / py_hat

        with np.errstate(divide='ignore', invalid='ignore'):
            t_to_detector = (y_detector_center - y0) / py_hat
            t_to_detector = np.where(py_hat > 0, t_to_detector, np.inf)

        # Decay probability between production and detector
        # P(not decay before detector) = exp(-t_to_detector / L_decay)
        P_survive_to_detector = np.exp(-t_to_detector / L_decay_lab)

        # P(decay in detector volume) ≈ P(survive to detector) × (detector_thickness / L_decay)
        # Simplified: assume detector thickness is small compared to L_decay
        detector_thickness = self.detector.y_max - self.detector.y_min

        P_decay_in_detector = P_survive_to_detector * (detector_thickness / L_decay_lab)

        # Clamp to [0, 1]
        P_decay_in_detector = np.clip(P_decay_in_detector, 0, 1)

        return P_decay_in_detector

    def compute_efficiency_per_parent(self, csv_path: Path,
                                      output_dir: Optional[Path] = None) -> Dict:
        """
        Compute ε_geom(parent, cτ) for all parents in the sample.

        Parameters:
        -----------
        csv_path : Path
            Path to Pythia CSV file
        output_dir : Path, optional
            Directory to save efficiency maps (default: same as CSV)

        Returns:
        --------
        Dict
            {parent_pdg: {ctau: eps_geom}}
        """
        # Load data
        df = self.load_csv(csv_path)

        # Group by parent_id
        parent_groups = df.groupby('parent_id')

        efficiency_map = {}

        print(f"Computing efficiency for {len(parent_groups)} parent species...")

        for parent_pdg, group_df in parent_groups:
            print(f"  Parent PDG {parent_pdg}: {len(group_df)} events")

            eps_vs_ctau = {}

            for ctau in self.ctau_scan_points:
                # Compute decay probability for this parent at this cτ
                P_decay = self.decay_probability_in_detector(group_df, ctau)

                # Weighted average over events
                weights = group_df['weight'].values
                eps_geom = np.average(P_decay, weights=weights)

                eps_vs_ctau[ctau] = eps_geom

            efficiency_map[parent_pdg] = eps_vs_ctau

        # Save to disk
        if output_dir is None:
            output_dir = csv_path.parent

        # Extract mass and flavour from filename
        # Expected: HNL_mass_1.0_muon_Meson.csv
        stem = csv_path.stem
        output_file = output_dir / f"{stem}_efficiency_map.pkl"

        with open(output_file, 'wb') as f:
            pickle.dump(efficiency_map, f)

        print(f"Saved efficiency map to {output_file}")

        # Also save as CSV for inspection
        self._save_efficiency_csv(efficiency_map, output_dir / f"{stem}_efficiency_map.csv")

        return efficiency_map

    def _save_efficiency_csv(self, efficiency_map: Dict, output_path: Path):
        """Save efficiency map as CSV for easy inspection."""
        rows = []
        for parent_pdg, eps_vs_ctau in efficiency_map.items():
            for ctau, eps in eps_vs_ctau.items():
                rows.append({
                    'parent_pdg': parent_pdg,
                    'ctau_m': ctau,
                    'eps_geom': eps
                })

        df = pd.DataFrame(rows)
        df.to_csv(output_path, index=False)
        print(f"Saved efficiency CSV to {output_path}")

    def load_efficiency_map(self, pkl_path: Path) -> Dict:
        """Load pre-computed efficiency map from pickle."""
        with open(pkl_path, 'rb') as f:
            return pickle.load(f)


def process_all_csvs(csv_dir: Path, output_dir: Optional[Path] = None):
    """
    Batch process all CSV files in a directory.

    Parameters:
    -----------
    csv_dir : Path
        Directory containing Pythia CSV files
    output_dir : Path, optional
        Output directory for efficiency maps
    """
    calculator = PerParentEfficiency()

    csv_files = sorted(csv_dir.glob("HNL_mass_*_*.csv"))

    print(f"Found {len(csv_files)} CSV files to process")

    for i, csv_path in enumerate(csv_files, 1):
        print(f"\n[{i}/{len(csv_files)}] Processing {csv_path.name}...")
        try:
            calculator.compute_efficiency_per_parent(csv_path, output_dir)
        except Exception as e:
            print(f"ERROR processing {csv_path.name}: {e}")
            continue

    print("\nAll files processed.")


# Example usage
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python per_parent_efficiency.py <csv_file_or_directory>")
        sys.exit(1)

    path = Path(sys.argv[1])

    if path.is_file():
        # Single file
        calc = PerParentEfficiency()
        calc.compute_efficiency_per_parent(path)
    elif path.is_dir():
        # Batch process directory
        process_all_csvs(path)
    else:
        print(f"ERROR: {path} is not a valid file or directory")
        sys.exit(1)
