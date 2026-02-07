

SIGMA_CCBAR_PB = 24.0 * 1e9

SIGMA_BBBAR_PB = 500.0 * 1e6


FRAG_C_D0     = 0.59
FRAG_C_DPLUS  = 0.24
FRAG_C_DS     = 0.10
FRAG_C_LAMBDA = 0.06

FRAG_B_B0     = 0.40
FRAG_B_BPLUS  = 0.40
FRAG_B_BS     = 0.10
FRAG_B_LAMBDA = 0.10

SIGMA_KAON_PB = 5.0 * 1e10

SIGMA_KL_PB = SIGMA_KAON_PB * 0.5

SIGMA_W_PB = 2.0 * 1e8
SIGMA_Z_PB = 6.0 * 1e7


def get_parent_sigma_pb(parent_pdg: int) -> float:
    pid = abs(int(parent_pdg))

    if pid == 321:
        return SIGMA_KAON_PB
    if pid == 130:
        return SIGMA_KL_PB

    if pid == 421:
        return SIGMA_CCBAR_PB * FRAG_C_D0 * 2
    if pid == 411:
        return SIGMA_CCBAR_PB * FRAG_C_DPLUS * 2
    if pid == 431:
        return SIGMA_CCBAR_PB * FRAG_C_DS * 2
    if pid == 4122:
        return SIGMA_CCBAR_PB * FRAG_C_LAMBDA * 2

    if pid == 511:
        return SIGMA_BBBAR_PB * FRAG_B_B0 * 2
    if pid == 521:
        return SIGMA_BBBAR_PB * FRAG_B_BPLUS * 2
    if pid == 531:
        return SIGMA_BBBAR_PB * FRAG_B_BS * 2
    if pid == 541:
        return SIGMA_BBBAR_PB * 0.001 * 2
    if pid == 5122:
        return SIGMA_BBBAR_PB * FRAG_B_LAMBDA * 2

    if pid == 5232:
        return SIGMA_BBBAR_PB * 0.03 * 2
    if pid == 5332:
        return SIGMA_BBBAR_PB * 0.01 * 2

    if pid == 15:
        return (SIGMA_CCBAR_PB * FRAG_C_DS * 2) * 0.055

    if pid == 24:
        return SIGMA_W_PB
    if pid == 23:
        return SIGMA_Z_PB

    print(f"[WARNING] Unknown parent PDG {pid} in cross-section lookup. Returning 0.")
    return 0.0


def get_parent_tau_br(parent_pdg: int) -> float:
    pid = abs(int(parent_pdg))

    if pid == 431:
        return 0.053

    if pid == 511:
        return 0.023

    if pid == 521:
        return 0.023

    if pid == 531:
        return 0.023

    return 0.0


def get_sigma_summary():
    common_parents = [
        (321, "K±"),
        (130, "K_L"),
        (421, "D0"),
        (411, "D+"),
        (431, "Ds+"),
        (511, "B0"),
        (521, "B+"),
        (531, "Bs0"),
        (541, "Bc+"),
        (4122, "Λc+"),
        (5122, "Λb0"),
        (15, "τ"),
        (24, "W±"),
        (23, "Z"),
    ]

    print("=" * 70)
    print("LHC Production Cross-Sections (14 TeV, PBC Standard)")
    print("=" * 70)
    print(f"Base rates:")
    print(f"  σ(ccbar) = {SIGMA_CCBAR_PB:.2e} pb = {SIGMA_CCBAR_PB/1e9:.1f} mb")
    print(f"  σ(bbbar) = {SIGMA_BBBAR_PB:.2e} pb = {SIGMA_BBBAR_PB/1e6:.1f} μb")
    print()
    print("Parent meson cross-sections:")
    print(f"{'Parent':<10} {'PDG':>6} {'σ (pb)':>15} {'σ (nb)':>15}")
    print("-" * 70)

    for pdg, name in common_parents:
        sigma_pb = get_parent_sigma_pb(pdg)
        sigma_nb = sigma_pb / 1e3
        print(f"{name:<10} {pdg:>6} {sigma_pb:>15.3e} {sigma_nb:>15.3e}")

    print("=" * 70)


if __name__ == "__main__":
    get_sigma_summary()
