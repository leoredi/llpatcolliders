# Production System Comparison: Old vs New

## TL;DR: Two Different Production Systems

You now have **TWO completely independent production systems** in the same directory:

### OLD System (Original, before Nov 26)
- ❌ **No longer used** (replaced by new system)
- Files: `main_hnl_single.cc`, `make.sh`, old `.cmnd` templates
- Used shell script (`make.sh`) as build + run system
- Had physics bugs (neutral meson decays, mixed SoftQCD+HardQCD)

### NEW System (Gold Standard, after Nov 26)
- ✅ **Current active system**
- Files: `main_hnl_production.cc`, `Makefile`, `cards/*.cmnd`, `run_benchmarks.sh`
- Uses proper Makefile for building
- Fixed all physics bugs
- Publication-quality

---

## Why Two Build Systems?

### Old System: `make.sh` (Shell Script Approach)

**File:** `make.sh` (from git history, no longer in repo)

**What it did:**
```bash
./make.sh all         # Compile + run all mass points
./make.sh electron    # Compile + run electron coupling
./make.sh muon        # Compile + run muon coupling
```

**How it worked:**
1. Auto-compile Pythia if needed
2. Compile `main_hnl_single.cc` with g++ directly
3. Loop over mass points
4. Run simulation for each mass
5. Output to `../output/csv/simulation/`

**Problems:**
- Shell script doing double duty (build + orchestration)
- Hard to maintain
- Not standard Unix workflow
- Physics bugs in the code

---

### New System: `Makefile` + `run_benchmarks.sh` (Professional Approach)

**Build system:** `Makefile` (standard Unix build tool)

**What it does:**
```bash
make                  # Just compile the binary
make clean            # Remove compiled files
make test             # Quick test run
```

**Orchestration:** `run_benchmarks.sh` (separate from building)

**What it does:**
```bash
./run_benchmarks.sh              # Run all benchmarks
./run_benchmarks.sh 50000        # Run with 50k events/point
```

**How it works:**
1. **Build phase (Makefile):**
   - Compile `main_hnl_production.cc` → `main_hnl_production` binary
   - Handle Pythia linking
   - One-time operation

2. **Run phase (run_benchmarks.sh):**
   - Loop over mass points × flavors
   - Call `./main_hnl_production` for each
   - Collect results

**Advantages:**
- **Separation of concerns:** Build ≠ Run
- **Standard Unix:** Everyone knows `make`
- **Efficient:** Rebuild only when source changes
- **Professional:** Industry standard approach
- **Fixed physics:** Publication-quality code

---

## What Files Are Actually Used?

### For Production (Current System)

**Essential files:**
```
production/
├── main_hnl_production.cc      # Source code (C++)
├── Makefile                     # Build instructions
├── run_benchmarks.sh            # Mass scan orchestration
├── cards/                       # Pythia configuration
│   ├── hnl_Kaon.cmnd
│   ├── hnl_Dmeson.cmnd
│   ├── hnl_Bmeson.cmnd
│   └── hnl_EW.cmnd
└── pythia/pythia8315/           # Pythia library (bundled)
```

**Generated files:**
```
main_hnl_production              # Compiled binary (executable)
output/*.csv                     # Simulation results
```

**Documentation:**
```
README.md                        # User guide
PHYSICS.md                       # Physics explanation
BUGFIXES.md                      # What we fixed today
analyze_hnl_output.py           # Analysis helper script
```

**Not used anymore (can be deleted):**
```
main_hnl_single.cc              # Old source (replaced)
main_hnl_single                 # Old binary (replaced)
make.sh                         # Old build script (replaced by Makefile)
hnl_*_Template.cmnd             # Old configs (replaced by cards/)
```

---

## Typical Workflow

### One-time setup:
```bash
cd production
make                              # Compile once
```

### Run single test:
```bash
./main_hnl_production 2.6 muon 10000
```

### Run full benchmark:
```bash
./run_benchmarks.sh              # All mass points × flavors
```

### Output location:
```
production/output/HNL_*.csv      # CSV files with HNL kinematics
```

---

## Why Makefile is Better

### Old approach (make.sh):
```bash
# Compile and run mixed together
./make.sh muon
# → Compiles
# → Runs 38 mass points
# → Can't just rebuild without re-running everything
```

**Problem:** Want to just recompile after fixing a bug?
**Answer:** Have to edit make.sh to skip the mass loop. Messy!

### New approach (Makefile + script):
```bash
# Step 1: Build (once)
make

# Step 2: Test a single point
./main_hnl_production 2.6 muon 1000

# Step 3: Run full scan
./run_benchmarks.sh
```

**Want to rebuild after fixing a bug?**
```bash
make clean
make
# Done! Binary updated. No accidental runs.
```

**Want to re-run without recompiling?**
```bash
./run_benchmarks.sh
# Just runs. No compilation.
```

---

## Unix Philosophy

**"Do one thing and do it well"**

- **Makefile:** Knows how to compile C++ with Pythia
- **main_hnl_production:** Knows how to simulate one mass point
- **run_benchmarks.sh:** Knows which mass points to run

Each component has a **single responsibility**.

Old `make.sh` tried to do all three → hard to maintain.

---

## Summary Table

| Aspect | Old System (make.sh) | New System (Makefile) |
|--------|----------------------|-----------------------|
| **Build tool** | Shell script | Makefile (industry standard) |
| **Source file** | main_hnl_single.cc | main_hnl_production.cc |
| **Config files** | Template .cmnd files | cards/*.cmnd directory |
| **Orchestration** | Built into make.sh | Separate run_benchmarks.sh |
| **Physics** | Bugs (neutral 2-body) | Fixed (only physical channels) |
| **Rebuild efficiency** | Recompiles every time | Only when source changes |
| **Testability** | Hard (mixed concerns) | Easy (separate binary) |
| **Documentation** | Minimal | Complete (README/PHYSICS) |
| **Status** | ❌ Obsolete | ✅ Current |

---

## Bottom Line

**You're using the Makefile system now because:**

1. ✅ **Standard:** Everyone knows `make`
2. ✅ **Efficient:** Only rebuild when needed
3. ✅ **Flexible:** Test single points or run full scans
4. ✅ **Clean:** Separates building from running
5. ✅ **Professional:** Publication-ready code
6. ✅ **Fixed:** No more physics bugs

**The old `make.sh` system was replaced entirely.**
You don't need both - the new system is strictly better!

---

## Optional Cleanup

If you want to remove the old system files (recommended):

```bash
cd production
rm -f main_hnl_single main_hnl_single.cc make.sh
rm -f hnl_*_Template.cmnd  # Old config templates
```

Keep:
- `Makefile` ✅
- `main_hnl_production.cc` ✅
- `run_benchmarks.sh` ✅
- `cards/` directory ✅
- Everything else in current structure ✅
