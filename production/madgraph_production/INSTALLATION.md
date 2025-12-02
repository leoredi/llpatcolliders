# MadGraph HNL Pipeline - Installation & Setup

This guide covers building and testing the Docker-based MadGraph pipeline for HNL electroweak production.

---

## Prerequisites

**On your host machine, you need:**
- Docker Desktop (macOS/Windows) or Docker Engine (Linux)
- Git (to clone the repository)
- ~2 GB disk space for Docker image

**That's it!** Everything else (MadGraph, Python, Fortran compilers) runs inside the container.

---

## Installation Steps

### 1. Clone Repository (If Not Already Done)

```bash
cd /path/to/your/projects
git clone <your-repo-url> llpatcolliders
cd llpatcolliders
```

### 2. Build Docker Image

```bash
cd production/madgraph
docker build -t mg5-hnl .
```

**What this does:**
- Downloads Ubuntu 22.04 base image
- Installs Python 3, Fortran compiler (gfortran), build tools
- Extracts MadGraph 3.6.6 from local tarball (`MG5_aMC_v3.6.6.tar.gz`)
- Installs Pythia8 and LHAPDF6 via MadGraph
- Installs Python packages: numpy, pandas, pylhe
- Sets up environment variables

**Build time**: ~5-10 minutes (first time only)

**Expected output:**
```
[+] Building 350.2s (10/10) FINISHED
 => [1/6] FROM docker.io/library/ubuntu:22.04
 => [2/6] RUN apt-get update && apt-get install...
 => [3/6] COPY MG5_aMC_v3.6.6.tar.gz /opt/
 => [4/6] RUN tar -xzf /opt/MG5_aMC_v3.6.6.tar.gz...
 => [5/6] RUN pip3 install --no-cache-dir...
 => [6/6] RUN echo "install pythia8..." && ...
 => exporting to image
 => => naming to docker.io/library/mg5-hnl:latest
```

### 3. Verify Image

```bash
docker images | grep mg5-hnl
```

**Expected output:**
```
mg5-hnl    latest    <image-id>   5 minutes ago   1.52GB
```

---

## Quick Test

### 1. Launch Container

```bash
# From llpatcolliders repository root:
docker run --rm -it -v "$(pwd)":/work mg5-hnl bash
```

**Explanation:**
- `--rm`: Remove container when it exits (clean up)
- `-it`: Interactive terminal
- `-v "$(pwd)":/work`: Mount current directory (`llpatcolliders/`) into container at `/work`
- `mg5-hnl`: Image name
- `bash`: Start bash shell

You should see:
```
root@<container-id>:/work#
```

### 2. Verify MadGraph

Inside the container:

```bash
# Check MadGraph version
/opt/MG5_aMC_v3_6_6/bin/mg5_aMC --version

# Expected output:
# MG5_aMC 3.6.6
```

### 3. Test HeavyN Model

Still inside container:

```bash
cd /work/production/madgraph

# Start MadGraph interactive
/opt/MG5_aMC_v3_6_6/bin/mg5_aMC
```

In the MadGraph prompt:
```
MG5_aMC> import model SM_HeavyN_CKM_AllMasses_LO
MG5_aMC> display particles
```

**Expected output** (should show):
```
...
n1 (PDG: 9900012)
n2 (PDG: 9900014)
n3 (PDG: 9900016)
...
```

Type `exit` to quit MadGraph.

### 4. Run Test Point

Inside container, from `/work/production/madgraph`:

```bash
python3 scripts/run_hnl_scan.py --test
```

**Expected behavior:**
- Creates process directory `work/hnl_muon_15.0GeV/`
- Generates 1000 events (takes ~2-5 minutes)
- Converts LHE → CSV
- Outputs CSV to `csv/muon/HNL_mass_15.0GeV_muon_EW.csv`

**Expected log snippets:**
```
[INFO] Starting MadGraph HNL Production Pipeline
[INFO] Test mode: 15 GeV muon, 1000 events
✓ Process directory created
✓ Cards written
✓ Events generated
✓ Cross-section extracted: 45.3 ± 0.2 pb
✓ Converted 1000 events to CSV
✓ Updated summary CSV
```

### 5. Verify Output

Still in container:

```bash
# Check CSV was created
ls -lh csv/muon/

# Expected:
# HNL_mass_15.0GeV_muon_EW.csv  (~100 KB for 1000 events)

# Check first few lines
head -5 csv/muon/HNL_mass_15.0GeV_muon_EW.csv

# Expected:
# event_id,parent_pdgid,hnl_pdgid,mass_hnl_GeV,weight,...
# 1,24,9900012,15.0,1.234e-03,...
# ...
```

### 6. Exit Container

```bash
exit
```

The container is automatically removed (`--rm` flag), but **files in the mounted volume persist** on your host machine.

Check on your host:
```bash
# From llpatcolliders/
ls production/madgraph/csv/muon/

# Should show:
# HNL_mass_15.0GeV_muon_EW.csv
```

---

## What's Installed in the Docker Image?

### System Packages
- **OS**: Ubuntu 22.04 LTS
- **Compilers**: gcc, g++, gfortran (for Fortran matrix elements)
- **Build tools**: make, wget, curl
- **Python**: Python 3.10

### MadGraph Setup
- **Location**: `/opt/MG5_aMC_v3_6_6/`
- **Version**: 3.6.6
- **Executable**: `/opt/MG5_aMC_v3_6_6/bin/mg5_aMC`
- **UFO Model**: `SM_HeavyN_CKM_AllMasses_LO` (included in image)

### Additional Tools (Installed via MadGraph)
- **Pythia8**: Event generator (parton shower, hadronization)
- **LHAPDF6**: Parton distribution functions

### Python Packages
- **numpy**: Numerical arrays
- **pandas**: CSV handling
- **pylhe**: LHE file parsing

### Environment Variables (set in container)
```bash
MG5_DIR=/opt/MG5_aMC_v3_6_6
MG5_PATH=/opt/MG5_aMC_v3_6_6/bin/mg5_aMC
PATH=/opt/MG5_aMC_v3_6_6/bin:$PATH
```

---

## Troubleshooting Installation

### Issue 1: "docker: command not found"

**Cause**: Docker not installed on host

**Solution**:
- **macOS**: Install Docker Desktop from https://www.docker.com/products/docker-desktop
- **Linux**:
  ```bash
  curl -fsSL https://get.docker.com | sh
  sudo usermod -aG docker $USER  # Allow non-root access
  # Log out and back in
  ```
- **Windows**: Install Docker Desktop (requires WSL 2)

### Issue 2: "permission denied while trying to connect to Docker daemon"

**Cause**: User not in `docker` group (Linux)

**Solution**:
```bash
sudo usermod -aG docker $USER
# Log out and back in
```

### Issue 3: "no space left on device"

**Cause**: Not enough disk space for Docker image

**Solution**:
```bash
# Check Docker disk usage
docker system df

# Clean up old images/containers
docker system prune -a

# Requires ~2 GB free space
```

### Issue 4: Build fails at "install pythia8"

**Symptom**:
```
ERROR: Could not install pythia8
```

**Cause**: Network issue downloading Pythia8 from external server

**Solution**:
```bash
# Retry build (MadGraph will resume download)
docker build -t mg5-hnl .

# If persistent, check network connection
# Pythia8 is downloaded from http://home.thep.lu.se/~torbjorn/pythia8/
```

### Issue 5: "MG5_aMC_v3.6.6.tar.gz: no such file"

**Cause**: Building from wrong directory

**Solution**:
```bash
# MUST build from production/madgraph/ where Dockerfile lives
cd llpatcolliders/production/madgraph
docker build -t mg5-hnl .
```

### Issue 6: Container starts but `/work` is empty

**Symptom**: `ls /work` shows nothing inside container

**Cause**: Wrong directory mounted

**Solution**:
```bash
# Run from llpatcolliders repository root, NOT from production/madgraph
cd /path/to/llpatcolliders
docker run --rm -it -v "$(pwd)":/work mg5-hnl bash

# Inside container, you should see:
ls /work/production/madgraph
```

---

## Updating the Pipeline

If you modify Python scripts (`run_hnl_scan.py`, `lhe_to_csv.py`) or card templates, you **don't need to rebuild** the Docker image. Just re-run the container:

```bash
# From llpatcolliders/
docker run --rm -it -v "$(pwd)":/work mg5-hnl bash

# Your changes are immediately visible inside /work
```

**Rebuild only if** you change:
- `Dockerfile` (e.g., different MadGraph version)
- System dependencies (e.g., new apt packages)
- UFO model installation

---

## Advanced: Running Without Interactive Shell

For automated workflows (e.g., on a cluster), you can run commands directly:

```bash
# From llpatcolliders repository root:

# Test run
docker run --rm -v "$(pwd)":/work mg5-hnl \
  bash -c "cd /work/production/madgraph && python3 scripts/run_hnl_scan.py --test"

# Full muon scan
docker run --rm -v "$(pwd)":/work mg5-hnl \
  bash -c "cd /work/production/madgraph && python3 scripts/run_hnl_scan.py --flavour muon"
```

**Tip**: For long runs, use `nohup` or `screen` on the host to keep the container running if you disconnect:

```bash
nohup docker run --rm -v "$(pwd)":/work mg5-hnl \
  bash -c "cd /work/production/madgraph && python3 scripts/run_hnl_scan.py" \
  > madgraph_run.log 2>&1 &
```

---

## Next Steps

After successful installation and test:

1. **Read main README**: See `README.md` for usage examples
2. **Understand physics**: Check physics details section
3. **Run production scan**: Start with single flavour (`--flavour muon`)
4. **Integrate with analysis**: See `../../analysis_pbc_test/` for downstream pipeline

---

## Summary Checklist

Installation is complete when you can check all boxes:

- [ ] Docker installed and running
- [ ] `mg5-hnl` image built successfully (~1.5 GB)
- [ ] Container starts with `docker run --rm -it -v "$(pwd)":/work mg5-hnl bash`
- [ ] MadGraph responds to `--version` flag
- [ ] HeavyN model loads (`import model SM_HeavyN_CKM_AllMasses_LO`)
- [ ] Test run completes (`python3 scripts/run_hnl_scan.py --test`)
- [ ] CSV output appears in `production/madgraph/csv/muon/`

**If all boxes checked**: ✅ Installation successful! Proceed to production runs.

---

## Getting Help

**Installation issues:**
- Check Docker logs: `docker logs <container-id>`
- Check disk space: `df -h`
- Check Docker version: `docker --version` (requires 20.10+)

**MadGraph issues:**
- Check MadGraph log: `cat work/hnl_muon_15.0GeV/Events/run_01/*.log`
- Test MadGraph directly: `/opt/MG5_aMC_v3_6_6/bin/mg5_aMC --help`

**For other issues**: See main project documentation at `../../CLAUDE.md`
