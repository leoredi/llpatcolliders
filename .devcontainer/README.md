# Devcontainer Setup for LLP Physics Analysis

This devcontainer provides a complete environment for Long-Lived Particle (LLP) physics analysis with ROOT and Pythia8.

## What's Included

- **ROOT**: Latest stable version built from source with Pythia8 support
- **Pythia8**: Version 8.312 with shared libraries
- **Python 3**: With PyROOT bindings for ROOT
- **Scientific Python packages**: numpy, pandas, matplotlib, scipy, trimesh, shapely, tqdm, jupyter

## How to Use

### Option 1: Using VS Code (Recommended)

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop)
2. Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) in VS Code
3. Open this repository in VS Code
4. When prompted, click "Reopen in Container"
   - Or use Command Palette (F1) → "Dev Containers: Reopen in Container"
5. Wait for the container to build (30-60 minutes for first build, cached afterwards)

### Option 2: Using Command Line

```bash
# Build the container
cd .devcontainer
docker build -t llp-analysis .

# Run the container
docker run -it -v $(pwd)/..:/workspace llp-analysis
```

## Testing the Installation

Once the container is running, test the setup:

```bash
# Run the automated test script
bash .devcontainer/test_setup.sh
```

This will verify:
- ✓ Environment variables (PYTHIA8, PYTHIA8DATA, ROOTSYS)
- ✓ Pythia8 installation and libraries
- ✓ ROOT installation and can execute
- ✓ PyROOT (Python bindings)
- ✓ ROOT-Pythia8 integration
- ✓ Python scientific packages
- ✓ Ability to create Pythia8 objects in Python

### Manual Testing

You can also test components individually:

```bash
# Test ROOT
root --version
root -b -q -e "cout << \"ROOT works!\" << endl;"

# Test Pythia8
echo $PYTHIA8
ls -la $PYTHIA8/lib/libpythia8.so

# Test PyROOT
python3 -c "import ROOT; print('PyROOT version:', ROOT.__version__)"

# Test ROOT + Pythia8
root -b -q -e "gSystem->Load(\"libEGPythia8\"); auto p = new TPythia8(); cout << \"Pythia8 loaded!\" << endl;"

# Test Python packages
python3 -c "import numpy, pandas, matplotlib, scipy, trimesh, shapely, tqdm; print('All packages OK')"
```

## Running the Analysis

After the container is set up, you can run the physics analysis:

```bash
# Run Pythia simulation (from pythiaStuff directory)
cd pythiaStuff
./main144 -c higgsLL.cmnd

# Run analysis scripts
cd ..
python3 decayProbPerEvent.py
python3 neutral3D.py
python3 neutralv2.py

# Or use Jupyter notebooks
jupyter notebook Investigation.ipynb
```

## Environment Variables

The following environment variables are automatically set:

- `PYTHIA8=/opt/pythia8` - Pythia8 installation directory
- `PYTHIA8DATA=/opt/pythia8/share/Pythia8/xmldoc` - Pythia8 XML data files
- `ROOTSYS=/opt/root` - ROOT installation directory
- `PATH` - Includes ROOT binaries
- `LD_LIBRARY_PATH` - Includes ROOT and Pythia8 libraries
- `PYTHONPATH` - Includes ROOT Python bindings

## Troubleshooting

### Container build fails
- Ensure you have enough disk space (build requires ~10GB)
- Check Docker has enough memory allocated (recommend 8GB+)
- Check internet connection (downloads ~2GB of data)

### ROOT/Pythia8 not found after container starts
- Ensure you've opened a new terminal after container creation
- Run: `source /opt/root/bin/thisroot.sh`
- Check environment variables: `echo $ROOTSYS $PYTHIA8`

### Permission issues
- The container runs as user `vscode` (UID 1000) by default
- If you have permission issues, check your local user UID matches

## Build Times

First build (from scratch):
- Pythia8 compilation: ~5-10 minutes
- ROOT compilation: ~30-50 minutes
- Total: ~40-60 minutes

Subsequent builds:
- Docker caches layers, so rebuilds are much faster
- Only changed layers are rebuilt

## Customization

To customize the environment:

1. **Add Python packages**: Edit `Dockerfile` and add to the `pip3 install` line
2. **Change ROOT version**: Modify `--branch latest-stable` in the git clone command
3. **Change Pythia version**: Update the wget URL for Pythia8
4. **Add VS Code extensions**: Edit `devcontainer.json` under `customizations.vscode.extensions`

## Resources

- [ROOT Documentation](https://root.cern/doc/master/)
- [Pythia8 Manual](https://pythia.org/latest-manual/)
- [Dev Containers Documentation](https://code.visualstudio.com/docs/devcontainers/containers)
