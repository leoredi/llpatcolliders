#!/bin/bash
# Test script to validate devcontainer setup with ROOT and Pythia8

set -e  # Exit on error

echo "=========================================="
echo "Testing Devcontainer Setup"
echo "=========================================="
echo ""

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counter
TESTS_PASSED=0
TESTS_FAILED=0

# Function to print test results
test_result() {
    if [ $1 -eq 0 ]; then
        echo -e "${GREEN}✓ PASSED${NC}: $2"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ FAILED${NC}: $2"
        ((TESTS_FAILED++))
    fi
    echo ""
}

# Test 1: Check environment variables
echo "Test 1: Checking environment variables..."
if [ -n "$PYTHIA8" ] && [ -n "$PYTHIA8DATA" ] && [ -n "$ROOTSYS" ]; then
    echo "  PYTHIA8: $PYTHIA8"
    echo "  PYTHIA8DATA: $PYTHIA8DATA"
    echo "  ROOTSYS: $ROOTSYS"
    test_result 0 "Environment variables are set"
else
    test_result 1 "Environment variables are missing"
fi

# Test 2: Check Pythia8 installation
echo "Test 2: Checking Pythia8 installation..."
if [ -f "$PYTHIA8/lib/libpythia8.so" ]; then
    echo "  Pythia8 library found at: $PYTHIA8/lib/libpythia8.so"
    test_result 0 "Pythia8 library exists"
else
    test_result 1 "Pythia8 library not found"
fi

# Test 3: Check Pythia8 data files
echo "Test 3: Checking Pythia8 data files..."
if [ -d "$PYTHIA8DATA" ] && [ -f "$PYTHIA8DATA/Index.xml" ]; then
    echo "  Pythia8 data directory: $PYTHIA8DATA"
    test_result 0 "Pythia8 data files exist"
else
    test_result 1 "Pythia8 data files not found"
fi

# Test 4: Check ROOT installation
echo "Test 4: Checking ROOT installation..."
if command -v root &> /dev/null; then
    ROOT_VERSION=$(root --version 2>&1 | head -n 1)
    echo "  $ROOT_VERSION"
    test_result 0 "ROOT is installed and in PATH"
else
    test_result 1 "ROOT command not found"
fi

# Test 5: Check ROOT libraries
echo "Test 5: Checking ROOT libraries..."
if [ -f "$ROOTSYS/lib/libCore.so" ]; then
    echo "  ROOT libraries found in: $ROOTSYS/lib/"
    test_result 0 "ROOT libraries exist"
else
    test_result 1 "ROOT libraries not found"
fi

# Test 6: Test ROOT can start
echo "Test 6: Testing ROOT executable..."
if root -b -q -e "cout << \"ROOT works!\" << endl;" &> /dev/null; then
    test_result 0 "ROOT runs successfully"
else
    test_result 1 "ROOT failed to run"
fi

# Test 7: Check Python ROOT bindings (PyROOT)
echo "Test 7: Checking PyROOT (Python ROOT bindings)..."
if python3 -c "import ROOT; print('PyROOT version:', ROOT.__version__)" 2>&1 | grep -q "PyROOT"; then
    PYROOT_VERSION=$(python3 -c "import ROOT; print(ROOT.__version__)" 2>&1)
    echo "  PyROOT version: $PYROOT_VERSION"
    test_result 0 "PyROOT is available"
else
    test_result 1 "PyROOT import failed"
fi

# Test 8: Check ROOT Pythia8 support
echo "Test 8: Checking ROOT Pythia8 integration..."
TEST_OUTPUT=$(root -b -q -e "gSystem->Load(\"libEGPythia8\"); cout << \"Pythia8 library loaded\" << endl;" 2>&1)
if echo "$TEST_OUTPUT" | grep -q "Pythia8 library loaded"; then
    echo "  ROOT can load Pythia8 library"
    test_result 0 "ROOT has Pythia8 support"
else
    echo "  Output: $TEST_OUTPUT"
    test_result 1 "ROOT cannot load Pythia8 library"
fi

# Test 9: Check Python packages
echo "Test 9: Checking Python scientific packages..."
MISSING_PACKAGES=""
for pkg in numpy pandas matplotlib scipy trimesh shapely tqdm jupyter; do
    if ! python3 -c "import $pkg" 2>/dev/null; then
        MISSING_PACKAGES="$MISSING_PACKAGES $pkg"
    fi
done

if [ -z "$MISSING_PACKAGES" ]; then
    echo "  All required packages installed: numpy, pandas, matplotlib, scipy, trimesh, shapely, tqdm, jupyter"
    test_result 0 "Python packages installed"
else
    echo "  Missing packages:$MISSING_PACKAGES"
    test_result 1 "Some Python packages missing"
fi

# Test 10: Simple Pythia8 test via PyROOT
echo "Test 10: Testing Pythia8 via PyROOT..."
cat > /tmp/test_pythia.py << 'PYEOF'
import ROOT
try:
    pythia = ROOT.TPythia8()
    print("Pythia8 object created successfully")
    exit(0)
except Exception as e:
    print(f"Failed to create Pythia8 object: {e}")
    exit(1)
PYEOF

if python3 /tmp/test_pythia.py 2>&1 | grep -q "successfully"; then
    test_result 0 "Can create Pythia8 objects in Python"
else
    test_result 1 "Cannot create Pythia8 objects in Python"
fi

# Summary
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo -e "${GREEN}Passed: $TESTS_PASSED${NC}"
echo -e "${RED}Failed: $TESTS_FAILED${NC}"
echo "=========================================="

if [ $TESTS_FAILED -eq 0 ]; then
    echo -e "${GREEN}All tests passed! ✓${NC}"
    echo "Your devcontainer is ready for LLP physics analysis!"
    exit 0
else
    echo -e "${YELLOW}Some tests failed. Please check the output above.${NC}"
    exit 1
fi
