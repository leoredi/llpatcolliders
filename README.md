# llpatcolliders

HNL sensitivity projections for a transverse LLP detector near CMS at the HL-LHC.

## Run

```bash
conda env create -f environment.yml
conda activate llpatcolliders

cd production/pythia_production
PYTHIA8=$(pwd)/pythia8315 make
./run_parallel_production.sh all both
./run_parallel_production.sh all direct hardccbar 10
./run_parallel_production.sh all direct hardbbbar 10
./run_parallel_production.sh all direct hardBc 15

# optional EW production (inside Docker image: mg5-hnl)
# see production/madgraph_production/

cd ../../analysis_pbc
python limits/combine_production_channels.py
python limits/run.py --parallel --workers 12

cd ../
python money_plot/plot_money_island.py
```

## Read next

- Physics assumptions and formulas: `PHYSICS.md`
- Code structure and execution details: `CODING.md`
