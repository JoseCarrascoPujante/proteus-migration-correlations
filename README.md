# Systemic and kinematic properties of *Amoeba proteus* migration correlate strongly across environments: code and data

Code and data for the article

> Carrasco-Pujante J, Malaina I, Pérez-Yarza G, Camino-Pontes B, Cortés JM, López JI, De la Fuente IM.
> Systemic and kinematic properties of *Amoeba proteus* migration correlate strongly across environments.

Running `main` in MATLAB reproduces every result and table of the article from the cell trajectories; its last step calls R for the two-lines test (see Requirements).

## Data

`data/coordinates.mat` contains the struct `coordinates`, with one field per experimental scenario
(`NoStimuli`, `Galvanotaxis`, `Chemotaxis`). Each holds two 3600 × 50 matrices, `x` and `y`: the
position, in mm, of each of the 50 cells every 0.5 s for 30 min, relative to its initial position.

## Requirements

- MATLAB R2026a with the Statistics and Machine Learning and Curve Fitting Toolboxes. With the
  Parallel Computing Toolbox, three of the loops run in parallel, with the same results.
- R (4.6.1 was used) with Uri Simonsohn's package statuser (0.4.0 was used), whose `twolines()` runs
  the two-lines test of `main`'s last step: `install.packages("statuser")`.

## Usage

From this folder, in MATLAB, run `main` (about 5 min with two parallel workers). The functions are in the MATLAB namespace folder `+corrs_proteus/` and are called as `corrs_proteus.<function>`; only this folder's parent needs to be on the path, which `main` sets. The results are written to `data/derived/`, where the repository already holds those of the run behind the article:

| Step | Output | In the article |
|---|---|---|
| 1 to 3 | `results.mat`, `lillieforsTests.mat` | The eleven metrics; normality tests |
| 4 | `pairwiseDependence.csv` | Table S13 |
| 5 | `stimulusEffectSizes.csv` | Table S12 |
| 6 | `dependenceContrasts.csv` | Table S15 |
| 7 to 9 | `GoodnessOfFitResults.mat` and `.xlsx`, `aicBorda_*.csv`, `perContestDecision_exp1.csv` | Candidate forms and their AICc ranking |
| 10 | `table1_modeledCouplings.csv` | Table 1 |
| 11 | `aiccBordaRanking.csv` | Table S14 |
| 12 | `perTrajectory_<metric>.csv` | Tables S1 to S11 |
| 13 | Printed in the command window | Departure of DFA γ from β/2 + 1 (Results) |
| 14 | Printed in the command window | Robustness of the every-frame metrics (Materials and Methods) |
| 15 | Printed in the command window | Side of the net displacement of every cell (Results) |
| 16 | Printed in the command window | Speed memory against temporally independent speeds (Results) |
| 17 | `samplingIntervalSweep.csv` | Fig S4 and the values the Results quote from it |
| 18 | `pcaComponents.csv`, `pcaMetrics.csv`, `participationRatio.csv` | Fig S7 |
| 19 | `dependenceDendrogram.csv` | Fig S8 |
| 20 | `twolines_check.csv` | Two-lines test |

The permutation tests and bootstrap resamples use fixed seeds, so every run gives the same results.

## License

Code under the MIT License, data under CC BY 4.0 (see `LICENSE`).

## Citation

Please cite the article; citation metadata for this repository are in `CITATION.cff`.
