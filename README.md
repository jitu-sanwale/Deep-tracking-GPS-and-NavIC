# Deep tracking of GPS and NavIC signals: IF simulator and software receiver

MATLAB code for:

- **Simulator.** A GPS L1 C/A and NavIC L5 SPS intermediate-frequency (IF) signal simulator for high-performance aircraft. It uses a nonlinear F-16 six-degree-of-freedom model with an autopilot, a navigation-grade IMU model, and an attitude-dependent antenna pattern.
- **Receiver.** A software receiver with three tracking architectures, each run with GPS alone, NavIC alone, or GPS + NavIC:
  - classical scalar tracking (**CT**);
  - vector tracking (**VT**): a vector DLL with cascaded carrier PLLs and a kinematic EKF;
  - INS-aided ultratight tracking (**UT**): the same vector DLL and cascaded PLLs with an 18-state INS error-state EKF.

The code generates all results of the paper:

> J. Sanwale, M. Kothari, and H. Hablani, *INS-Aided Deep Tracking of GPS and NavIC Signals Under Simulated Aircraft Maneuvers*.

The paper's supplementary appendices are in [`supplementary/Supplementary_Appendices.pdf`](supplementary/Supplementary_Appendices.pdf).

## Requirements

- **MATLAB:** R2023a or later, with the Signal Processing, Communications, and Curve Fitting Toolboxes.
- **Tested setup:** the paper results were produced with MATLAB R2025b on Windows 10, on a PC with 8 GB of RAM.
- **Disk space:** about **16 GB** for the full paper reproduction:
  - 9.1 GB for the IF dataset (600 MB per 60-s file and band);
  - about 6 GB for the receiver outputs.

No other software or downloads are needed. The broadcast ephemerides are included in the `Simulator` folder.

## Quick start

From the MATLAB command window, in the repository folder:

```matlab
init
```

The menu offers:

1. generate an IF dataset (simulator);
2. process an IF dataset (receiver);
3. both: generate, then process.

Both parts prompt for their choices. Press Enter to accept the default shown in brackets.

## Reproducing the paper results

The paper uses eight 60-s profiles: a level-flight baseline and seven maneuvers. They are named `baseline_cruise`, `turn`, `loop_pushover`, `break_gonset`, `jinking`, `rolls`, `splitS_immelmann`, and `departure_spin`. The receiver processes each in three configurations (GPS alone, NavIC alone, GPS + NavIC) with all three architectures, which makes 24 receiver runs.

### Interactive

Run `init` and choose `3`.

- **Simulator:** choose `9` (all eight predefined maneuvers) and both bands.
- **Receiver:** choose `all of them`, configuration `4` (all three), and receiver modes `4` (CT, VT, and UT).

### Non-interactive (recommended for the long run)

Set the data and output folders, then run each part with its settings struct. With these structs defined, neither part asks any questions.

```matlab
data = 'D:\deeptracking\data';            % IF dataset, truth, and IMU files
out  = 'D:\deeptracking\receiver_output'; % receiver results

SIM = struct('select', {{'baseline_cruise','turn','loop_pushover','break_gonset', ...
                         'jinking','rolls','splitS_immelmann','departure_spin'}}, ...
             'gnss', 'both', 'seconds', 60, 'dataDir', data);
run(fullfile('Simulator', 'init.m'));

RX = struct('dataDir', data, 'outDir', out, 'select', 'all', ...
            'config', {{'GPS', 'NavIC', 'GPS+NavIC'}}, 'mode', 4);
run(fullfile('Receiver', 'init.m'));
```

Or run it from a terminal without the desktop, for example by saving the lines above as `run_all.m`:

```
matlab -batch "run_all"
```

### Run time

These times were measured on the test PC:

| Stage | Time |
|---|---|
| Simulator: 8 profiles, GPS + NavIC IF files | about 2.5 h |
| Receiver: GPS alone, per profile | about 15 min |
| Receiver: NavIC alone, per profile | about 15 min |
| Receiver: GPS + NavIC, per profile | about 35–40 min |
| **Full reproduction** | **about 13 h** |

The receiver can be interrupted and restarted. Completed runs are kept; a run is complete when it has `truth_comparison.mat`. On restart the receiver continues with the remaining runs. To recompute everything, set `RX.rerun = true`.

### Where the paper's tables and figures come from

When all 24 runs are available, the receiver writes the full report to `<outDir>\report_all_cases\`:

| Paper | File in `report_all_cases` |
|---|---|
| Table 3, NavIC satellites and geometry | `table_navic_satellites.tex` |
| Figure 2, DOP values | `geometry_dop.txt` |
| Figure 3, antenna gain vs. attitude | `antenna_attitude.png` |
| Table 5, maneuver profiles | `table_profiles.tex` |
| Figure 4, load factor and roll rate | `profiles_dynamics.png` |
| Table 6, NavIC L5 alone | `table_navic_alone.tex` |
| Figure 5, error comparison (log scale) | `navic_compare.png` |
| Figure 6, NavIC-only turn errors | `turn_error.png` |
| Table 7, GPS L1 alone | `table_gps_alone.tex` |
| Figure 7, UT filter diagnostics (break) | `break_ut_diag.png` |
| Table 8, NavIC L5 + GPS L1 and inter-system bias | `table_navic_gps.tex` |
| Statistics quoted in the text | `results_statistics.txt` |

**Not in the report:**
- Tables 2 and 4 list fixed IMU and receiver parameters. The values are set in `Simulator/imuErrorModel_iNAT.m`, `Simulator/generateIMU_iNAT_FSLG01.m`, and `Receiver/initParameters.m`.
- The sky plot of Figure 2 is drawn by the simulator (`Simulator/skyPlot.m`).
- Figures 6 and 7 in the paper are the same data redrawn at print size.

**Statistics window.** Errors are evaluated against the simulator truth, in each maneuver window from the maneuver start (no earlier than 3.1 s) to 1 s after its end (no later than 59.7 s). For the baseline the window is 3.1–59.7 s. The window times are listed in Supplementary Appendix C.

**Bitwise reproducibility.** The simulator uses fixed random seeds. With the same MATLAB release, a fresh run regenerates the maneuver IF files bit for bit. The receiver is deterministic.

## Using the code for your own cases

### Simulator (`Simulator/init.m`)

**Choices:**
- `1`–`8`: a predefined maneuver.
- `9`: all eight predefined maneuvers.
- `10`: a custom maneuver built from blocks (turn, bank, roll, pull/push, loop, and others) and flown by the same 6-DOF model.
- `11`: your own trajectory from a CSV file.

Several choices can be entered at once, for example `2 4 10`.

Non-interactive examples:

```matlab
SIM.select = {'turn','break_gonset'};  SIM.gnss = 'both';  run('Simulator/init.m')
SIM.select = {'custom_myturn'};  SIM.custom.custom_myturn = {{'turn',6,8,3,1},{'roll',360,270}};  run('Simulator/init.m')
SIM.select = {'custom_track'};   SIM.csv.custom_track = 'C:\data\track.csv';  run('Simulator/init.m')
```

`SIM` fields:

| Field | Meaning |
|---|---|
| `select` | profile names |
| `gnss` | `'GPS'`, `'NavIC'`, or `'both'` |
| `seconds` | record length; default 60 |
| `overwrite` | replace existing files; default false |
| `dataDir` | output folder |
| `custom.<name>` | block list for a custom maneuver |
| `csv.<name>` | CSV trajectory file |
| `description.<name>` | description of the profile |

The header of `Simulator/customManeuverBlocks.m` documents the block list. `Simulator/importCustomTrajectory.m` documents the CSV format.

**Output in the data folder, for each profile `<name>`:**

| File | Content |
|---|---|
| `FighterProfiles\Fighter_<name>.csv` | 1-kHz truth trajectory |
| `FighterProfiles\Fighter_<name>_IMU_200Hz.mat` (and `.csv`) | 200-Hz IMU record with truth |
| `IFdata\gps_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin` | GPS L1 C/A IF |
| `IFdata\navic_0IF_Fs5MHz_8bit_IQ_Fighter_<name>.bin` | NavIC L5 IF |

The IF files are 5-MHz complex baseband, stored as 8-bit I/Q.

### Receiver (`Receiver/init.m`)

The receiver reads any data folder with the layout above. The GPS IF file is its timing reference, so it is needed in every configuration, including NavIC alone.

**Selection:**
- one number processes a single case;
- several numbers, for example `2 4`, process a combination;
- `all of them`, with all three configurations, also produces the full paper report.

`RX` fields:

| Field | Meaning |
|---|---|
| `dataDir` | input data folder |
| `outDir` | output folder; default `Receiver\output` |
| `select` | profile names or `'all'` |
| `config` | a subset of `{'GPS','NavIC','GPS+NavIC'}` |
| `mode` | `1` CT, `2` VT, `3` UT, `4` all; default 4 |
| `rerun` | recompute existing results; default false |

```matlab
RX.select = {'turn'};  RX.config = {'NavIC'};  run('Receiver/init.m')
```

**Outputs:**
- **Every case,** in `<outDir>\results_<dataset>\`:
  - plots: positions, velocities, C/N0, clock, errors against truth, and UT diagnostics;
  - error CSVs at 10 Hz;
  - `truth_comparison.mat`;
  - the case results table, `case_table.md` and `case_table.tex`, also printed in the command window.
- **Several cases,** in `<outDir>\results_FighterSummary*\`: a cross-maneuver summary table and bar chart.
- **All cases:** the full report described above.

**Warm start.** The maneuvers share the first 12 s of the baseline flight. Their runs therefore reuse the baseline's decoded ephemeris (a warm start). If `baseline_cruise` is in the data folder, it is processed first, and it is added to the selection automatically when needed.

## Repository layout

| Path | Content |
|---|---|
| `init.m` | entry menu: simulate, process, or both |
| `Simulator/` | IF signal generation, F-16 model and autopilot (`F16Aircraft/`), IMU model, antenna patterns, broadcast ephemerides (RINEX) |
| `Receiver/` | receiver entry and parameters |
| `Receiver/acqtckpos/` | acquisition, CT/VT/UT tracking and positioning, navigation-message decoding, truth comparison |
| `Receiver/geo/` | coordinate, Earth-rotation, and ionosphere functions |
| `Receiver/common/` | functions shared with the simulator (antenna, troposphere, trajectory spline) |
| `Receiver/report/` | case tables and the full paper report |
| `supplementary/` | supplementary appendices of the paper |

## Credits and licence

**Receiver.** Builds on GPSSDR_vt by B. Xu and L.-T. Hsu (*GPS Solutions* 23, 46, 2019).

**Simulator.** Builds on GPS_L1CA_signal_generator, the open-source IF signal simulator of Y. Li (Beijing Information Science and Technology University), as modified by J. Sun. Both upstream codes have been extended substantially for NavIC.

**Other sources:**
- Both upstream codes include routines from SoftGNSS (K. Borre, D. Akos, D. Plausinaitis, and others).
- The F-16 aerodynamic data are from F16Model.jl, under the MIT licence; see `Simulator/F16Aircraft/LICENSE_F16Model_jl.txt`.
- The broadcast ephemerides are International GNSS Service (IGS) data, from https://igs.org.

The original notices are kept in the file headers.

This code is released under the GNU General Public License v3.0 or later; see [`LICENSE`](LICENSE).
