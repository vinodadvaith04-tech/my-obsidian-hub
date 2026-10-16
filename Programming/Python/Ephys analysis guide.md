# Ephys analysis: how the script works and what the readouts mean

Oct 14, 2026 · @Adv

## Start here

This guide has two parts: how `ephys_analysis.py` works section by section (Part 1), and what every number it reports means (Part 2).

The script does one job. It opens every `.abf` file, works out from the protocol name stored *inside* the file what kind of recording it is, measures that recording, and writes the results as CSV tables, one quick-look PNG per recording, and one summary row per cell.

Two habits make editing safe:

- Copy the file first (`ephys_analysis_v1.py`) and edit the copy.
- Test changes on a single day in seconds: `python ephys_analysis.py "<data folder>" --only 20261007 --no-plots --out test_out`.

In Part 1 I refer to sections by function name, not line number. Use your editor's search (Ctrl+F) for the name in backticks.

## Part 1: How the script works

### The flow, from files to results

The script reads top to bottom as a library of functions, and `main()` at the bottom calls them in this order:

1. Collect every `.abf` under the data folder (ignoring `_analysis`).
2. For each file, open it (`Rec`), and skip it if an identical copy was already seen (same recording time and file size).
3. Work out the recording type from the protocol name in the header (`classify`).
4. Run the matching analysis function and write its numbers into a table row.
5. Draw the quick-look PNG for that recording.
6. Group files into cells using the cell part of the new file name.
7. Build one summary row per cell (first cap check, first HCN, first firing run, and so on), and flag problems.
8. Merge your metadata sheet if you gave one, then write all CSVs and `qc_log.txt`.

### Map of the file

The file is ordered the same way you would read it. Each row is a block you can find by searching for its name.

| Block | Functions | What it does |
| --- | --- | --- |
| Settings | `CFG`, colours, `Vc` | All thresholds and windows; the LJP correction |
| Reading files | `classify`, `Skip`, `need_both`, `check_cc`, `Rec`, `find_edges`, `smooth` | Open a file, identify its channels, refuse unreliable ones |
| Cap check | `analyze_capcheck`, `_fit_transient`, `plot_capcheck` | Ra, Rm, Cm, tau, holding current |
| AP detection | `detect_aps` | Finds every action potential and measures it |
| Step helpers | `step_geometry`, `sweep_step_levels`, `base_ss` | Find the current step and its baseline and steady-state voltages |
| HCN, current clamp | `analyze_hcn_cc`, `plot_hcn_cc` | Sag, rebound, Rin, tau\_m |
| HCN, voltage clamp | `analyze_hcn_vc`, `plot_hcn_vc` | Ih amplitude and tail current |
| Firing | `analyze_firing`, `plot_firing` | Rheobase, AP features, f-I, adaptation, Rin |
| AHP protocols | `analyze_ahp`, `plot_ahp` | AHP after a single evoked AP |
| Gap-free | `analyze_gapfree`, `plot_gapfree` | Basic QC only |
| Bookkeeping | `NAME_RE`, `parse_name`, `write_csv` | File names to cell keys, writing CSVs |
| `main()` | the loop and the summary block | Ties everything together |

### The top of the file: settings

Most tweaks you will want live in one place, the `CFG` dictionary near the top. Changing a number there changes it everywhere the script uses it.

| Setting | Default | What it controls |
| --- | --- | --- |
| `LJP_MV` | 15 | Subtracted from every absolute voltage (RMP, threshold, AP peak). Overridden by `--ljp` on the command line. |
| `SS_WINDOW_S` | 0.100 | Steady-state voltage = mean of the last 100 ms of a current step. |
| `BASE_WINDOW_S` | 0.100 | Baseline Vm = mean of the 100 ms before a step. |
| `DVDT_THRESH` | 20 | AP threshold criterion in mV/ms. Higher (30, 50) gives a more depolarised threshold. |
| `AP_MIN_PEAK_MV` | -10 | An event counts as an AP only if it peaks above this (recorded mV). |
| `AP_MIN_AMP_MV` | 25 | And only if it rises at least this far above its threshold. Filters noise and artefacts. |
| `REBOUND_WINDOW_S` | 0.5 | How long after the HCN step the rebound is searched for. |
| `RIN_MAX_DEPOL_MV` | 15 | Sub-threshold steps used for Rin must stay within this many mV of baseline. |
| `RA_MAX_MOHM` | 30 | Cap checks above this Ra are flagged. |
| `RA_DRIFT_MAX_PCT` | 20 | Flag if Ra changed by more than this between a cell's first and last cap check. |
| `RM_MIN_MOHM` | 100 | Flag a very leaky cell. |
| `IHOLD_MAX_ABS_PA` | 200 | Flag a large holding current at the cap check. |
| `BIAS_MAX_ABS_PA` | 50 | Flag cells that had current injected at baseline in current clamp. |

Flags never delete data. They only add text to `qc_flags` in the summary.

### Liquid junction potential: `Vc()`

`Vc(v)` returns `v - LJP`. It is applied to absolute voltages only: RMP, steady state, peak, AP threshold, AP peak, rebound peak. Differences such as sag amplitude, AP amplitude and half-width are unaffected, so they are the same in Clampfit.

If you want to change the sign or apply it differently, this one function is the only place to do it.

### Colours and plot style

`C1`, `C2`, `C3` are the three plot colours (blue, orange, green). `INK`, `GRID` and `SURF` are text, gridline and background colours. `plt.rcParams.update(...)` sets fonts and the clean look. Change these to restyle every PNG at once.

### Numbers that are not in `CFG`

Some windows are written directly inside the functions. Search for them if you want to change them.

- Cap check: baseline = 12 ms to 1 ms before the step; steady state = last 5 ms of the step; transient peak searched in the first 2 ms.
- HCN voltage clamp: instantaneous current = 8 to 13 ms after step onset; tail current = 8 to 28 ms after step offset.
- AP measures: fAHP window 5 ms and AHP window 50 ms after the peak; half-width search 6 ms after the peak.
- Firing: depolarisation block rule, adaptation reference sweep (at least 5 spikes), f-I slope uses the first 6 spiking sweeps.
- AHP protocols: fAHP 0 to 5 ms, mAHP 5 to 100 ms, sAHP 0.1 to 1.5 s after the AP peak.

### Reading files and deciding what they are

This block turns an `.abf` file into numpy arrays and decides whether it can be trusted. Everything else depends on it.

**`classify(protocol)`** reads the protocol name saved inside the file (for example `CC_spikefreq_0-300pA_25pAsteps`) and returns a label: `CapCheck`, `GapFree`, `HCN-CC`, `HCN-VC`, `Firing25pA`, `Firing10pA`, `AHP` or `Other`. To add a new protocol, add a case here and a branch in `main()`.

**`Rec`** is the class that holds one file. It reads the data into arrays shaped (sweeps, samples), stores the sampling rate `fs`, a time axis `t`, and the recording time. It then picks the voltage channel (`r.V`) and current channel (`r.I`) by their **units**, not by their position. That matters because the channel order differs between your rigs and sessions.

**`Skip`, `need_both`, `check_cc`** are the safety checks. Your rigs sometimes saved files with both channels labelled pA, or with a flat current channel. Measuring those would produce nonsense, so the script raises `Skip`. The file is listed in `qc_log.txt` as SKIPPED and the script carries on.

- `need_both` requires one mV channel and one pA channel.
- `check_cc` additionally requires that Vm looks like a membrane potential (between -200 and +100 mV, median between -130 and +20 mV, some movement) and that the current channel really steps (no more than 2000 pA).
- The cap check skips itself if the command step is under 2 mV.

**`find_edges(x, fs)`** finds the start and end sample of the biggest pulse in a trace. It looks for where the signal departs from its starting level by more than half of its maximum deviation. Nearly every protocol uses it to find where the step happens, so the code never assumes a fixed time.

**`smooth(x, n)`** is a running average of `n` samples. It is used where a single noisy sample could give a false minimum, for example the sag peak (1 ms) and the AHP trough.

**`nominal(x, unit)`** rounds a measured current to the nearest protocol step (25, 10 or 100 pA), so a measured 25.2 pA is reported as 25. Both values are kept in the tables (`step_pA` and `step_pA_meas`).

### Why units, not channel numbers?

In some of your recordings the amplifier labelled both channels pA, or the order was swapped. Choosing channels by units, and skipping when they do not make sense, is what keeps those files from silently contaminating your averages.

### The analysis functions, one by one

Every `analyze_*` function takes a `Rec`, returns plain numbers, and does not draw or write files. The matching `plot_*` function draws the PNG. This split means you can change a measurement without touching a plot, and the reverse.

#### Cap check: `analyze_capcheck`

1. Average the 10 sweeps, which cuts the noise by about 3.
2. Find the 5 mV command step with `find_edges` on the command channel.
3. Measure the baseline current (12 ms to 1 ms before the step). That is the holding current, `Ihold_pA`.
4. `_fit_transient` fits `I(t) = Iss + A·exp(-t/tau)` from the peak of the capacitive spike to the end of the step, using scipy's `curve_fit`. It then extrapolates back to the exact moment of the step to get the current at t=0.
5. Ra = step / that extrapolated current. Total resistance = step / steady-state current. Rm = total minus Ra. Cm = tau × (Ra+Rm) / (Ra × Rm). The ×1000 factors only convert units.
6. Two cross-checks: the same fit on the *off* transient (`tau_off_ms`, should match `tau_ms`), and Cm from the charge under the transient (`Cm_charge_pF`, should match `Cm_pF`).
7. Flags are added for high Ra, low Rm, large holding current, or a poor fit.

#### Finding action potentials: `detect_aps`

This is the engine for firing, rebound and AHP analysis.

1. Compute dV/dt as a central difference over 0.2 ms.
2. A candidate AP starts where dV/dt crosses `DVDT_THRESH` (20 mV/ms) upward.
3. The peak is the highest point within 3 ms. The candidate is kept only if the peak is above `AP_MIN_PEAK_MV` and at least `AP_MIN_AMP_MV` above the threshold.
4. After accepting an AP, the next 1.5 ms are skipped so one AP is never counted twice.
5. For each accepted AP it measures: threshold (voltage at the crossing), peak, amplitude, half-width (linear interpolation at half height), fAHP and AHP depth, time of the AHP trough, and the maximum rise and decay slopes.

Voltages inside this function are raw (recorded). The LJP correction is applied by the caller through `Vc()`.

#### Step protocols: `step_geometry`, `sweep_step_levels`, `base_ss`

These three small helpers are shared by HCN-CC, firing and AHP.

- `step_geometry` finds the step onset and offset using the sweep with the largest step, then reuses those times for all sweeps.
- `sweep_step_levels` returns the baseline current (median of the 100 ms before) and the step size (mean of the middle 80% of the step, minus baseline). Subtracting the baseline is what makes cells held with bias current analyse correctly.
- `base_ss` returns the baseline voltage and the steady-state voltage for a sweep.

#### HCN in current clamp: `analyze_hcn_cc`

For each of the six hyperpolarising steps it finds the **peak** (lowest smoothed voltage during the step), the **steady state** (last 100 ms), and from those the sag amplitude, ratio and percentage. It then looks in the 500 ms after the step for the highest voltage (rebound) and counts APs there. Finally it fits an exponential to the smallest hyperpolarising step (2 to 150 ms after onset) to estimate tau\_m.

#### HCN in voltage clamp: `analyze_hcn_vc`

Finds the step in the command trace, then for each sweep measures the current just after the capacitive spike has decayed (`I_inst`, 8 to 13 ms after onset), the current at the end of the step (`I_ss`, last 50 ms) and the current after the step ends (`I_tail`, 8 to 28 ms after offset). Ih amplitude = I\_ss minus I\_inst. The first sweep steps to -60 mV (no step) and is ignored in the summary.

#### Firing: `analyze_firing`

Per sweep: count APs inside the step, rate (count / step duration), latency to the first AP, ISIs, adaptation, a depolarisation block flag, and the full AP measures of the first AP.

Across sweeps:

- **Rheobase** = the first sweep with at least one AP.
- **Rin** = slope of steady-state voltage against current over sweeps with no AP and a steady state within 15 mV of baseline (needs at least 3 points).
- **Max rate** = highest rate among sweeps not flagged as blocked.
- **f-I slope** = linear fit over the first up to 6 spiking sweeps.
- **Adaptation reference** = the first sweep with at least 5 spikes.
- **RMP estimate** = median baseline voltage over all sweeps.

The block flag fires when a sweep has APs, ends depolarised above -25 mV (recorded), and the last AP came before 70% of the step.

#### AHP protocols: `analyze_ahp`

For each sweep it finds the current pulse, takes the pre-pulse baseline, detects the AP the pulse evoked, and measures troughs after the peak: fAHP (0 to 5 ms), mAHP (5 to 100 ms) and sAHP (0.1 to 1.5 s, using 50 ms smoothing). Values are *baseline minus trough*, so positive means a real dip below baseline. The pulse itself creates a voltage artefact, so the dV/dt threshold is meaningless here. AP amplitude and half-width are therefore measured from the baseline instead.

#### Gap-free: `analyze_gapfree`

First some QC: duration, median holding current (whole file, first 5 s, last 5 s), drift in pA per minute, and noise (SD after removing each 0.5 s chunk's mean). Then the event detection, in the same function: it calls `detect_psc` (which uses `slow_baseline`, `matched_filter_criterion` and `measure_event`) and summarises the events with `psc_stats`. In plain words, for every recording it does five things:

1. **Removes the slow baseline.** The holding current wanders by tens of pA over seconds. A running median over 0.4 s is subtracted, so what is left is a flat trace with the fast events on top.
2. **Slides a template over the trace** (matched filter, the method behind Clampfit's template search). A model event (quick rise, exponential decay) is fitted at every position; where it fits clearly better than noise, an event is called. It uses a fast template for inward events (EPSC-like) and a slower one for outward events (IPSC-like), plus extra templates: a short one so that an event right next to another is not missed, and a very slow one for IPSCs that rise over 10 to 30 ms.
3. **Measures every candidate** on a 2 kHz low-passed copy: baseline just before it, peak, amplitude, 10 to 90 % rise time, half-width and decay time constant.
4. **Throws out what is not a synaptic event:** candidates smaller than the detection limit, the same peak found twice (or two peaks with no clear dip between them), and anything above 250 pA (unclamped action potentials or artefacts, plus 0.3 s after such an event), and unstable stretches where the baseline swings far more than in the rest of the file (membrane or seal movement).
5. **Joins the holding.** `load_holdings` and `lookup_hold` read `holdings.csv` (one row per gap-free file, matched on date and run number) so each recording knows its true and commanded holding, the internal solution and how sure the holding is. `main()` then adds the nearest cap checks (Ra), builds the QC flags, and writes three tables: one row per recording, one row per event, and one row per cell and holding.

#### Plot functions

`plot_capcheck`, `plot_hcn_cc`, `plot_hcn_vc`, `plot_firing`, `plot_ahp` and `plot_gapfree` build the PNGs with matplotlib. If you want to change what a figure shows, they are self-contained and safe to edit.

### Bookkeeping and `main()`

This block links individual files into cells and writes everything out. It is the part you will touch most when you want a new column in the summary.

#### File names to cell keys: `NAME_RE` and `parse_name`

`NAME_RE` is a regular expression that matches your new naming convention, for example `20261007_C1_CapCheck_0005` or `20261012_F_DIV13_S1_C5_TTX_CapCheck_0034`. `parse_name` uses it to pull out date, sex, DIV, slice, cell, condition, tags and run number.

The **cell key** is the file name with the protocol and run number cut off (`20261007_C1`). All files sharing a key are treated as one cell. If a file name does not match the pattern, it still runs, but its key starts with `UNRENAMED_`. If you ever change the naming convention, `NAME_RE` is where to update it.

The rig (Sacha or Gerda) is taken from the first word of the top-level folder. If the same key appeared on two rigs, the rig name is added to the key so the cells stay separate.

#### The main loop

Inside `main()` the per-file loop does, in order:

1. **Duplicate check.** Same recording time and same file size means the second copy is logged and ignored.
2. **Re-routing.** If the header says HCN or firing but the data are a short +500 pA pulse, the file is treated as an AHP recording and logged.
3. **Dispatch.** An `if / elif` chain calls the right `analyze_*`, appends rows to the right table and, unless `--no-plots`, calls the plot.
4. **Error handling.** A `Skip` is logged as SKIPPED. Any other error is logged as ANALYSIS FAILED with a short traceback, and the loop carries on. One bad file never stops the run.

#### The tables

Every recording writes rows into dictionaries stored in `T` (for example `T["capcheck"]`, `T["firing_sweeps"]`). At the end `write_csv` saves them to `_analysis\tables\`. These are the detailed, unfiltered numbers.

#### The per-cell summary

After the loop, files are grouped by cell key and sorted by recording time. For each cell the code builds one row (`row[...] = ...`). Rules worth knowing:

- **Cap check:** `first` uses the earliest cap check, `last` the latest. Ra drift = (last - first) / first, flagged if above 20% over more than a minute.
- **HCN, firing:** the *first* file of that protocol is used (`by[kind][0]`). Repeats stay in the tables.
- **Firing columns** are prefixed `f25_` or `f10_` for the two protocols.
- **AHP:** the first AHP file that contains at least one AP.
- **RMP estimate:** median of the baseline voltages from the firing and HCN recordings.
- **`qc_flags`:** all warnings for that cell joined into one text cell.

#### Metadata

`metadata_template.csv` is written once, listing every cell key with empty columns (`rmp_manual_mV`, `hold_exc_mV`, `hold_inh_mV`, `notes`, plus condition, DIV and sex). It is never overwritten. Fill it in, then run with `--metadata metadata_template.csv`. Rows are matched by `cell_key`. A manual RMP replaces the estimate in `rmp_mV` and `rmp_source` says which was used. Metadata rows that match no cell are reported in `qc_log.txt`.

#### Outputs

`summary.csv`, the tables, the PNGs in `quicklook\`, and `qc_log.txt` all go to `_analysis` (or to `--out`).

### Cheat sheet: common tweaks and where to make them

| I want to... | Search for | Change |
| --- | --- | --- |
| Use a different LJP, or none | `--ljp` on the command line, or `LJP_MV` | A number. `--ljp 0` gives raw Clampfit-style voltages. |
| Use a stricter or looser AP threshold | `DVDT_THRESH` | 10 to 50 mV/ms are all used in the literature. Report which one you used. |
| Change the steady-state window | `SS_WINDOW_S` | Seconds. Affects sag, Rin and steady-state voltage. |
| Change QC cut-offs | `RA_MAX_MOHM`, `RM_MIN_MOHM`, `RA_DRIFT_MAX_PCT` | Numbers in `CFG`. |
| Report sag at one fixed step only | `hcn_sag_amp_mV_` in `main()` | Replace the loop over steps with one named step. |
| Add a new column to the summary | the block that starts `# HCN CC` or `# firing` in `main()` | Add `row["my_column"] = value`. Columns appear in the CSV automatically. |
| Add a new measurement per sweep | the `row = OrderedDict(...)` in the matching `analyze_*` | Add a key. It lands in the table CSV automatically. |
| Use more than the first firing run | `by[kind][0]` | Pick `[-1]`, or loop over all of them. |
| Use a different cap check as "first" | `caps[0]`, `caps[-1]` | Index into the list (it is sorted by time). |
| Change the Rin range | `RIN_MAX_DEPOL_MV` | mV above baseline allowed for the fit. |
| Change the adaptation reference | `ref = [x for x in rows if x["n_spikes"] >= 5 ...` in `analyze_firing` | The spike count, or pick a current instead. |
| Change the depolarisation block rule | `block_flag` in `analyze_firing` | The -25 mV level or the 70% point. |
| Change AHP windows | the three `trough(...)` calls in `analyze_ahp` | Start and end times in seconds. |
| Change the rebound window | `REBOUND_WINDOW_S` | Seconds after the step. |
| Restyle figures | `C1`, `C2`, `C3`, `plt.rcParams`, or inside a `plot_*` | Colours, sizes, labels. |
| Add a new protocol | `classify`, then a new `elif r.kind == ...` in `main()` | Write an `analyze_` function that returns numbers. |
| Change the file naming pattern | `NAME_RE` | Keep the named groups (`date`, `cell`, ...). |
| Re-run only part of the data | `--only text` | Any text that is in the file path. |
| Skip the plots to run faster | `--no-plots` | A flag, no edit needed. |

After any edit, run on one day with `--no-plots`. If a column goes empty or a count looks wrong, `qc_log.txt` says which file failed and why.

## Part 2: What the readouts mean

The ranges below are rough guides for whole-cell recordings from pyramidal neurons in slices, to give you a sense of scale. Your own practice cells (for example Cm of 90 to 160 pF, Ra of 15 to 28 MΩ) already fall inside them. Treat them as orientation, not as pass/fail rules.

### Voltage clamp: the cap check

The cap check clamps the cell at -60 mV, jumps 5 mV for 20 ms and watches the current. The cell behaves like a simple circuit: a resistor (the pipette opening, Ra) in series with a parallel resistor and capacitor (the membrane, Rm and Cm). Every readout below comes from that picture.

**Ra (access or series resistance, MΩ).** How hard it is for current to get from the amplifier into the cell through the opening you made in the membrane. Typically 10 to 25 MΩ. It matters most of all quality measures:

- It makes the voltage the cell really sees differ from the command (error = current × Ra).
- It slows and filters fast currents, which shrinks and widens synaptic events, so sEPSC and sIPSC amplitudes depend on it.
- It tends to creep upwards as the opening partly reseals. High Ra means unreliable synaptic data from that cell.

*How Ra is derived.* The script averages the 10 sweeps to cut noise. At the very instant of the 5 mV step the membrane acts like a short circuit, so only Ra limits the current. The script fits a decaying curve to the current spike, runs the curve back to the exact moment of the step to read the current at time zero, and divides the step by that current. For cell 1 on 7 Oct: 5.06 mV ÷ 238 pA = 21.3 MΩ.

**Rm (membrane resistance, MΩ).** How tightly the membrane and seal hold the cell electrically together. A healthy pyramidal cell is often 100 to 500 MΩ. Low Rm means a leaky cell (damage, poor seal) or a very large one. It is the voltage-clamp counterpart of Rin in current clamp.

*How Rm is derived.* Once the spike has faded, the current settles to a flat level. The script averages the last 5 ms of the step to get that steady current, and divides the step by it. That gives the total resistance, Ra plus Rm, so Rm is the total minus Ra. For cell 1: 5.06 mV ÷ 36.7 pA = 138.1 MΩ, minus 21.3 MΩ Ra, gives 116.8 MΩ.

**Cm (membrane capacitance, pF).** The size of the cell's membrane. About 1 µF per cm² of membrane, so 100 pF is roughly 10,000 µm² of surface. Bigger cells with more dendrites have higher Cm. Cm is also the first quick sign of structural differences between groups, which is relevant to a question about pruning. Note that it only measures membrane that is electrically close to the pipette.

*How Cm is derived.* The decay speed of the spike depends on how long the capacitor takes to discharge through the resistance around it, which is Ra and Rm side by side. So Cm = tau ÷ (Ra and Rm in parallel). For cell 1: Ra and Rm in parallel are 18.0 MΩ, so 2.34 ms ÷ 18.0 MΩ = 130 pF. The script cross-checks this with a second method: the total charge under the spike (its area) divided by the step, with a small correction. That value is `Cm_charge_pF`.

**tau (ms).** How fast the cell charges after the voltage step. It is roughly Ra × Cm. A small Ra gives a fast clamp. If `tau_ms` and `tau_off_ms` disagree by more than 30%, the cell is not behaving like a single compartment and the numbers are less reliable.

*How tau is derived.* The fit finds the time the spike takes to fall to about 37% of its starting height above the steady level. The fit runs from the top of the spike to the end of the step. The script repeats the fit on the spike that follows when the step ends (the *off* transient). That gives `tau_off_ms`, which should be close to `tau_ms`.

**Ihold (holding current, pA).** The steady current the amplifier injects to keep the cell at -60 mV command. Near 0 to -100 pA is comfortable. A strongly negative value means the cell is depolarised or leaky and the amplifier is fighting to hold it. In the summary this is `cap_first_Ihold_pA`. Above 200 pA in magnitude it is flagged.

*How Ihold is derived.* It is the average current in the averaged trace from 12 ms to 1 ms before the step, when nothing is changing. Every other current in the cap check is measured relative to this.

**hold\_cmd\_mV and hold\_cmd\_LJPcorr\_mV.** The commanded holding voltage, before and after subtracting the LJP. With K-gluconate and a 15 mV LJP, a -60 mV command is closer to -75 mV at the membrane.

*How they are derived.* `hold_cmd_mV` is the average command voltage in the same pre-step window. The LJP-corrected version is that number minus 15 mV (or whatever `--ljp` you set).

**Ra drift (%).** Ra at the last cap check compared to the first. More than about 20% means the later recordings in that cell were made under worse conditions.

*How Ra drift is derived.* The cell's cap checks are put in time order. Drift = (Ra of the last − Ra of the first) ÷ Ra of the first × 100. It needs at least two cap checks, and it is only flagged if they are more than a minute apart.

**Technical columns.** `Iss_pA` (steady-state current in the step), `I_peak_pA` (extrapolated peak current), `fit_rms_pA` (misfit of the exponential, small is good), `Cm_charge_pF` (an independent Cm estimate from the charge, should be within about 10% of `Cm_pF`).

*How they are derived.* `Iss_pA` is the average of the last 5 ms of the step minus Ihold. `I_peak_pA` is the fitted spike extrapolated back to the moment of the step. `fit_rms_pA` is the typical gap between the fitted curve and the real trace. `Cm_charge_pF` is the area under the spike divided by the step, with a correction for the current that leaks through Rm.

### Voltage clamp: HCN (Ih)

This protocol steps the cell from -60 mV to -70, -80, -90 and -100 mV for 500 ms.

**Ih (the HCN current).** HCN channels open slowly when the cell is hyperpolarised and let a mixed Na and K inward current through. In the trace you see an instantaneous jump (leak) followed by a slow extra inward drift. It is the signal behind sag in current clamp.

- **`I_inst_pA`** is the current just after the capacitive spike, before HCN channels have opened.
- **`I_ss_pA`** is the current at the end of the step.
- **`Ih_amp_pA` = `I_ss` minus `I_inst`.** More negative means more Ih. It grows with the size of the step.
- **`tail_amp_pA`** is the current just after returning to -60 mV minus the pre-step current. The channels close slowly, so the tail reflects how much Ih had opened. It confirms the Ih measure from a second angle.

*How the Ih numbers are derived.* For each step the script averages the current in three short windows. `I_inst` uses 8 to 13 ms after the step starts: long enough for the capacitive spike to fade, short enough that HCN has barely opened. `I_ss` uses the last 50 ms of the step. Ih = `I_ss` minus `I_inst`. `I_tail` uses 8 to 28 ms after the step ends. `tail_amp_pA` is that minus the current averaged over the 50 ms before the step began. Every current is in pA and negative means inward.

What this does not do: it does not isolate Ih from other slow currents. For a pharmacologically clean value you would compare with and without an HCN blocker such as ZD7288.

### Current clamp: resting properties

Here you inject known currents and read the voltage. Ohm's law (V = I × R) turns the response into the cell's properties.

**Baseline Vm and RMP estimate (mV).** The voltage with 0 pA injected, taken before each step. It is set by a balance of leak currents, mostly K+, and is a health indicator: roughly -65 to -75 mV (LJP-corrected) is normal for a pyramidal cell, and a value like -45 mV suggests a damaged or depolarised cell. The summary reports it as `rmp_est_mV`. Your manually recorded RMP overrides it as `rmp_mV`.

*How it is derived.* Before each current step, the script averages the voltage over the 100 ms leading up to it (stopping 2 ms before the step). The RMP estimate is the median of these baselines across the cell's firing runs and its first HCN run, minus the 15 mV LJP. Using the median stops one odd sweep from moving it.

**Bias current (`bias_pA`).** The current that was already being injected before the step. In a good cell it is about 0 (the few pA you see are monitor offset). If it is large (for example -290 pA), the cell was being held at a more negative voltage by the amplifier, and everything measured in that run is for a held cell, not a resting one. This is flagged in `qc_flags`.

*How it is derived.* The current-monitor channel is averaged over the 100 ms before each step, and the median across sweeps is reported. That is what the amplifier was injecting before the step started. The size of each step is measured relative to it, so held cells are still analysed correctly.

**Rin (input resistance, MΩ).** How much the voltage changes per unit of injected current: the slope of steady-state voltage against current. A high Rin means the cell is easy to excite and often small or with few open channels at rest. Rin from the firing protocol uses small depolarising steps. `hcn_Rin_m100_MOhm` uses the -100 pA step instead. They differ because membrane conductance changes with voltage (rectification).

*How Rin is derived.* From every sweep that has no AP and stays within 15 mV of baseline, the script takes two numbers: the step size in pA and the rise in steady-state voltage in mV (the 0 pA sweep counts as the origin). It draws the best straight line through those points. The slope, converted from mV per pA to MΩ, is Rin. For example, a 25 pA step that moves the voltage by 2.8 mV gives 2.8 ÷ 25 = 0.112 GΩ, which is 112 MΩ. It needs at least 3 points. `hcn_Rin_m100_MOhm` is simpler: the steady-state change at the -100 pA step divided by -100 pA.

**tau\_m (membrane time constant, ms).** How quickly the voltage follows a current step, roughly Rin × Cm. It sets how long inputs are integrated. Here it is estimated from the smallest hyperpolarising step, so treat it as approximate: sag and dendrites bend the curve.

*How tau\_m is derived.* The script takes the smallest hyperpolarising step (-100 pA) and looks at the voltage from 2 ms to 150 ms after the step starts. It fits a curve that eases from the starting voltage toward a steady level, and tau\_m is the time constant of that curve: the time to cover about 63% of the way. Because HCN and dendrites bend the curve, it is only an estimate.

### Current clamp: sag and rebound (HCN)

When the cell is pushed to a very negative voltage, HCN channels open slowly and let inward current in, so the voltage drifts back up toward rest. That drift is the sag.

- **`peak_mV`.** The lowest voltage reached early in the step, before HCN opens.
- **`ss_mV`.** The voltage at the end of the step, after HCN opened.
- **`sag_amp_mV`.** `ss_mV` minus `peak_mV`. Larger means more Ih. About 3 to 12 mV is typical for large steps.
- **`sag_ratio`.** (peak change) / (steady-state change). 1 means no sag. 1.3 means the peak was 30% deeper than the end of the step.
- **`sag_pct`.** The same thing as a percentage: how much of the early deflection was recovered.
- **`rebound_amp_mV`.** After the step ends, the voltage overshoots baseline before settling. It happens because the opened HCN channels take time to close. Bigger sag usually means bigger rebound.
- **`rebound_spikes`.** The number of APs fired in that rebound. Usually 0 unless the rebound reaches threshold.

*How the sag and rebound numbers are derived.* For each step: baseline is the voltage averaged over the 100 ms before; peak is the lowest point during the step after a 1 ms running average, which stops one noisy sample from faking a minimum; steady state is the average of the last 100 ms. Then `sag_amp_mV` = steady state minus peak. The two deflections are peak minus baseline and steady state minus baseline, and `sag_ratio` is the first divided by the second. `sag_pct` = 100 × (1 − 1 ÷ ratio). For cell 1 at -200 pA: baseline -60.8, peak -85.1, steady state -78.7 mV. The deflections are -24.3 and -17.9 mV, so sag is 6.4 mV and the ratio is 1.36. For rebound, the script takes the 500 ms after the step ends, smooths 1 ms, takes the highest point and subtracts baseline. `rebound_spikes` counts APs in that window with the same AP detector used for firing.

Sag size depends on the step, so the summary has one column per step (`_m100` to `_m600`). HCN makes the cell less sensitive to slow input and affects dendritic integration, so it is a common readout when comparing groups of cells.

### Current clamp: single action potential

These values come from the first AP at rheobase and are in the summary with the prefix `f25_` or `f10_`.

**Rheobase (pA).** The smallest current step that fires at least one AP. It combines how far the cell is from threshold and how much voltage a pA produces (Rin). A lower rheobase means a more excitable cell. Resolution is 25 pA in the 25 pA protocol and 10 pA in the 10 pA protocol.

*How it is derived.* The script goes through the sweeps from the smallest current upward and stops at the first one where the AP detector finds at least one AP during the step. That sweep's step size is the rheobase. The step is the measured current rounded to the protocol's step (25 or 10 pA), so a measured 25.2 pA is reported as 25.

**AP threshold (mV).** The voltage at which dV/dt first reaches 20 mV/ms. It is where voltage-gated Na+ channels take over. More negative means easier to fire. The absolute value depends on the criterion, so always say which one you used.

*How threshold is derived.* The script computes how fast the voltage is changing (dV/dt, in mV per ms) over a 0.2 ms window. The first time that speed climbs past 20 mV/ms is taken as the start of an AP. The threshold is the voltage at that moment, minus the LJP. This is a speed criterion, not a visual one, so it is the same for every spike and every cell.

**AP peak and amplitude (mV).** Peak is the top of the spike (LJP-corrected). Amplitude is peak minus threshold. Both depend on how many Na+ channels are available. They get smaller in unhealthy cells and during long trains.

*How peak and amplitude are derived.* The peak is the highest voltage within 3 ms after the threshold point. Amplitude is the peak minus the threshold voltage. A candidate only counts as an AP if the peak is above -10 mV (recorded) and at least 25 mV above its threshold.

**Half-width (ms).** The width of the AP at half its amplitude. It reflects how fast K+ channels repolarise the cell. Longer means slower repolarisation. It depends on temperature.

*How half-width is derived.* The script works out the halfway voltage, which is the threshold plus half the amplitude. It finds the time the voltage crosses that level on the way up and again on the way down, using straight-line interpolation between samples so the answer is finer than the 0.05 ms sampling. Half-width is the time between the two crossings.

**Max rise and max decay (mV/ms).** The steepest upstroke and downstroke. Rise reflects Na+ channel density and availability. Decay reflects K+ repolarisation.

*How they are derived.* Max rise is the largest dV/dt between threshold and peak. Max decay is the most negative dV/dt in the 4 ms after the peak.

**fAHP and AHP (mV).** After the spike, the voltage often dips below threshold. `ap_fahp_mV` is the depth within 5 ms, and `ap_ahp_mV` within 50 ms, both measured down from threshold. Positive means a dip. Negative means no dip. Fast AHPs shape how soon the next spike can fire.

*How they are derived.* After the peak, the script looks for the lowest voltage in a short window. It smooths with a 0.5 ms average, so one noisy sample cannot fake a dip. The fAHP window is 5 ms and the AHP window is 50 ms, ending sooner if the next AP starts. The value is the threshold minus that lowest voltage, so a dip below threshold is positive.

**Latency (ms).** Time from step onset to the first spike near rheobase. A long latency (hundreds of ms) points to slowly inactivating K+ currents that delay firing.

*How it is derived.* The time from the start of the current step to the peak of the first AP in that sweep.

### Current clamp: repetitive firing

**f-I curve.** Firing rate (spikes per second) against injected current. It is the cell's input-output function. The figure shows it for each recording. The tables hold the rate for every step.

*How it is derived.* In each sweep the script counts the APs inside the step window and divides by the step length (0.5 s) to get spikes per second. Plotting that against the step current gives the curve.

**f-I slope (Hz/pA).** The gain: how many extra spikes per second you get for each extra pA, fitted over the first up to six spiking steps. Steeper means a more sensitive cell.

*How it is derived.* The script takes the first up to six sweeps that contain spikes and draws a straight line through (current, rate). The slope of that line is the gain in Hz per pA. It needs at least 3 spiking sweeps, otherwise it is left blank.

**Max rate (Hz) and the current where it occurs.** The highest rate reached without depolarisation block.

*How it is derived.* The highest rate among all sweeps that were not flagged as blocked, and the current of that sweep.

**Depolarisation block.** At high currents a cell can stop firing because Na+ channels inactivate while the voltage sits depolarised. Sweeps like this are flagged (`depol_block_seen`) and left out of the max rate. A block at low current in a cell with a depolarised RMP is usually a sign of poor health.

*How it is derived.* A sweep is flagged when it contains at least one AP, the voltage averaged over the last 50 ms of the step is above -25 mV (recorded), and the last AP happened before 70% of the step had passed. In words: the cell fired, then stopped early while sitting depolarised.

**ISI, adaptation ratio and adaptation index.** ISI is the time between spikes. Cells slow down during a step because slow K+ currents build up (spike-frequency adaptation).

- **`adapt_ratio`.** Last ISI divided by first ISI. 1 means no adaptation, 2 means the last interval is twice the first.
- **`adapt_index`.** (last - first) / (last + first). 0 is no adaptation. A value of 0.27 means the last ISI is about 1.7 times the first. Values near 1 mean strong adaptation.
- **`isi_cv`.** Spread of ISIs, a measure of regularity.

*How the adaptation numbers are derived.* ISIs are the gaps between the peak times of neighbouring APs in the step. `adapt_ratio` = last ISI ÷ first ISI. `adapt_index` = (last − first) ÷ (last + first). `isi_cv` = standard deviation of the ISIs ÷ their mean. All three need at least 3 ISIs (4 spikes). In the summary, the script uses the first sweep that has at least 5 spikes.

The summary reports adaptation at the first step that has at least 5 spikes (`adapt_ref_pA`), so cells are compared at a similar stage of firing.

### AHP protocols

These use a brief (2 to 7 ms) strong current pulse to make exactly one AP, then watch the recovery.

- **`ap_amp_from_base_mV` and `ap_halfwidth_ms`.** Measured from the pre-pulse baseline instead of a threshold, because the pulse itself distorts the early voltage.
- **fAHP (0 to 5 ms), mAHP (5 to 100 ms), sAHP (0.1 to 1.5 s).** Depth of the dip below baseline in three time windows: fast, medium and slow after-hyperpolarisation. They come from different channel types (broadly BK and Kv for fast, SK, M and HCN for medium, Ca2+-dependent K+ for slow).
- **Sign.** Positive means the voltage dips below baseline. Negative means it stays above baseline, which is an after-depolarisation, not an AHP. The Gerda cell \`20260918\_C1\` shows this.
- **`_t_ms` columns.** When the trough occurs after the peak. If it equals the window edge, the voltage was still moving and there is no true trough in that window.

*How the AHP protocol numbers are derived.* The baseline is the voltage averaged over the 100 ms before the pulse. The AP detector looks for one AP between pulse start and 15 ms after pulse end. Amplitude from baseline = AP peak minus baseline. Half-width is measured at the level halfway between baseline and peak. The three troughs are the lowest voltage in each window after the peak (5 ms, 5 to 100 ms, and 0.1 to 1.5 s). The first two use a 1 ms running average, and the slow one uses a 50 ms average because slow AHPs are small and noisy. Each value is baseline minus trough. The `_t_ms` columns are the time of the trough after the peak.

### Quality-control columns and gap-free numbers

**`qc_flags`.** One text cell per cell in `summary.csv` listing everything that looks off: high Ra, low Rm, large holding current, Ra drift, a held-with-bias-current cell, or skipped files. An empty cell means nothing was flagged. A flag is a reason to look at the quick-look PNG, not an automatic exclusion.

*How it is derived.* While building each cell's row, the script checks a list of rules (Ra above 30 MΩ, Rm below 100 MΩ, holding current above 200 pA, Ra drift above 20%, bias current above 50 pA, skipped files). Each rule that triggers adds a short phrase. The cut-offs are in `CFG`.

**`n_files_skipped`.** How many of that cell's files were not analysed because of wrong channel units or implausible traces. The file names and reasons are in `qc_log.txt`.

**`n_capchecks`, `Ra_drift_span_min`.** How many cap checks the cell has, and the minutes between its first and last one. Drift over a short span is less worrying than the same drift over an hour.

**Gap-free QC (`gapfree_qc.csv`).**

- **`duration_s`.** Nominally 120 s. Runs under 60 s are logged as probably aborted.
- **`Ihold_median_pA`, `Ihold_start_pA`, `Ihold_end_pA`.** The holding current through the run. A big difference between start and end means the cell was changing.
- **`drift_pA_per_min`.** Slope of the holding current. Large slow drift means an unstable recording.
- **`noise_sd_pA`.** Baseline noise. This is the number that decides which small synaptic events you will be able to see: events much smaller than about 3 to 4 times this value will be lost. Gerda is noisier, so expect higher values there.

*How the gap-free numbers are derived.* Duration is the number of samples divided by the sampling rate. The three holding-current values are medians of the whole file, the first 5 s and the last 5 s. Drift is a straight-line fit through the current averaged in 10 ms blocks, expressed per minute. Noise is found by cutting the trace into 0.5 s pieces, subtracting each piece's own average, and taking the standard deviation of what is left. That removes slow drift but keeps fast noise, and it also includes any synaptic events, so it overestimates pure noise in a busy cell.

### Voltage clamp, gap-free: sEPSC and sIPSC

These are the spontaneous synaptic currents in the 2-minute recordings (`psc_recordings.csv`, one row per file; `psc_by_cell_holding.csv`, files of one cell at the same holding pooled; `psc_events.csv`, every single event). The holding matters because it decides what you can see. With the K-gluconate internal, chloride sits near -90 mV. At a true holding of -75 to -35 mV that is below the holding, so GABA-A (inhibitory) currents flow **outward** (upward), while glutamatergic (excitatory) currents flow **inward** (downward). The script therefore uses direction as the label: inward events are called sEPSCs and outward events sIPSCs. Close to -90 mV inhibitory events shrink to nothing (no driving force); close to 0 mV excitatory events disappear. With the high-chloride (KCl) internal that logic breaks: chloride currents are inward too, so inward events are a mix of excitatory and inhibitory, and the flag says so.

The holding comes from `holdings.csv`, not from the .abf file (the header does not store it). `hold_true_mV` is the real membrane voltage after the 15 mV junction-potential correction; `hold_cmd_mV` is what was typed into Clampex. Until 6 Oct they differ by 15 mV (a typed -60 was really -75); from 6 Oct, and on 25 Sep, the typed value was already corrected.

**Inward and outward events (`inward_n`, `outward_n`).** How many events of each direction were accepted. Inward is excitatory (EPSC-like), outward is inhibitory (IPSC-like) for the K-gluconate cells.

*How they are derived.* Each direction has its own template: inward rises with 1.2 ms and decays with 5 ms, outward rises with 2 ms and decays with 9 ms. At every sample the template (scaled and offset) is fitted to the following 30 ms (inward) or 55 ms (outward) and the script computes scale divided by the standard error of the fit. Where that number exceeds 4, there is an event. The extra templates work the same way: inward 1.2 ms rise and 4 ms decay in a 12 ms window (a neighbouring event does not spoil the fit); outward 3 ms and 8 ms in 20 ms, and 6 ms and 25 ms in 90 ms (slow IPSCs). A detection from an extra template is kept only if no earlier one lies within 4 ms. Two detections closer than 3 ms are one, and two peaks of the same direction closer than 12 ms (inward) or 30 ms (outward) count as separate events only if the trace dips between them by at least the smallest accepted amplitude. Outward events can rise slowly, so their baseline is taken just before the foot of the rise (the lowest point in the 40 ms before the peak). Then each candidate is measured (below) and kept only if its amplitude is at least 3 pA and at least 4.5 times the noise (inward) or 5.5 times (outward). Those two limits come from a test on phase-randomised copies of your own recordings (same noise spectrum, no real events): between 3.5 and 4.5 times the noise about two thirds of the inward detections and all of the outward ones were noise; between 4.5 and 5.5 times about 7% of the inward ones were noise, and outward needed 5.5.

*How well it works.* 321 of your 341 Clampfit screenshots (25 Sep, 7 Oct, 8 Oct) could be matched automatically to their recording: the run number comes from the Clampfit title, the axes from the gridlines, and the trace is matched to the raw data. In them I counted 378 event-shaped excursions and the script finds 331 (88%): 274 of 308 inward and 57 of 70 outward. The misses are mostly events under about 12 pA, close to the noise limit, plus about 9 giant spikes (above 250 pA) that are excluded on purpose. Synthetic events added to your real traces are recovered about 85 to 100% of the time above 10 pA, with amplitudes about 5 to 10% too large. On noise-only copies of your traces the false detections are about 7% of inward events between 4.5 and 5.5 times the noise and about 5% of outward events between 5.5 and 7 times the noise. Events close to the limit should be treated with caution, and the 10 pA columns are the safer numbers.

**Frequency (`inward_freq_Hz`, `outward_freq_Hz`).** Events per second. The most direct readout of how many active synapses the cell receives, and how often they fire.

*How it is derived.* Accepted events divided by the analysable time: the file duration minus any stretch masked after a very large event (`masked_s`). Frequency depends on the detection limit, so it is only comparable between recordings with similar noise; use `inward_freq_ge10pA_Hz` (next but one) when they differ.

**Amplitude (`amp_mean_pA`, `amp_median_pA`, `amp_sd_pA`).** How big a typical event is: how many receptors open at one synapse, and how well the current gets from the synapse to the pipette. The median is safer than the mean because a few large events pull the mean up.

*How it is derived.* For each event, the baseline is the median of a 6 ms window that ends 5 ms before the peak (slow-rising events start before the detector marks them). The peak is the highest point (lowest for inward) of the 2 kHz low-passed trace within 10 ms (inward) or 15 ms (outward) of the detection. Amplitude = peak minus baseline, always reported as a positive number. The numbers are summarised over all accepted events. Tests with fake events put into your real traces showed the amplitude comes out about 5 to 15 % high, because picking the highest point of a noisy trace adds a little noise.

**Rise time (`rise_10_90_ms_median`).** How fast the current builds up. Fast rise means the synapse is close to the soma (little filtering by the dendrite); slow rise means a distant synapse, or a poor access resistance that filters everything.

*How it is derived.* Walking back from the peak on the low-passed trace, the script finds the last time the current was at 10 % of the amplitude and at 90 %, interpolating between samples. Rise = the time between them. Shown as a median, and only when the recording has at least 10 events.

**Decay (`decay_tau_ms_median`).** How fast the current dies away. For EPSCs this follows the AMPA receptor and glutamate clearance; for IPSCs it follows GABA-A receptor kinetics (often slower, around 8 to 15 ms).

*How it is derived.* From the peak, a single exponential (amplitude times exp(-t/tau)) is fitted over 5 time-constants or until the next event, whichever is shorter. The fit must explain at least half the variance, otherwise that event gets no tau. The median of the good fits is reported.

**Half-width (`halfwidth_ms_median`).** The duration of the event at half its height. A simple width that mixes rise and decay and needs no fit.

*How it is derived.* The time between the rising and the falling crossing of half the amplitude, on the low-passed trace.

**Events of 10 pA or more (`n_ge10pA`, `freq_ge10pA_Hz`).** The same count and frequency, but only for events at least 10 pA. A noisy recording cannot see small events, so its total frequency is artificially low; counting only events everybody can see gives a fairer comparison between the two rigs. The 10 pA limit is `COMMON_FLOOR_PA` in the settings; recordings noisier than about 3 pA cannot reach it.

**Noise and the smallest detectable event (`noise_hp_lp2k_pA`, `min_detectable_pA`).** How quiet the baseline is, and the smallest event the script can accept in this recording.

*How they are derived.* Noise = 1.4826 times the median absolute deviation of the baseline-subtracted, 2 kHz low-passed trace (robust: events hardly change it). Smallest detectable = the larger of 3 pA and 4.5 times that noise for inward events (min\_detectable\_pA), and the larger of 3 pA and 5.5 times that noise for outward events (min\_detectable\_outward\_pA). A recording with noise of 4 pA cannot detect inward events under about 18 pA.

**Analysed seconds, masked seconds, excluded large events (`analysed_s`, `masked_s`, `n_excluded_large`).** How much of the file went into the numbers. Two things are masked. A very large inward transient (above 250 pA) is almost always a cell that fired an unclamped action potential, not a synaptic event: 0.3 s after it is not analysed. And any 0.5 s stretch whose baseline wobble (the median absolute deviation) is more than twice the file's typical value counts as unstable, and it plus its neighbouring stretches are not analysed, because slow waves of 20 pA would otherwise be mistaken for events. The large transients are counted in `n_excluded_large` and left out of all statistics.

**Ra before and after (`Ra_before_MOhm`, `Ra_after_MOhm`).** The access resistance from the cap checks just before and just after this recording. High Ra (above 30 MOhm) shrinks and slows events and makes the voltage clamp poor, so amplitude and kinetics from such files should be treated with care.

**Holding confidence (`hold_confidence`).** `high` = the journal lists exactly as many holdings as there are files; `medium` = the number of files was off by one or two and the script dropped the shortest, or the holding is the default (-60 mV typed, so -75 true) for days without notes; `low` = short or aborted file on a day without notes. A blank holding means the order could not be established; those rows are analysed but cannot be assigned to a holding.

**Flags (`psc_flags`).** One text per recording. It lists short duration, Ra above 30 MOhm, holding current above 200 pA, drift above 25 pA/min, noise above 4 pA, more than 10 % of the file masked, few events (kinetics not reported), events that cluster at the detection limit (many are probably noise), a holding that is unsure or unknown, and KCl internal. Nothing is removed because of a flag; you decide.

**How much to trust it.** The detector was tested by putting fake events of known size into real recordings (including the noisy Gerda ones). Above 10 pA it found 89 to 100 % of the events at noise levels of 1.6 to 2.4 pA, between 6 and 10 pA it found about half, and it added no extra detections beyond what the untouched recording already had. It has not yet been compared with events marked by hand or with Clampfit; that comparison is the next step.

### How the readouts relate to your thesis question

You are testing whether sex hormone exposure changes microglial pruning and excitatory/inhibitory balance. The readouts fall into three groups, and the first two are the context that makes the third interpretable.

- **Is the cell healthy and comparable?** Ra, Ra drift, Rm, Ihold, RMP, bias current. If groups differ here, differences elsewhere may be technical.
- **Is the cell itself different?** Cm (size), Rin, tau\_m, sag (HCN), AP features, rheobase, f-I. These describe intrinsic excitability. They can change without any change in synapses, and they change how synaptic input is converted into spikes.
- **Is the synaptic input different?** sEPSC and sIPSC frequency and amplitude (not yet implemented). Frequency changes usually point to a change in the number or release probability of synapses, which is the pruning question. Amplitude changes usually point to postsynaptic receptor changes.

Ra and noise limit what you can say in the third group. A recording with Ra of 40 MΩ or noise of 8 pA is a weak candidate for event detection, even if its intrinsic properties look fine.

Ask me for any value in the summary and I can show you where it comes from in the traces.
