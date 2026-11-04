# CryoFlowTools

A small collection of command-line utilities for cryo-ET and RELION
workflows on HPC clusters.

## Tools

| Tool | Location | Description |
|------|----------|-------------|
| `gpuview` | `hpc/` | Summarises Slurm GPU-node resources and running jobs |
| `jobdiag` | `hpc/` | Read-only Slurm job diagnostics for pending reasons, resource fit, and scheduler state |
| `missali_report` | `missalignment/` | Summarises MissAlignment `*_alignment_loss.json` results |
| `class_distribution` | `relion/` | Tracks RELION Class3D class distributions across iterations |
| `compare_class3d_jobs` | `relion/` | Compares particle class assignments between two RELION Class3D `*_data.star` files |
| `split_star_by_tomo` | `relion/` | Splits a RELION STAR file into one file per `rlnTomoName` or `rlnMicrographName` |
| `merge_star_files` | `relion/` | Merges compatible RELION STAR files while preserving particle metadata and other blocks |
| `m_progressive` | `m/` | Prepares RELION-to-M inputs and generates a progressive refinement shell script |
| `m_resolution` | `m/` | Summarises Warp/M species refinement version ancestry and resolution history |

## Installation

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

The scripts are standalone executables. Make them executable if needed:

```bash
chmod +x hpc/gpuview hpc/jobdiag relion/class_distribution relion/compare_class3d_jobs relion/split_star_by_tomo relion/merge_star_files missalignment/missali_report m/m_resolution m/m_progressive
```

`m_resolution` and `m_progressive` use only the Python standard library.
For `m_progressive`, use Python 3.9+ and Bash; load RELION's `relion_mask_create`
if generating a mask. Running the generated script also requires a working
Warp/M environment with `MTools`, `MCore`, and CUDA dependencies.

## Usage

### gpuview

Display Slurm GPU partition status (run on a login node):

```bash
hpc/gpuview
hpc/gpuview gpu_strubi
```

### jobdiag

Diagnose Slurm jobs using scheduler-visible state, requested resources, candidate-node capacity, and pending reasons.

```bash
hpc/jobdiag 34378338
hpc/jobdiag 34378338 --verbose
hpc/jobdiag 34378338 --no-color
hpc/jobdiag 34378338 --partition gpu_interactive
```

`jobdiag` handles both running and pending jobs, including resource shortages, node states, priority waits, scheduled starts, job holds, dependencies, partition configuration, and multi-partition submissions.

`--partition PARTITION` evaluates the existing resource request against an alternative partition without modifying or resubmitting the job. Use `--verbose` to show detailed evidence, classifications, and commands executed.

`jobdiag` is strictly read-only: it inspects Slurm state but never submits, cancels, holds, releases, requeues, or modifies jobs.

> Scheduler-visible evidence can be incomplete. When the exact blocker cannot be established reliably, `jobdiag` reports it as unknown rather than guessing.

### class_distribution

Track class populations across iterations of a RELION Class3D job:

```bash
relion/class_distribution Class3D/job047
```

### compare_class3d_jobs

Compare particle class assignments between two Class3D data STAR files:

```bash
relion/compare_class3d_jobs \
    Class3D/job020/run_it025_data.star \
    Class3D/job025/run_it022_data.star
```

### split_star_by_tomo

Split a RELION STAR file into one STAR file per tomogram or micrograph.
By default, `rlnTomoName` is used; if absent, `rlnMicrographName` is
used as a fallback:

```bash
relion/split_star_by_tomo run_data.star
relion/split_star_by_tomo run_data.star -o by_tomo
```

Explicitly group by `rlnMicrographName`:

```bash
relion/split_star_by_tomo run_data.star --group-by rlnMicrographName
```

Strip the `.tomostar` suffix from grouping-column values in the output
files (useful for ArtiaX):

```bash
relion/split_star_by_tomo run_data.star --strip-tomostar-suffix
```

### merge_star_files

Merge compatible STAR files (e.g. those produced by `split_star_by_tomo`)
back into a single file. Non-particle blocks such as `data_optics` are
preserved when they are identical across all inputs:

```bash
relion/merge_star_files \
    --input \
    by_tomo/L3_pos01_ts_004.star \
    by_tomo/L3_pos01_ts_007.star \
    by_tomo/L3_pos01_ts_006.star \
    --output merged_L3.star
```

### missali_report

Summarise MissAlignment alignment-loss results:

```bash
missalignment/missali_report tiltseries/
```

### Warp/M — Refinement Resolution History

`m_resolution` summarises how the global resolution of a Warp/M species
has changed across refinement versions. It is read-only and uses only the
Python standard library.

The path may be an M project directory, a species directory, or a single
`.species` file (default: current directory). The tool reads the active
`.species` XML and the archived copies under `versions/`, then rebuilds
the version ancestry by following `PreviousVersion` pointers back from the
active `Version`. Ordering comes from these pointers, not from filenames,
timestamps, or resolution.

For each version in the active lineage it reports `GlobalResolution` and
`DenoisedAtResolution` (Å), the change from the previous step, and a trend
(`better` / `worse` / `same`). The best historical version is the one with
the lowest `GlobalResolution` in the active lineage. Archived versions
outside this lineage are counted, and listed with `--all`. Missing
ancestors or cycles are reported as warnings.

```bash
m/m_resolution M_FV3_Nucleosome_job028_compg022_V100
m/m_resolution . --species FV3_Nucleosome_Core_299c059f
m/m_resolution . --all --show-xml-paths
m/m_resolution . --csv m_history.csv
```

The positional `path` is the input (M project, species directory, or
`.species` file). `--show-xml-paths` is an output option: it prints the
source `.species` XML file under each version row.

Other options: `--species` (choose a species when the project has more
than one), `--all` (also list versions outside the active lineage),
`--csv FILE` (export the table; existing files are never overwritten).

### Warp/M — Progressive Refinement

`m_progressive` prepares inputs from a completed RELION Refine3D job and
generates a shell script for progressive Warp/M refinement.

```bash
m/m_progressive Refine3D/job001 \
    --settings warp/tiltseries.settings \
    --pre-m-xml snapshots/pre_m_xml \
    --out m_runs/progressive_01 \
    --mask masks/binary_mask.mrc
```

Use a verified original Pre-M XML snapshot and a new output directory.
Omit `--mask` to generate one with RELION. Inspect the mask and generated
script before running refinement in your Warp/M environment on a GPU node
or workstation:

```bash
bash m_runs/progressive_01/run_m_progressive.sh
# To select GPU 1 instead of the default GPU 0:
GPU_ID=1 bash m_runs/progressive_01/run_m_progressive.sh
```

No SLURM submission is required. Keep generated outputs outside Git because
they contain local paths and research metadata. See `--help` for all options.

## License

See individual scripts for authorship and attribution.
