# CryoFlowTools

A small collection of command-line utilities for cryo-ET and RELION
workflows on HPC clusters.

## Tools

| Tool | Location | Description |
|------|----------|-------------|
| `gpuview` | `hpc/` | Summarises Slurm GPU-node resources and running jobs |
| `class_distribution` | `relion/` | Tracks RELION Class3D class distributions across iterations |
| `compare_class3d_jobs` | `relion/` | Compares particle class assignments between two RELION Class3D `*_data.star` files |
| `split_star_by_tomo` | `relion/` | Splits a RELION STAR file into one file per `rlnTomoName` or `rlnMicrographName` |
| `merge_star_files` | `relion/` | Merges compatible RELION STAR files while preserving particle metadata and other blocks |
| `missali_report` | `missalignment/` | Summarises MissAlignment `*_alignment_loss.json` results |
| `m_progressive` | `m/` | Prepares RELION-to-M inputs and generates a progressive refinement shell script |
| `m_resolution` | `m/` | Summarises Warp/M species refinement version ancestry and resolution history |

## Installation

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

The scripts are standalone executables. Make them executable if needed:

```bash
chmod +x hpc/gpuview relion/class_distribution relion/compare_class3d_jobs relion/split_star_by_tomo relion/merge_star_files missalignment/missali_report m/m_resolution m/m_progressive
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

`m_progressive` prepares inputs from a completed RELION Refine3D job and generates
`run_m_progressive.sh` for nine stages of M setup and progressive refinement.
It copies and checksums original Pre-M XML metadata, adjusts Warp settings,
links tilt stacks, validates a binary mask, and writes `input_summary.json`.
Preparation does not run M refinement and refuses an existing output directory.

```bash
m/m_progressive Refine3D/job001 \
    --settings warp/tiltseries.settings \
    --pre-m-xml snapshots/pre_m_xml \
    --out m_runs/progressive_01 \
    --mask masks/binary_mask.mrc
```

Use a completed job containing `RELION_JOB_EXIT_SUCCESS`, `run_data.star`,
`run_class001.mrc`, and both `run_half{1,2}_class001_unfil.mrc` files. Supply a
verified original Pre-M XML snapshot with names matching the tilt-stack folders.
Omit `--mask` to generate one with RELION; mask parameters are dataset-dependent.
Use `--sym C1 --diameter 200` if these cannot be read from `note.txt`, and
`--tiltstack PATH` to override tilt-stack discovery. See `--help` for all options.

Inspect the generated script, mask, and summary before executing on a GPU node
or workstation. No SLURM submission or environment variables are required.
Choose GPU 0 (default) or GPU 1 for a fresh output:

```bash
GPU_ID=0 bash m_runs/progressive_01/run_m_progressive.sh
# Alternatively, run that fresh output on GPU 1:
GPU_ID=1 bash m_runs/progressive_01/run_m_progressive.sh
```

`GPU_ID` selects MCore's device among the GPUs visible to the process. To run
only setup and baseline, then continue the remaining stages:

```bash
GPU_ID=0 START_STEP=1 END_STEP=4 bash m_runs/progressive_01/run_m_progressive.sh
GPU_ID=0 START_STEP=5 END_STEP=9 bash m_runs/progressive_01/run_m_progressive.sh
```

Optional runtime overrides include `CPU_MEMORY=1` and `CTF_BATCH=4`. Completed
steps are marked; inspect partial results before resuming after a failure.
Keep generated outputs outside Git: they contain local paths and research XML.
STAR matching is a best-effort preview; verify source/particle matching in MTools.
End-to-end M refinement validation is still in progress.

Basic synthetic checks (no refinement software or GPUs):
`python3 -m unittest discover -s tests -v`.

## License

See individual scripts for authorship and attribution.
