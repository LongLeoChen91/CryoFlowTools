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

## Installation

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

The scripts are standalone executables. Make them executable if needed:

```bash
chmod +x hpc/gpuview relion/class_distribution relion/compare_class3d_jobs relion/split_star_by_tomo relion/merge_star_files missalignment/missali_report
```

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

## License

See individual scripts for authorship and attribution.
