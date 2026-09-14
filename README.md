# CryoFlowTools

A small collection of command-line utilities for cryo-ET and RELION
workflows on HPC clusters.

## Tools

| Tool | Location | Description |
|------|----------|-------------|
| `gpuview` | `hpc/` | Summarises Slurm GPU-node resources and running jobs |
| `class_distribution` | `relion/` | Tracks RELION Class3D class distributions across iterations |
| `compare_class3d_jobs` | `relion/` | Compares particle class assignments between two RELION Class3D `*_data.star` files |
| `split_star_by_tomo` | `relion/` | Splits a RELION STAR file into one file per tomogram (`rlnTomoName`) |
| `missali_report` | `missalignment/` | Summarises MissAlignment `*_alignment_loss.json` results |

## Installation

Install the Python dependencies:

```bash
pip install -r requirements.txt
```

The scripts are standalone executables. Make them executable if needed:

```bash
chmod +x hpc/gpuview relion/class_distribution relion/compare_class3d_jobs relion/split_star_by_tomo missalignment/missali_report
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

Split a RELION STAR file into one STAR file per tomogram:

```bash
relion/split_star_by_tomo run_data.star
relion/split_star_by_tomo run_data.star -o by_tomo
```

### missali_report

Summarise MissAlignment alignment-loss results:

```bash
missalignment/missali_report tiltseries/
```

## License

See individual scripts for authorship and attribution.
