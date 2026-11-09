# CryoFlowTools

A small collection of command-line utilities for cryo-ET and RELION
workflows on HPC clusters.

## Tools

| Tool | Location | Description |
|------|----------|-------------|
| `gpuview` | `hpc/` | Summarises Slurm GPU-node resources and running jobs |
| `nodeplan` | `hpc/` | Read-only Slurm node availability check: state, free GPUs, blockers, and predicted schedule |
| `jobdiag` | `hpc/` | Read-only Slurm job diagnostics for pending reasons, resource fit, and scheduler state |
| `missali_report` | `missalignment/` | Summarises MissAlignment `*_alignment_loss.json` results |
| `class_distribution` | `relion/` | Tracks RELION Class3D class distributions across iterations |
| `compare_class3d_jobs` | `relion/` | Compares particle class assignments between two RELION Class3D `*_data.star` files |
| `class_pool_merge` | `relion/` | Validates and merges manually assigned Class3D pools with frozen-pool preservation and particle lineage |
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
chmod +x hpc/gpuview hpc/nodeplan hpc/jobdiag relion/class_distribution relion/compare_class3d_jobs relion/class_pool_merge relion/split_star_by_tomo relion/merge_star_files missalignment/missali_report m/m_resolution m/m_progressive
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

### nodeplan

Check whether a node can be used now, and if not, why:

```bash
hpc/nodeplan compg027
hpc/nodeplan compg037 --verbose
```

The default output is a short summary: node state, free GPUs, the main blocker (`DRAIN`, `DOWN`, `PLANNED`, or fully allocated), the next job Slurm predicts for the node, and an availability verdict of `NOW`, `NO`, or `UNKNOWN`. `--verbose` adds node resources, running jobs, planned jobs, and the commands executed.

A node name Slurm does not define is reported as `NOT DEFINED` (exit code 2), which is different from a drained or down node; `sinfo` omits nodes that belong to no partition, so it cannot confirm a name.

`nodeplan` is strictly read-only and uses only `scontrol show` and `squeue`.

> Start times are Slurm estimates that change every scheduling cycle, not guarantees. `nodeplan` never treats `PLANNED` as a reservation window and reports availability as `UNKNOWN` unless the evidence establishes it.

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

### class_pool_merge

Route manually selected RELION Class3D classes into configurable particle pools
across multiple jobs, retaining frozen pools and cumulative particle lineage.
Requires Python 3.11+, `pandas`, and `starfile` from the existing requirements.
The tool validates the researcher's assignments; it does not interpret class morphology.

Save a TOML configuration such as `round2.toml`, using your actual project paths:

```toml
round = 2
previous_round = 1
project_dir = "/path/to/relion/project"
output_dir = "Round2_split"
pools = ["Good", "Junk"]
id_column = "rlnImageName"

[jobs.job002]
star = "Class3D/job002/run_it025_data.star"
input_star = "previous_input.star"
input_pool = "Mix"
k = 2

[jobs.job002.classes]
Good = [1]
Junk = [2]
```

Add `[jobs.<name>]` sections for additional jobs. Every class in `1..k`, including
zero-particle classes, must belong to exactly one pool. Assignments accept integer
lists or strings such as `"1-3, 7"`; `[]` and `""` assign no classes.
The model path is derived from `*_data.star`, or supplied as `model_star`.
Optional `expected_particles` on jobs/frozen pools and top-level `expected_total`
are checked against measured counts. Paths resolve against `project_dir`, which
defaults to the current directory; `--project-dir` and `--outdir` override it and
the destination respectively.

```bash
python relion/class_pool_merge --config round2.toml --dry-run
python relion/class_pool_merge --config round2.toml
```

Dry-run is read-only, including when the destination already exists. Normal
execution requires each job's `input_star` to match its complete output identity
set, regardless of row order. If evidence is unavailable, execution refuses by
default; `--allow-unverified-inputs` explicitly permits those missing checks and
records the incomplete verification in reports. It never bypasses a failed check
on a configured input. No existing result directory is overwritten.

To carry a pool forward, declare its name in `pools` and add, for example:

```toml
[frozen.Archive]
star = "Round1_split/particles_Round1_Archive_250.star"
mode = "copy_only"
expected_particles = 250
```

`copy_only` preserves that STAR byte-for-byte and forbids new class assignments
to the pool. Explicit `mode = "append"` allows new assignments: existing frozen
rows come first, followed by selected rows in configured job order, preserving
each source's row order. Set optional `input_pool` when carrying a previous-round
pool into a differently named output pool. Multiple frozen pools are supported.
Their identities must be unique and disjoint from every other active or frozen source.

Outputs include `particles_<round>_<pool>_<count>.star` for every pool, including
header-only STAR tables for empty pools, plus `manifest.json`, `summary.txt`,
the effective `config.toml`, `provenance.json`, and `provenance.csv`. The manifest
records input paths, assignments, counts, software versions, parameters, and
verification results. JSON provenance uses schema version 1 with identity tuples
and per-particle history events; CSV provides cumulative flat lineage fields.

Set top-level `previous_provenance` to a prior JSON provenance file or historical
CSV, and supply `previous_round`. Historical CSV must contain the selected identity
columns and `<previous_round_prefix>_pool` (for example, `round1_pool`). Every
particle must be covered exactly once, with pool labels consistent with the input
sources. All previous fields are retained and the original provenance file is
archived alongside the reports. Unsupported formats and conflicting lineage fail
explicitly; absent previous provenance is recorded as unavailable.

Before publication, the tool checks optics-group references, compatible optics
definitions and other STAR blocks, output particle schemas, and exact identity-set
conservation. Compatible optics groups may be combined without renumbering; other
metadata comparisons allow harmless row ordering and tight floating-point
serialization differences. Mixed optics/no-optics inputs are rejected. Inputs
without optics are supported when none references an optics group. Files are
written to a sibling staging directory, read back to verify their identities,
counts and metadata, then published atomically without replacing an existing
destination. Where the filesystem supports it, publication uses a no-replace
directory rename. On BMRC GPFS, the final path is an atomically created directory
symlink to a hidden sibling `.<output>.data-<suffix>` containing the complete
result. The CLI reports this fallback explicitly. Atomic publication means that
readers see either no result or the complete result; it does not guarantee
persistence after power loss or a filesystem crash.

For GPFS results, keep the relative link and its hidden backing directory together
in the same parent when moving them. To make an independent copy, copy the contents
through the link into a new destination, for example
`cp -aL Round2_split/. /path/to/new-copy/`. To archive only the result, archive its
contents by changing into it, for example
`tar -czf /path/to/round2.tar.gz -C Round2_split .`; archiving the symlink alone
does not include its data. Check the copied or archived contents before removing
the source. To delete a result, first record and verify the target with
`readlink Round2_split`, remove the link without a trailing slash, then remove
only that verified backing directory when no longer needed. Deleting only the
link leaves the data; deleting only the backing directory leaves a broken link.

Caught failures clean up unpublished storage. Abrupt termination (such as SIGKILL)
or node failure can leave an unpublished hidden directory. Inspect such directories
before cleanup, and never delete storage referenced by a result link. Processing
is in memory; plan RAM for the source tables, selected rows, and provenance.

Choose an explicit identity known to be stable in your workflow, such as
`rlnImageName` for an appropriate subtomogram workflow, or a composite list such
as `["rlnTomoName", "rlnTomoParticleId"]`. `id_column = "auto"` requires exactly
one common named candidate across all inputs; it does not establish stability
outside those files. IDs retain exact string spelling, including leading zeros;
paths and coordinates are not normalized. Coordinate identities require explicit
selection. Embedded quotes, unquoted hashes, and multiline particle tokens that
the installed parser cannot preserve safely are rejected. Repeated block names in
particle STAR files are rejected before parsing can discard particles or metadata.
Particle column sets must agree within each merged pool; no columns are silently
dropped.

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
