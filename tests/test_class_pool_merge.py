"""Synthetic tests for class_pool_merge config, STAR validation, and recovery."""

from __future__ import annotations

from dataclasses import replace
import errno
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import textwrap
from typing import Any
import unittest
from unittest import mock


CLI = Path(__file__).resolve().parents[1] / "relion" / "class_pool_merge"
MODULE = runpy.run_path(str(CLI))
load_config = MODULE["load_config"]
parse_class_spec = MODULE["parse_class_spec"]
parse_model_star = MODULE["parse_model_star"]
validate_job_star = MODULE["validate_job_star"]
validate_star_inputs = MODULE["validate_star_inputs"]
format_validation_summary = MODULE["format_validation_summary"]


def write_synthetic_model_star(
    path: Path,
    distributions: dict[int, float],
    nr_classes_header: int | None = None,
) -> None:
    """Write a minimal RELION Class3D run_it*_model.star file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    k = nr_classes_header if nr_classes_header is not None else len(distributions)
    lines = [
        "# version 50001",
        "",
        "data_model_general",
        "",
        f"_rlnNrClasses {k}",
        "",
        "data_model_classes",
        "",
        "loop_",
        "_rlnReferenceImage #1",
        "_rlnClassDistribution #2",
        "_rlnEstimatedResolution #3",
    ]
    for cls_id, dist in distributions.items():
        res = "inf" if dist == 0.0 else "28.500000"
        lines.append(
            f"Class3D/job/run_it025_class{cls_id:03d}.mrc {dist:.6f} {res}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_synthetic_data_star(
    path: Path,
    class_numbers: list[Any],
    image_prefix: str = "subtomo/tomo01/part",
) -> None:
    """Write a minimal RELION Class3D run_it*_data.star file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# version 50001",
        "",
        "data_optics",
        "",
        "loop_",
        "_rlnOpticsGroup #1",
        "_rlnOpticsGroupName #2",
        "_rlnImagePixelSize #3",
        "1 opticsGroup1 10.000000",
        "",
        "data_particles",
        "",
        "loop_",
        "_rlnImageName #1",
        "_rlnClassNumber #2",
        "_rlnOpticsGroup #3",
    ]
    for idx, cls_val in enumerate(class_numbers, start=1):
        lines.append(f"{image_prefix}_{idx:05d}.mrc {cls_val} 1")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_particle_table(path: Path, columns: list[str], rows: list[list[Any]]) -> None:
    """Write literal synthetic identifiers, including STAR null/empty values."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["data_particles", "", "loop_"]
    lines.extend(f"_{column} #{index}" for index, column in enumerate(columns, 1))
    for row in rows:
        tokens = []
        for value in row:
            token = "?" if value is None else str(value)
            if not token or any(char.isspace() for char in token):
                token = f'"{token}"'
            tokens.append(token)
        lines.append(" ".join(tokens))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class ClassSpecTests(unittest.TestCase):
    """Test compact range string and integer list class specifications."""

    def test_valid_integer_lists_and_empty_specs(self):
        self.assertEqual(parse_class_spec([]), [])
        self.assertEqual(parse_class_spec(""), [])
        self.assertEqual(parse_class_spec("   "), [])
        self.assertEqual(parse_class_spec([1, 4, 6, 12]), [1, 4, 6, 12])

    def test_valid_compact_range_strings(self):
        self.assertEqual(
            parse_class_spec("1-5, 7, 9-12"),
            [1, 2, 3, 4, 5, 7, 9, 10, 11, 12],
        )
        self.assertEqual(
            parse_class_spec("3-3, 8 - 10"),
            [3, 8, 9, 10],
        )

    def test_invalid_class_specs_rejected(self):
        invalid_specs = [
            None,
            12,
            [0, 1, 2],
            [-1, 2],
            [1, "2"],
            [True, 2],
            [1, 2, 2],
            "0-3",
            "5-2",
            "1-4, 3",
            "1-3, , 5",
            "1-3,",
            "1..3",
            "class1",
        ]
        for spec in invalid_specs:
            with self.subTest(spec=spec):
                with self.assertRaises(ValueError):
                    parse_class_spec(spec)


class ConfigValidationTests(unittest.TestCase):
    """Test TOML configuration parsing and validation rules."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="class_pool_merge_2a_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_toml(self, content: str, name: str = "config.toml") -> Path:
        path = self.root / name
        path.write_text(textwrap.dedent(content).lstrip(), encoding="utf-8")
        return path

    def test_valid_multi_job_and_frozen_pools_config(self):
        cfg_path = self.write_toml(
            """
            round = "Round15"
            previous_round = 14
            project_dir = "/tmp/relion_proj"
            output_dir = "Round15_split"
            pools = ["Nucleosome", "Mix", "Junk", "Discard", "Stacking"]
            id_column = "rlnImageName"
            expected_total = 147888

            [frozen.Nucleosome]
            star = "Round14_split/particles_Round14_Nucleosome_50099.star"
            mode = "copy_only"
            expected_particles = 50099

            [jobs.job046]
            star = "Class3D/job046/run_it025_data.star"
            input_pool = "Mix"
            expected_particles = 23905
            k = 16

            [jobs.job046.classes]
            Mix = "2-4, 6-10, 12-14, 16"
            Junk = [1, 5, 11, 15]
            Discard = []
            Stacking = ""

            [jobs.job047]
            star = "Class3D/job047/run_it025_data.star"
            input_pool = "Junk"
            expected_particles = 14324
            k = 32

            [jobs.job047.classes]
            Mix = [6, 21, 25, 30]
            Junk = [1, 4, 7, 9, 14, 17, 19, 22, 28, 29, 31]
            Discard = [2, 3, 5, 8, 10, 11, 12, 13, 15, 16, 18, 20, 23, 24, 26, 27, 32]
            Stacking = []

            [jobs.job048]
            star = "Class3D/job048/run_it025_data.star"
            input_pool = "Discard"
            expected_particles = 59560
            k = 32
            note = "Job048 Class12 has weak stacking-like feature, manually assigned Junk"

            [jobs.job048.classes]
            Mix = []
            Junk = [6, 12, 22, 29]
            Discard = "4, 5, 7-11, 13-21, 24-28, 30-32"
            Stacking = [1, 2, 3, 23]
            """
        )
        cfg = load_config(cfg_path)
        self.assertEqual(cfg.round_label, "Round15")
        self.assertEqual(cfg.round_prefix, "round15")
        self.assertEqual(cfg.previous_round, "Round14")
        self.assertEqual(
            cfg.pools,
            ("Nucleosome", "Mix", "Junk", "Discard", "Stacking"),
        )
        self.assertEqual(cfg.expected_total, 147888)
        self.assertEqual(cfg.frozen["Nucleosome"].mode, "copy_only")
        self.assertEqual(cfg.frozen["Nucleosome"].source_label, "FROZEN_ROUND14")
        self.assertEqual(len(cfg.jobs), 3)
        self.assertEqual(cfg.jobs["job046"].classes["Discard"], [])
        self.assertEqual(cfg.jobs["job046"].classes["Stacking"], [])
        self.assertEqual(cfg.jobs["job048"].classes["Stacking"], [1, 2, 3, 23])

    def test_valid_frozen_append_mode_allows_active_assignments(self):
        cfg_path = self.write_toml(
            """
            round = 16
            pools = ["Nucleosome", "Mix", "Junk"]

            [frozen.Nucleosome]
            star = "Round15_split/particles_Round15_Nucleosome_50099.star"
            mode = "append"
            expected_particles = 50099

            [jobs.job050]
            star = "Class3D/job050/run_it025_data.star"
            k = 8
            expected_particles = 10000

            [jobs.job050.classes]
            Nucleosome = [1, 2]
            Mix = "3-5"
            Junk = "6-8"
            """
        )
        cfg = load_config(cfg_path)
        self.assertEqual(cfg.round_label, "Round16")
        self.assertEqual(cfg.frozen["Nucleosome"].mode, "append")
        self.assertEqual(cfg.jobs["job050"].classes["Nucleosome"], [1, 2])

    def test_copy_only_frozen_pool_rejects_active_assignments(self):
        cfg_path = self.write_toml(
            """
            round = "Round14"
            pools = ["Nucleosome", "Mix", "Junk"]

            [frozen.Nucleosome]
            star = "Round13_split/particles_Round13_Nucleosome_50099.star"
            mode = "copy_only"

            [jobs.job042]
            star = "Class3D/job042/run_it025_data.star"
            k = 4

            [jobs.job042.classes]
            Nucleosome = [1]
            Mix = [2, 3]
            Junk = [4]
            """
        )
        with self.assertRaisesRegex(ValueError, "copy_only"):
            load_config(cfg_path)

    def test_missing_required_fields_rejected(self):
        cases = {
            "missing_round": """
                pools = ["Good", "Junk"]
                [jobs.job001]
                star = "a.star"
                k = 2
                [jobs.job001.classes]
                Good = [1]
                Junk = [2]
            """,
            "missing_pools": """
                round = "Round1"
                [jobs.job001]
                star = "a.star"
                k = 2
                [jobs.job001.classes]
                Good = [1]
                Junk = [2]
            """,
            "missing_jobs": """
                round = "Round1"
                pools = ["Good", "Junk"]
            """,
            "missing_job_star": """
                round = "Round1"
                pools = ["Good", "Junk"]
                [jobs.job001]
                k = 2
                [jobs.job001.classes]
                Good = [1]
                Junk = [2]
            """,
            "missing_job_k": """
                round = "Round1"
                pools = ["Good", "Junk"]
                [jobs.job001]
                star = "a.star"
                [jobs.job001.classes]
                Good = [1]
                Junk = [2]
            """,
            "missing_job_classes": """
                round = "Round1"
                pools = ["Good", "Junk"]
                [jobs.job001]
                star = "a.star"
                k = 2
            """,
        }
        for name, toml_str in cases.items():
            with self.subTest(case=name):
                path = self.write_toml(toml_str, f"{name}.toml")
                with self.assertRaises(ValueError):
                    load_config(path)

    def test_duplicate_pools_and_unknown_target_pools_rejected(self):
        dup_pools = self.write_toml(
            """
            round = "Round1"
            pools = ["Mix", "Junk", "Mix"]
            [jobs.job001]
            star = "a.star"
            k = 2
            [jobs.job001.classes]
            Mix = [1]
            Junk = [2]
            """,
            "dup_pools.toml",
        )
        with self.assertRaisesRegex(ValueError, "Duplicate pool names"):
            load_config(dup_pools)

        unknown_job_pool = self.write_toml(
            """
            round = "Round1"
            pools = ["Mix", "Junk"]
            [jobs.job001]
            star = "a.star"
            k = 2
            [jobs.job001.classes]
            Mix = [1]
            TypoPool = [2]
            """,
            "unknown_job_pool.toml",
        )
        with self.assertRaisesRegex(ValueError, "unknown target pool"):
            load_config(unknown_job_pool)

        unknown_frozen_pool = self.write_toml(
            """
            round = "Round1"
            pools = ["Mix", "Junk"]
            [frozen.UnknownFrozen]
            star = "frozen.star"
            [jobs.job001]
            star = "a.star"
            k = 2
            [jobs.job001.classes]
            Mix = [1]
            Junk = [2]
            """,
            "unknown_frozen_pool.toml",
        )
        with self.assertRaisesRegex(ValueError, "not declared in 'pools'"):
            load_config(unknown_frozen_pool)

    def test_invalid_frozen_mode_rejected(self):
        bad_mode = self.write_toml(
            """
            round = "Round1"
            pools = ["Nucleosome", "Mix", "Junk"]
            [frozen.Nucleosome]
            star = "frozen.star"
            mode = "auto_switch"
            [jobs.job001]
            star = "a.star"
            k = 2
            [jobs.job001.classes]
            Mix = [1]
            Junk = [2]
            """,
            "bad_mode.toml",
        )
        with self.assertRaisesRegex(ValueError, "mode must be one of"):
            load_config(bad_mode)

    def test_missing_duplicate_and_out_of_range_classes_rejected(self):
        missing_cls = self.write_toml(
            """
            round = "Round1"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "a.star"
            k = 5
            [jobs.job001.classes]
            Good = "1-2"
            Junk = [4, 5]
            """,
            "missing_cls.toml",
        )
        with self.assertRaisesRegex(ValueError, "missing class assignments"):
            load_config(missing_cls)

        dup_across_pools = self.write_toml(
            """
            round = "Round1"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "a.star"
            k = 4
            [jobs.job001.classes]
            Good = "1-3"
            Junk = "3-4"
            """,
            "dup_across_pools.toml",
        )
        with self.assertRaisesRegex(ValueError, "duplicate class assignments across pools"):
            load_config(dup_across_pools)

        out_of_range = self.write_toml(
            """
            round = "Round1"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "a.star"
            k = 4
            [jobs.job001.classes]
            Good = "1-3"
            Junk = [4, 5]
            """,
            "out_of_range.toml",
        )
        with self.assertRaisesRegex(ValueError, "out of range"):
            load_config(out_of_range)

    def test_expected_total_mismatch_rejected(self):
        bad_total = self.write_toml(
            """
            round = "Round1"
            pools = ["Good", "Junk"]
            expected_total = 500
            [jobs.job001]
            star = "a.star"
            k = 2
            expected_particles = 300
            [jobs.job001.classes]
            Good = [1]
            Junk = [2]
            """,
            "bad_total.toml",
        )
        with self.assertRaisesRegex(ValueError, "does not match expected_total"):
            load_config(bad_total)

    def test_cli_execution_and_overrides(self):
        write_synthetic_model_star(
            self.root / "Class3D/job039/run_it025_model.star",
            {1: 0.4, 2: 0.3, 3: 0.3, 4: 0.0},
        )
        write_synthetic_data_star(
            self.root / "Class3D/job039/run_it025_data.star",
            [1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        )
        cfg_path = self.write_toml(
            """
            round = "Round13"
            pools = ["Nucleosome", "Mix", "Junk", "Discard"]
            [jobs.job039]
            star = "Class3D/job039/run_it025_data.star"
            k = 4
            [jobs.job039.classes]
            Nucleosome = [1, 2]
            Mix = [3]
            Junk = [4]
            Discard = []
            """,
            "cli_ok.toml",
        )
        res = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--config",
                str(cfg_path),
                "--project-dir",
                str(self.root),
                "--outdir",
                str(self.root / "custom_out"),
                "--dry-run",
            ],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertIn(
            "Configuration and STAR input validation succeeded (--dry-run).",
            res.stdout,
        )
        self.assertIn("custom_out", res.stdout)
        self.assertIn("Checks not yet performed:", res.stdout)
        self.assertFalse((self.root / "custom_out").exists())

        # Execution requires previous-input evidence or an explicit override.
        no_dry_res = subprocess.run(
            [
                sys.executable,
                str(CLI),
                "--config",
                str(cfg_path),
                "--project-dir",
                str(self.root),
            ],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(no_dry_res.returncode, 0)
        self.assertIn("input_star", no_dry_res.stderr)
        self.assertFalse((self.root / "Round13_split").exists())

        explicit_out = self.root / "explicit_execution"
        write_res = subprocess.run(
            [sys.executable, str(CLI), "--config", str(cfg_path),
             "--project-dir", str(self.root), "--outdir", str(explicit_out),
             "--allow-unverified-inputs"],
            capture_output=True, text=True,
        )
        self.assertEqual(write_res.returncode, 0, write_res.stderr)
        self.assertTrue((explicit_out / "particles_Round13_Nucleosome_7.star").is_file())
        manifest = json.loads((explicit_out / "manifest.json").read_text())
        self.assertIs(manifest["execution"]["allow_unverified_inputs"], True)
        self.assertIsNone(manifest["validation"]["all_previous_inputs_verified"])

        bad_path = self.write_toml(
            """
            round = "Round13"
            pools = ["Good", "Junk"]
            [jobs.job039]
            star = "Class3D/job039/run_it025_data.star"
            k = 3
            [jobs.job039.classes]
            Good = [1]
            Junk = [2]
            """,
            "cli_bad.toml",
        )
        bad_res = subprocess.run(
            [sys.executable, str(CLI), "--config", str(bad_path)],
            capture_output=True,
            text=True,
        )
        self.assertNotEqual(bad_res.returncode, 0)
        self.assertIn("ERROR:", bad_res.stderr)


class StarInputValidationTests(unittest.TestCase):
    """Test Phase 2B-1A RELION data.star and model.star parsing and validation."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="class_pool_merge_2b1a_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_toml(self, content: str, name: str = "config.toml") -> Path:
        path = self.root / name
        path.write_text(textwrap.dedent(content).lstrip(), encoding="utf-8")
        return path

    def test_valid_multi_job_with_zero_particle_classes(self):
        write_synthetic_model_star(
            self.root / "Class3D/job046/run_it025_model.star",
            {1: 0.25, 2: 0.50, 3: 0.25, 4: 0.0},
        )
        write_synthetic_data_star(
            self.root / "Class3D/job046/run_it025_data.star",
            [1, 1, 2, 2, 2, 2, 3, 3],
            image_prefix="subtomo/tomo1/p",
        )

        write_synthetic_model_star(
            self.root / "Class3D/job047/run_it025_model.star",
            {1: 0.50, 2: 0.0, 3: 0.0, 4: 0.50},
        )
        write_synthetic_data_star(
            self.root / "Class3D/job047/run_it025_data.star",
            [1, 1, 1, 4, 4, 4],
            image_prefix="subtomo/tomo2/p",
        )
        write_synthetic_data_star(
            self.root / "Round14_split/particles_Round14_Nucleosome_10.star",
            [1] * 10,
            image_prefix="frozen/particle",
        )

        cfg_path = self.write_toml(
            f"""
            round = "Round15"
            project_dir = "{self.root}"
            pools = ["Nucleosome", "Mix", "Junk", "Discard", "Stacking"]
            expected_total = 24

            [frozen.Nucleosome]
            star = "Round14_split/particles_Round14_Nucleosome_10.star"
            mode = "copy_only"
            expected_particles = 10

            [jobs.job046]
            star = "Class3D/job046/run_it025_data.star"
            input_pool = "Mix"
            expected_particles = 8
            k = 4

            [jobs.job046.classes]
            Mix = [2, 3]
            Junk = [1]
            Discard = [4]

            [jobs.job047]
            star = "Class3D/job047/run_it025_data.star"
            input_pool = "Junk"
            expected_particles = 6
            k = 4

            [jobs.job047.classes]
            Mix = [1]
            Junk = []
            Discard = [2, 3]
            Stacking = [4]
            """
        )
        cfg = load_config(cfg_path)
        report = validate_star_inputs(cfg)

        self.assertEqual(report.active_particles, 14)
        self.assertEqual(report.frozen_expected_particles, 10)
        self.assertEqual(report.total_particles, 24)
        self.assertEqual(report.frozen["Nucleosome"].particle_count, 10)
        self.assertEqual(report.frozen_particles, 10)
        self.assertEqual(report.unique_particle_ids, 24)
        summary = format_validation_summary(cfg, report)
        self.assertIn("Frozen checks confirmed", summary)
        self.assertIn("Validated pool totals:", summary)
        self.assertEqual(report.jobs["job046"].empty_classes, (4,))
        self.assertEqual(report.jobs["job047"].empty_classes, (2, 3))
        self.assertEqual(
            report.pool_totals,
            {"Nucleosome": 10, "Mix": 9, "Junk": 2, "Discard": 0, "Stacking": 3},
        )

    def test_detected_k_mismatch_with_configured_k_rejected(self):
        write_synthetic_model_star(
            self.root / "Class3D/job001/run_it025_model.star",
            {1: 0.4, 2: 0.3, 3: 0.3, 4: 0.0, 5: 0.0},
        )
        write_synthetic_data_star(
            self.root / "Class3D/job001/run_it025_data.star",
            [1, 2, 3],
        )
        cfg_path = self.write_toml(
            f"""
            round = "Round1"
            project_dir = "{self.root}"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "Class3D/job001/run_it025_data.star"
            k = 4
            [jobs.job001.classes]
            Good = "1-2"
            Junk = "3-4"
            """
        )
        cfg = load_config(cfg_path)
        with self.assertRaisesRegex(ValueError, "configured k=4 does not match model STAR K=5"):
            validate_star_inputs(cfg)

    def test_model_header_nr_classes_disagreement_rejected(self):
        model_path = self.root / "Class3D/job001/run_it025_model.star"
        write_synthetic_model_star(
            model_path,
            {1: 0.5, 2: 0.5},
            nr_classes_header=3,
        )
        with self.assertRaisesRegex(ValueError, "_rlnNrClasses .* disagrees"):
            parse_model_star(model_path)

    def test_invalid_or_out_of_range_data_class_numbers_rejected(self):
        write_synthetic_model_star(
            self.root / "Class3D/job001/run_it025_model.star",
            {1: 0.5, 2: 0.5, 3: 0.0},
        )
        # Class 4 in data.star is outside K=3
        write_synthetic_data_star(
            self.root / "Class3D/job001/run_it025_data.star",
            [1, 2, 4],
        )
        cfg_path = self.write_toml(
            f"""
            round = "Round1"
            project_dir = "{self.root}"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "Class3D/job001/run_it025_data.star"
            k = 3
            [jobs.job001.classes]
            Good = [1, 2]
            Junk = [3]
            """
        )
        cfg = load_config(cfg_path)
        with self.assertRaisesRegex(ValueError, "outside 1..3"):
            validate_star_inputs(cfg)

        # Non-integer class number in data.star
        write_synthetic_data_star(
            self.root / "Class3D/job001/run_it025_data.star",
            [1, "2.5", 2],
        )
        with self.assertRaisesRegex(ValueError, "non-integer"):
            validate_star_inputs(cfg)

    def test_model_vs_data_empty_class_inconsistency_rejected(self):
        # Case A: model says class 3 has 0.0 distribution, but data.star has particles in class 3
        write_synthetic_model_star(
            self.root / "Class3D/job001/run_it025_model.star",
            {1: 0.5, 2: 0.5, 3: 0.0},
        )
        write_synthetic_data_star(
            self.root / "Class3D/job001/run_it025_data.star",
            [1, 2, 3],
        )
        cfg_path = self.write_toml(
            f"""
            round = "Round1"
            project_dir = "{self.root}"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "Class3D/job001/run_it025_data.star"
            k = 3
            [jobs.job001.classes]
            Good = [1, 2]
            Junk = [3]
            """
        )
        cfg = load_config(cfg_path)
        with self.assertRaisesRegex(ValueError, "0.0 distribution .* but contain particles"):
            validate_star_inputs(cfg)

        # Case B: data.star has 0 particles in class 2, but model says non-zero distribution
        write_synthetic_model_star(
            self.root / "Class3D/job001/run_it025_model.star",
            {1: 0.5, 2: 0.5, 3: 0.0},
        )
        write_synthetic_data_star(
            self.root / "Class3D/job001/run_it025_data.star",
            [1, 1, 1],
        )
        with self.assertRaisesRegex(ValueError, "0 particles .* but non-zero distribution"):
            validate_star_inputs(cfg)

    def test_particle_count_mismatch_rejected(self):
        write_synthetic_model_star(
            self.root / "Class3D/job001/run_it025_model.star",
            {1: 0.5, 2: 0.5},
        )
        write_synthetic_data_star(
            self.root / "Class3D/job001/run_it025_data.star",
            [1, 1, 2, 2],
        )
        # Job expected_particles mismatch (expected 5, actual 4)
        cfg_job_mismatch = self.write_toml(
            f"""
            round = "Round1"
            project_dir = "{self.root}"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "Class3D/job001/run_it025_data.star"
            k = 2
            expected_particles = 5
            [jobs.job001.classes]
            Good = [1]
            Junk = [2]
            """,
            "job_count_mismatch.toml",
        )
        with self.assertRaisesRegex(ValueError, "expected 5 particles, found 4"):
            validate_star_inputs(load_config(cfg_job_mismatch))

        # Global expected_total mismatch when per-job expected_particles is omitted
        cfg_total_mismatch = self.write_toml(
            f"""
            round = "Round1"
            project_dir = "{self.root}"
            pools = ["Good", "Junk"]
            expected_total = 10
            [jobs.job001]
            star = "Class3D/job001/run_it025_data.star"
            k = 2
            [jobs.job001.classes]
            Good = [1]
            Junk = [2]
            """,
            "total_count_mismatch.toml",
        )
        with self.assertRaisesRegex(ValueError, "Total validated particles .* does not match"):
            validate_star_inputs(load_config(cfg_total_mismatch))

    def test_missing_or_malformed_star_files_rejected(self):
        cfg_missing = self.write_toml(
            f"""
            round = "Round1"
            project_dir = "{self.root}"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "Class3D/job001/run_it025_data.star"
            k = 2
            [jobs.job001.classes]
            Good = [1]
            Junk = [2]
            """,
            "missing.toml",
        )
        with self.assertRaises(FileNotFoundError):
            validate_star_inputs(load_config(cfg_missing))

        # Write model.star without data.star
        write_synthetic_model_star(
            self.root / "Class3D/job001/run_it025_model.star",
            {1: 0.5, 2: 0.5},
        )
        with self.assertRaises(FileNotFoundError):
            validate_star_inputs(load_config(cfg_missing))

        # Write data.star missing rlnClassNumber column
        bad_data = self.root / "Class3D/job001/run_it025_data.star"
        bad_data.write_text(
            "data_particles\n\nloop_\n_rlnImageName #1\np1.mrc\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "missing required column 'rlnClassNumber'"):
            validate_star_inputs(load_config(cfg_missing))


class RecoveryValidationTests(unittest.TestCase):
    """Only checks actually performed may appear as verified results."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="class_pool_merge_recovery_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        write_synthetic_model_star(self.root / "run_model.star", {1: 0.5, 2: 0.5})
        write_synthetic_data_star(self.root / "run_data.star", [1, 2])
        config_path = self.root / "config.toml"
        config_path.write_text(textwrap.dedent(f"""
            round = "Round1"
            project_dir = "{self.root}"
            pools = ["Good", "Junk"]
            [jobs.job001]
            star = "run_data.star"
            k = 2
            [jobs.job001.classes]
            Good = [1]
            Junk = [2]
        """), encoding="utf-8")
        self.cfg = load_config(config_path)

    def test_identity_checks_run_while_unperformed_checks_remain_unknown(self):
        report = validate_star_inputs(self.cfg)
        self.assertEqual(report.active_particles, 2)
        self.assertEqual(report.total_particles, 2)
        self.assertEqual(report.frozen, {})
        self.assertEqual(report.frozen_particles, 0)
        self.assertEqual(report.id_columns, ("rlnImageName",))
        self.assertEqual(report.unique_particle_ids, 2)
        self.assertEqual(report.optics_groups, (1,))
        self.assertIs(report.has_optics_block, True)
        self.assertIsNone(report.all_previous_inputs_verified)
        self.assertIsNone(report.jobs["job001"].input_star_verified)
        self.assertEqual(report.jobs["job001"].optics_groups, (1,))
        summary = format_validation_summary(self.cfg, report)
        self.assertIn("Identity checks confirmed", summary)
        self.assertIn("optics", summary.lower())
        self.assertIn("not performed", summary.lower())

        # The formerly unverified duplicate and mismatched input now fail.
        data_path = self.cfg.jobs["job001"].star
        original = data_path.read_text(encoding="utf-8")
        data_path.write_text(
            original.replace("part_00002", "part_00001"),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "(?i)duplicate"):
            validate_star_inputs(self.cfg)
        data_path.write_text(original, encoding="utf-8")
        input_path = self.root / "previous.star"
        write_synthetic_data_star(input_path, [1], image_prefix="different/particle")
        job = replace(self.cfg.jobs["job001"], input_star=input_path)
        with self.assertRaisesRegex(ValueError, "(?i)missing.*unexpected"):
            validate_star_inputs(replace(self.cfg, jobs={job.name: job}))

    def test_frozen_expectations_are_checked_against_measured_counts(self):
        frozen_path = self.root / "frozen.star"
        write_synthetic_data_star(frozen_path, [1, 1, 2, 2], image_prefix="frozen/particle")
        for expected in (10, 4, None):
            with self.subTest(expected=expected):
                frozen = MODULE["FrozenPoolConfig"](
                    pool="Junk", star=frozen_path, mode="append", expected_particles=expected,
                )
                cfg = replace(self.cfg, frozen={"Junk": frozen})
                if expected == 10:
                    with self.assertRaisesRegex(ValueError, "(?i)expected 10.*found 4"):
                        validate_star_inputs(cfg)
                    continue
                report = validate_star_inputs(cfg)
                self.assertEqual(report.active_particles, 2)
                self.assertEqual(report.frozen_expected_particles, expected or 0)
                self.assertEqual(report.total_particles, 6)
                self.assertEqual(report.pool_totals, {"Good": 1, "Junk": 5})
                self.assertEqual(report.frozen["Junk"].particle_count, 4)
                self.assertEqual(report.frozen_particles, 4)
                self.assertEqual(report.unique_particle_ids, 6)
                summary = format_validation_summary(cfg, report)
                self.assertIn("Frozen checks confirmed", summary)
                self.assertIn("Validated pool totals:", summary)


class ParticleIdentityValidationTests(unittest.TestCase):
    """Identity conservation and disjointness using only synthetic STAR files."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="class_pool_merge_identity_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.job = self.make_job("job001", ["rlnImageName"], [["p1.mrc"], ["p2.mrc"]])
        self.cfg = MODULE["MergeConfig"](
            round_label="Round2", round_prefix="round2", project_dir=self.root,
            output_dir=self.root / "output", pools=("Good", "Old1", "Old2"),
            jobs={self.job.name: self.job}, id_column="rlnImageName",
        )

    def make_job(self, name, columns, rows):
        path = self.root / name / "run_data.star"
        write_particle_table(path, [*columns, "rlnClassNumber"], [[*row, 1] for row in rows])
        write_synthetic_model_star(path.with_name("run_model.star"), {1: 1.0})
        return MODULE["JobConfig"](name=name, star=path, k=1, classes={"Good": [1]})

    def with_job(self, columns, rows, **changes):
        job = self.make_job("job001", columns, rows)
        return replace(self.cfg, jobs={job.name: job}, **changes)

    def with_input(self, cfg, columns, rows):
        path = self.root / "previous.star"
        write_particle_table(path, columns, rows)
        job = replace(cfg.jobs["job001"], input_star=path)
        return replace(cfg, jobs={**cfg.jobs, job.name: job})

    def with_frozen(self, cfg, pool, ids, *, mode="copy_only", expected=None):
        path = self.root / f"{pool}.star"
        if mode == "append":
            write_particle_table(path, ["rlnImageName", "rlnClassNumber"], [[identity, 1] for identity in ids])
        else:
            write_particle_table(path, ["rlnImageName"], [[identity] for identity in ids])
        frozen = MODULE["FrozenPoolConfig"](
            pool=pool, star=path, mode=mode, expected_particles=expected,
        )
        return replace(cfg, frozen={**cfg.frozen, pool: frozen})

    def test_explicit_image_identity_overrides_other_candidates(self):
        cfg = self.with_job(
            ["rlnImageName", "rlnTomoParticleName"], [["p1.mrc", "1"], ["p2.mrc", "2"]],
        )
        report = validate_star_inputs(cfg)
        self.assertEqual(report.id_columns, ("rlnImageName",))
        self.assertEqual(report.unique_particle_ids, 2)

    def test_custom_identifiers_preserve_leading_zeros_and_numeric_text(self):
        rows = [[value] for value in ("0001", "1", "1e3", "1000", "NA")]
        cfg = self.with_job(["customParticleId"], rows, id_column="customParticleId")
        cfg = self.with_input(cfg, ["customParticleId"], list(reversed(rows)))
        report = validate_star_inputs(cfg)
        self.assertEqual(report.id_columns, ("customParticleId",))
        self.assertEqual(report.unique_particle_ids, 5)
        self.assertIs(report.jobs["job001"].input_star_verified, True)
        self.assertEqual(report.jobs["job001"].input_particle_count, 5)

    def test_embedded_quotes_fail_closed_but_star_quoted_plain_ids_work(self):
        # The installed parser rewrites apostrophes. Such IDs must not silently
        # match another spelling after parsing, whether quoted or unquoted.
        header = "data_particles\n\nloop_\n_rlnImageName #1\n_rlnClassNumber #2\n"
        for token in ("O'Brien.mrc", '"O\'Brien.mrc"', "'O\"Brien.mrc'"):
            with self.subTest(token=token):
                self.job.star.write_text(header + f"{token} 1\nplain.mrc 1\n", encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "(?i)quote|cannot preserve"):
                    validate_star_inputs(self.cfg)

        self.job.star.write_text(
            header + "'plain-one.mrc' 1\n\"plain-two.mrc\" 1\n", encoding="utf-8",
        )
        cfg = self.with_input(
            self.cfg, ["rlnImageName"], [["plain-two.mrc"], ["plain-one.mrc"]],
        )
        report = validate_star_inputs(cfg)
        self.assertEqual(report.unique_particle_ids, 2)
        self.assertIs(report.jobs["job001"].input_star_verified, True)

    def test_hash_tokens_cannot_be_truncated_and_quoted_hashes_match(self):
        # Put identity last: the parser can silently truncate this token without
        # losing a required class value, so ordinary class validation is insufficient.
        header = "data_particles\n\nloop_\n_rlnClassNumber #1\n_rlnImageName #2\n"
        self.job.star.write_text(header + "1 abc#one\n1 abc#two\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "unquoted '#'.*particle token"):
            validate_star_inputs(self.cfg)

        # An apostrophe rewritten as a closing quote can also expose a hash to
        # comment stripping and leave an apparently valid, truncated identity.
        self.job.star.write_text(header + '1 "abc\' #one"\n1 plain\n', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "embedded quote"):
            validate_star_inputs(self.cfg)

        self.job.star.write_text(
            header + "# whole-line comment\n1 'abc#one' # trailing comment\n1 \"abc#two\"\n",
            encoding="utf-8",
        )
        input_path = self.root / "previous.star"
        input_path.write_text(
            "data_particles\n\nloop_\n_rlnImageName #1\n"
            '"abc#two" # ordinary trailing comment\n\'abc#one\'\n',
            encoding="utf-8",
        )
        job = replace(self.job, input_star=input_path)
        report = validate_star_inputs(replace(self.cfg, jobs={job.name: job}))
        self.assertEqual(report.unique_particle_ids, 2)
        self.assertIs(report.jobs["job001"].input_star_verified, True)

        self.job.star.write_text(header + "1\n;multiline\nparticle\n;\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "unsupported multiline particle value"):
            validate_star_inputs(self.cfg)

    def test_missing_explicit_identity_column_rejected(self):
        cfg = replace(self.cfg, id_column="rlnTomoParticleName")
        with self.assertRaisesRegex(ValueError, "(?i)missing.*rlnTomoParticleName"):
            validate_star_inputs(cfg)

    def test_null_empty_and_whitespace_identifiers_rejected(self):
        for value in (None, ".", "", "   ", "NaN", "NULL"):
            with self.subTest(value=value):
                cfg = self.with_job(["rlnImageName"], [["p1.mrc"], [value]])
                with self.assertRaises(ValueError) as caught:
                    validate_star_inputs(cfg)
                self.assertIn("run_data.star", str(caught.exception))

    def test_duplicate_particle_ids_rejected(self):
        cfg = self.with_job(["rlnImageName"], [["p1.mrc"], ["p1.mrc"]])
        with self.assertRaisesRegex(ValueError, "(?i)duplicate") as caught:
            validate_star_inputs(cfg)
        self.assertIn("p1.mrc", str(caught.exception))

    def test_composite_identifiers_are_tuples_and_preserve_local_ids(self):
        columns = ["rlnTomoName", "rlnTomoParticleId"]
        rows = [["tomoA", "001"], ["tomoA", "1"], ["tomoB", "001"]]
        cfg = self.with_job(columns, rows, id_column=tuple(columns))
        cfg = self.with_input(cfg, columns, list(reversed(rows)))
        report = validate_star_inputs(cfg)
        self.assertEqual(report.id_columns, tuple(columns))
        self.assertEqual(report.unique_particle_ids, 3)
        self.assertIs(report.all_previous_inputs_verified, True)

    def test_composite_identity_loaded_from_toml(self):
        columns = ["rlnTomoName", "rlnTomoParticleId"]
        rows = [["tomoA", "001"], ["tomoA", "1"], ["tomoB", "001"]]
        self.make_job("job001", columns, rows)
        write_particle_table(self.root / "previous.star", columns, list(reversed(rows)))
        config_path = self.root / "composite.toml"
        config_path.write_text(textwrap.dedent(f"""
            round = "Round2"
            project_dir = "{self.root}"
            pools = ["Good"]
            id_column = ["rlnTomoName", "rlnTomoParticleId"]
            [jobs.job001]
            star = "job001/run_data.star"
            input_star = "previous.star"
            expected_particles = 3
            k = 1
            [jobs.job001.classes]
            Good = [1]
        """), encoding="utf-8")
        cfg = load_config(config_path)
        self.assertEqual(cfg.id_column, tuple(columns))
        report = validate_star_inputs(cfg)
        self.assertEqual(report.id_columns, tuple(columns))
        self.assertEqual(report.unique_particle_ids, 3)
        self.assertEqual(report.jobs["job001"].input_particle_count, 3)
        self.assertIs(report.all_previous_inputs_verified, True)

    def test_composite_keys_do_not_collide_when_components_contain_separators(self):
        columns = ["namespace", "localId"]
        rows = [["a:b", "c"], ["a", "b:c"], ["a|b", "c"], ["a", "b|c"]]
        cfg = self.with_job(columns, rows, id_column=tuple(columns))
        self.assertEqual(validate_star_inputs(cfg).unique_particle_ids, 4)

    def test_duplicate_and_incomplete_composite_identifiers_rejected(self):
        columns = ["rlnTomoName", "rlnTomoParticleId"]
        for rows in ([["tomoA", "1"], ["tomoA", "1"]], [["tomoA", "1"], ["tomoB", None]]):
            with self.subTest(rows=rows):
                cfg = self.with_job(columns, rows, id_column=tuple(columns))
                with self.assertRaises(ValueError):
                    validate_star_inputs(cfg)
        cfg = self.with_job(["rlnTomoName"], [["tomoA"]], id_column=tuple(columns))
        with self.assertRaisesRegex(ValueError, "(?i)missing.*rlnTomoParticleId"):
            validate_star_inputs(cfg)

    def test_auto_detection_requires_one_common_named_candidate(self):
        candidates = [
            (["rlnImageName"], [["p1.mrc"], ["p2.mrc"]]),
            (["rlnTomoParticleName"], [["tomoA/1"], ["tomoA/2"]]),
            (["rlnTomoName", "rlnTomoParticleId"], [["tomoA", "1"], ["tomoB", "1"]]),
        ]
        for columns, rows in candidates:
            with self.subTest(columns=columns):
                cfg = self.with_job(columns, rows, id_column="auto")
                self.assertEqual(validate_star_inputs(cfg).id_columns, tuple(columns))

    def test_auto_detection_rejects_ambiguous_candidates(self):
        for particle_names in (("one", "two"), ("same", "same")):
            with self.subTest(particle_names=particle_names):
                cfg = self.with_job(
                    ["rlnImageName", "rlnTomoParticleName"],
                    [["p1.mrc", particle_names[0]], ["p2.mrc", particle_names[1]]],
                    id_column="auto",
                )
                with self.assertRaisesRegex(ValueError, "id_column"):
                    validate_star_inputs(cfg)

    def test_auto_detection_does_not_guess_coordinate_identity(self):
        columns = ["rlnTomoName", "rlnCoordinateX", "rlnCoordinateY", "rlnCoordinateZ"]
        rows = [["tomoA", "1", "2", "3"], ["tomoA", "4", "5", "6"]]
        cfg = self.with_job(columns, rows, id_column="auto")
        with self.assertRaisesRegex(ValueError, "id_column"):
            validate_star_inputs(cfg)
        report = validate_star_inputs(replace(cfg, id_column=tuple(columns)))
        self.assertEqual(report.unique_particle_ids, 2)

    def test_auto_detection_requires_candidate_in_previous_and_frozen_sources(self):
        cfg = self.with_input(replace(self.cfg, id_column="auto"), ["customId"], [["1"], ["2"]])
        with self.assertRaisesRegex(ValueError, "id_column"):
            validate_star_inputs(cfg)
        cfg = self.with_frozen(replace(self.cfg, id_column="auto"), "Old1", ["f1.mrc"])
        write_particle_table(cfg.frozen["Old1"].star, ["customId"], [["f1"]])
        with self.assertRaisesRegex(ValueError, "id_column"):
            validate_star_inputs(cfg)

    def test_cross_job_overlap_rejected_with_source_names(self):
        second = self.make_job("job002", ["rlnImageName"], [["p2.mrc"], ["p3.mrc"]])
        cfg = replace(self.cfg, jobs={self.job.name: self.job, second.name: second})
        with self.assertRaisesRegex(ValueError, "(?i)overlap") as caught:
            validate_star_inputs(cfg)
        for value in ("job001", "job002", "p2.mrc"):
            self.assertIn(value, str(caught.exception))

    def test_previous_input_exact_conservation_ignores_order_and_class_column(self):
        cfg = self.with_input(self.cfg, ["rlnImageName"], [["p2.mrc"], ["p1.mrc"]])
        report = validate_star_inputs(cfg)
        self.assertEqual(report.active_particles, 2)
        self.assertIs(report.jobs["job001"].input_star_verified, True)
        self.assertEqual(report.jobs["job001"].input_particle_count, 2)
        self.assertIs(report.all_previous_inputs_verified, True)
        summary = format_validation_summary(cfg, report)
        self.assertIn("Previous input: identity match confirmed", summary)

    def test_previous_input_missing_and_unexpected_ids_reported(self):
        cfg = self.with_input(self.cfg, ["rlnImageName"], [["p1.mrc"], ["missing.mrc"]])
        with self.assertRaises(ValueError) as caught:
            validate_star_inputs(cfg)
        message = str(caught.exception)
        for value in ("job001", "missing", "unexpected", "missing.mrc", "p2.mrc"):
            self.assertIn(value, message)

    def test_previous_input_actual_count_mismatch_rejected(self):
        cfg = self.with_input(self.cfg, ["rlnImageName"], [["p1.mrc"]])
        with self.assertRaisesRegex(ValueError, "(?i)missing.*unexpected") as caught:
            validate_star_inputs(cfg)
        self.assertIn("p2.mrc", str(caught.exception))
        self.assertIn("1", str(caught.exception))
        self.assertIn("2", str(caught.exception))

    def test_previous_input_duplicates_and_missing_identifiers_rejected(self):
        cfg = self.with_input(self.cfg, ["rlnImageName"], [["p1.mrc"], ["p1.mrc"]])
        with self.assertRaisesRegex(ValueError, "(?i)duplicate"):
            validate_star_inputs(cfg)
        cfg = self.with_input(self.cfg, ["other"], [["p1.mrc"], ["p2.mrc"]])
        with self.assertRaisesRegex(ValueError, "(?i)missing.*rlnImageName"):
            validate_star_inputs(cfg)

    def test_missing_previous_input_file_rejected(self):
        job = replace(self.job, input_star=self.root / "missing_previous.star")
        with self.assertRaises(FileNotFoundError):
            validate_star_inputs(replace(self.cfg, jobs={job.name: job}))

    def test_partial_previous_verification_does_not_claim_all_jobs_checked(self):
        cfg = self.with_input(self.cfg, ["rlnImageName"], [["p2.mrc"], ["p1.mrc"]])
        second = self.make_job("job002", ["rlnImageName"], [["p3.mrc"]])
        cfg = replace(cfg, jobs={**cfg.jobs, second.name: second})
        report = validate_star_inputs(cfg)
        self.assertIs(report.jobs["job001"].input_star_verified, True)
        self.assertIsNone(report.jobs["job002"].input_star_verified)
        self.assertIsNone(report.all_previous_inputs_verified)
        self.assertIn("input_star not configured", format_validation_summary(cfg, report))

    def test_multiple_frozen_pools_and_append_use_actual_counts(self):
        cfg = self.with_frozen(self.cfg, "Old1", ["f1.mrc", "f2.mrc"], expected=2)
        cfg = self.with_frozen(cfg, "Old2", ["f3.mrc"])
        cfg = self.with_frozen(cfg, "Good", ["f4.mrc"], mode="append")
        report = validate_star_inputs(replace(cfg, expected_total=6))
        self.assertEqual(report.active_particles, 2)
        self.assertEqual(report.frozen_particles, 4)
        self.assertEqual(report.frozen_expected_particles, 2)
        self.assertEqual(report.total_particles, 6)
        self.assertEqual(report.unique_particle_ids, 6)
        self.assertEqual(report.pool_totals, {"Good": 3, "Old1": 2, "Old2": 1})
        self.assertEqual(report.frozen["Old1"].mode, "copy_only")
        self.assertEqual(report.frozen["Good"].mode, "append")
        self.assertTrue(all(item.optics_groups == () for item in report.frozen.values()))

    def test_frozen_count_mismatch_rejected(self):
        cfg = self.with_frozen(self.cfg, "Old1", ["f1.mrc"], expected=2)
        with self.assertRaisesRegex(ValueError, "(?i)expected 2.*found 1"):
            validate_star_inputs(cfg)

    def test_frozen_duplicates_and_invalid_identifiers_rejected(self):
        for ids in (["f1.mrc", "f1.mrc"], ["f1.mrc", None], ["f1.mrc", ""]):
            with self.subTest(ids=ids):
                cfg = self.with_frozen(self.cfg, "Old1", ids)
                with self.assertRaises(ValueError) as caught:
                    validate_star_inputs(cfg)
                self.assertIn("Old1", str(caught.exception))

    def test_active_frozen_overlap_rejected_in_both_modes(self):
        for mode, pool in (("copy_only", "Old1"), ("append", "Good")):
            with self.subTest(mode=mode):
                cfg = self.with_frozen(self.cfg, pool, ["p1.mrc"], mode=mode)
                with self.assertRaisesRegex(ValueError, "(?i)overlap") as caught:
                    validate_star_inputs(cfg)
                self.assertIn("job001", str(caught.exception))
                self.assertIn(pool, str(caught.exception))

    def test_frozen_pool_overlap_rejected(self):
        cfg = self.with_frozen(self.cfg, "Old1", ["f1.mrc"])
        cfg = self.with_frozen(cfg, "Old2", ["f1.mrc"])
        with self.assertRaisesRegex(ValueError, "(?i)overlap") as caught:
            validate_star_inputs(cfg)
        for value in ("Old1", "Old2", "f1.mrc"):
            self.assertIn(value, str(caught.exception))

    def test_missing_frozen_star_and_missing_id_column_rejected(self):
        cfg = self.with_frozen(self.cfg, "Old1", ["f1.mrc"])
        path = cfg.frozen["Old1"].star
        write_particle_table(path, ["other"], [["f1.mrc"]])
        with self.assertRaisesRegex(ValueError, "(?i)missing.*rlnImageName"):
            validate_star_inputs(cfg)
        path.unlink()
        with self.assertRaises(FileNotFoundError):
            validate_star_inputs(cfg)

    def test_global_expected_total_includes_frozen_without_expected_count(self):
        cfg = self.with_frozen(self.cfg, "Old1", ["f1.mrc"])
        with self.assertRaisesRegex(ValueError, "expected_total"):
            validate_star_inputs(replace(cfg, expected_total=2))
        self.assertEqual(validate_star_inputs(replace(cfg, expected_total=3)).total_particles, 3)


class MergeExecutionTests(unittest.TestCase):
    """Complete merge, metadata, provenance and publication checks on tiny inputs."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="class_pool_merge_execution_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.pd = MODULE["pd"]
        self.starfile = MODULE["starfile"]
        self.jobs = {}
        self.input_ids = {"a001", "a002", "a003", "b001", "b002"}
        for name, rows in (
            ("jobA", [["a001", 1, 0.123456789], ["a002", 2, -1.25], ["a003", 1, 2.5]]),
            ("jobB", [["b001", 2, 3.125], ["b002", 1, -0.000123456]]),
        ):
            folder = self.root / name
            folder.mkdir()
            output = folder / "run_it025_data.star"
            particles = self.pd.DataFrame(rows, columns=[
                "rlnImageName", "rlnClassNumber", "rlnOriginXAngst",
            ])
            self.starfile.write({"particles": particles}, output, float_format="%.12g")
            write_synthetic_model_star(folder / "run_it025_model.star", {1: 0.5, 2: 0.5, 3: 0.0})
            previous = folder / "previous.star"
            self.starfile.write({"particles": particles.iloc[::-1].drop(columns="rlnClassNumber")}, previous, float_format="%.12g")
            self.jobs[name] = MODULE["JobConfig"](
                name=name, star=output, k=3, classes={"Keep": [1], "Reject": [2], "Empty": [3]},
                input_pool=f"Input{name}", input_star=previous, expected_particles=len(rows),
            )
        self.cfg = MODULE["MergeConfig"](
            round_label="Round2", round_prefix="round2", project_dir=self.root,
            output_dir=self.root / "result", pools=("Keep", "Reject", "Archive", "Empty"),
            jobs=self.jobs, id_column="rlnImageName", expected_total=5,
        )

    def read_blocks(self, path):
        return self.starfile.read(path, always_dict=True, parse_as_string=["rlnImageName"])

    def rewrite(self, path, transform):
        blocks = self.read_blocks(path)
        transform(blocks)
        self.starfile.write(blocks, path, overwrite=True, float_format="%.12g")

    def merge(self, cfg=None, **kwargs):
        cfg = self.cfg if cfg is None else cfg
        return MODULE["execute_merge"](cfg, validate_star_inputs(cfg), **kwargs)

    def pool_path(self, pool, count, cfg=None):
        cfg = self.cfg if cfg is None else cfg
        return cfg.output_dir / f"particles_{cfg.round_label}_{pool}_{count}.star"

    def add_frozen(self, cfg, pool, ids, mode):
        path = self.root / f"frozen_{pool}.star"
        particles = self.pd.DataFrame(
            [[identity, 1, 7.125] for identity in ids],
            columns=["rlnImageName", "rlnClassNumber", "rlnOriginXAngst"],
        )
        self.starfile.write({"particles": particles}, path)
        # A custom comment proves copy_only preserves bytes beyond table values.
        path.write_text("# preserved frozen source comment\n" + path.read_text(), encoding="utf-8")
        frozen = MODULE["FrozenPoolConfig"](
            pool=pool, star=path, mode=mode, expected_particles=len(ids), source_label="CARRIED_ROUND1",
        )
        return replace(cfg, frozen={**cfg.frozen, pool: frozen}, expected_total=cfg.expected_total + len(ids))

    def add_optics(self, cfg=None):
        cfg = self.cfg if cfg is None else cfg
        optics = self.pd.DataFrame({
            "rlnOpticsGroup": [1, 2], "rlnOpticsGroupName": ["optics1", "optics2"],
            "rlnImagePixelSize": [1.25, 2.5], "rlnVoltage": [300., 300.],
        })
        paths = [path for job in cfg.jobs.values() for path in (job.star, job.input_star) if path]
        paths.extend(pool.star for pool in cfg.frozen.values())
        for path in paths:
            def transform(blocks):
                blocks["optics"] = optics.copy()
                blocks["particles"]["rlnOpticsGroup"] = 1
            self.rewrite(path, transform)
        return cfg

    def snapshot_inputs(self):
        return {path: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in self.root.rglob("*") if path.is_file()}

    def assert_no_partial_output(self, cfg=None):
        cfg = self.cfg if cfg is None else cfg
        self.assertFalse(cfg.output_dir.exists())
        self.assertFalse([path for path in cfg.output_dir.parent.iterdir()
                          if path.name.startswith(".") and path.is_dir()])

    def test_multi_job_merge_preserves_rows_metadata_counts_and_empty_pools(self):
        before = self.snapshot_inputs()
        output = self.merge()
        self.assertEqual(Path(output), self.cfg.output_dir)
        expected = {"Keep": ["a001", "a003", "b002"], "Reject": ["a002", "b001"], "Archive": [], "Empty": []}
        all_ids = []
        for pool, ids in expected.items():
            path = self.pool_path(pool, len(ids))
            table = self.read_blocks(path)["particles"]
            if not ids:
                self.assertEqual(len(table), 0)
                # The installed starfile reader drops column labels on empty
                # tables, so independently verify their written STAR headers.
                for column in ("rlnImageName", "rlnClassNumber", "rlnOriginXAngst"):
                    self.assertIn(f"_{column}", path.read_text())
                continue
            self.assertEqual(table["rlnImageName"].tolist(), ids)
            self.assertEqual(list(table.columns), ["rlnImageName", "rlnClassNumber", "rlnOriginXAngst"])
            all_ids.extend(table["rlnImageName"])
        self.assertEqual(set(all_ids), self.input_ids)
        self.assertEqual(len(all_ids), len(set(all_ids)))
        keep = self.read_blocks(self.pool_path("Keep", 3))["particles"]
        self.assertAlmostEqual(keep.iloc[0]["rlnOriginXAngst"], 0.123456789, places=8)
        for filename in ("manifest.json", "provenance.json", "provenance.csv", "summary.txt", "config.toml"):
            self.assertTrue((self.cfg.output_dir / filename).is_file(), filename)
        manifest = json.loads((self.cfg.output_dir / "manifest.json").read_text())
        self.assertEqual(manifest["counts"], {"active": 5, "frozen": 0, "total": 5, "unique": 5,
                                            "pools": {"Keep": 3, "Reject": 2, "Archive": 0, "Empty": 0}})
        self.assertIs(manifest["validation"]["input_output_conserved"], True)
        self.assertIs(manifest["validation"]["metadata_verified"], True)
        self.assertIs(manifest["validation"]["all_previous_inputs_verified"], True)
        self.assertEqual(set(manifest["validation"]["previous_inputs_verified"]), {"jobA", "jobB"})
        self.assertTrue(manifest["software"])
        for entry in manifest["inputs"]:
            self.assertEqual(entry["sha256"], hashlib.sha256(Path(entry["path"]).read_bytes()).hexdigest())
        for entry in manifest["outputs"].values():
            self.assertEqual(entry["sha256"], hashlib.sha256((self.cfg.output_dir / entry["file"]).read_bytes()).hexdigest())
        self.assertEqual(before, {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in before})
        effective = load_config(self.cfg.output_dir / "config.toml")
        self.assertEqual(effective.pools, self.cfg.pools)
        self.assertEqual(effective.jobs["jobA"].classes, {pool: self.jobs["jobA"].classes.get(pool, []) for pool in self.cfg.pools})

    def test_frozen_copy_only_is_byte_identical_and_append_preserves_source_order(self):
        cfg = self.add_frozen(self.cfg, "Archive", ["old001", "old002"], "copy_only")
        cfg = self.add_frozen(cfg, "Keep", ["old003"], "append")
        source_bytes = cfg.frozen["Archive"].star.read_bytes()
        self.merge(cfg)
        self.assertEqual(self.pool_path("Archive", 2, cfg).read_bytes(), source_bytes)
        keep = self.read_blocks(self.pool_path("Keep", 4, cfg))["particles"]
        self.assertEqual(keep["rlnImageName"].tolist(), ["old003", "a001", "a003", "b002"])
        provenance = json.loads((cfg.output_dir / "provenance.json").read_text())
        records = {tuple(row["identity"]): row for row in provenance["particles"]}
        self.assertEqual(records[("old003",)]["history"][-1]["source_type"], "frozen")
        self.assertEqual(records[("old003",)]["pool"], "Keep")

    def test_particle_column_mismatch_rejected_without_output(self):
        self.rewrite(self.jobs["jobB"].star, lambda blocks: blocks["particles"].__setitem__("extraColumn", 42))
        with self.assertRaisesRegex(ValueError, "(?i)column"):
            self.merge()
        self.assert_no_partial_output()

    def test_duplicate_star_blocks_cannot_hide_particles_or_metadata(self):
        self.add_optics()
        for job in self.jobs.values():
            for path in (job.star, job.input_star):
                self.rewrite(path, lambda blocks: blocks.__setitem__("notes", {"label": "retained"}))
        validate_star_inputs(self.cfg)
        source = self.jobs["jobA"].star
        original = source.read_text(encoding="utf-8")
        blocks = self.read_blocks(source)
        prefix = self.root / "hidden_block.star"
        for block_name in ("particles", "optics", "notes"):
            with self.subTest(block=block_name):
                hidden = blocks[block_name].copy()
                if block_name == "particles":
                    hidden.loc[0, "rlnImageName"] = "silently_discarded_particle"
                elif block_name == "optics":
                    hidden.loc[0, "rlnImagePixelSize"] = 999.0
                else:
                    hidden["label"] = "silently_discarded_metadata"
                self.starfile.write({block_name: hidden}, prefix, overwrite=True)
                # The last block is valid: a reader that overwrites repeated
                # names would conceal the conflicting first block completely.
                source.write_text(prefix.read_text(encoding="utf-8") + "\n" + original, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, f"duplicate STAR block 'data_{block_name}'"):
                    validate_star_inputs(self.cfg)
                self.assert_no_partial_output()
        source.write_text(original, encoding="utf-8")

    def test_block_header_checks_ignore_quoted_values_and_multiline_metadata(self):
        path = self.root / "quoted_blocks.star"
        path.write_text(
            'data_particles\n\nloop_\n_rlnImageName #1\n_rlnClassNumber #2\n'
            '"data_particles" 1\n"data_optics" 2\n', encoding="utf-8",
        )
        _, particles = MODULE["read_data_star"](path)
        self.assertEqual(particles["rlnImageName"].tolist(), ["data_particles", "data_optics"])
        path.write_text(
            'data_notes\n\n_label\n;metadata text\ndata_notes\ndata_particles\n;\n\n'
            + path.read_text(encoding="utf-8"), encoding="utf-8",
        )
        # Metadata text must not be mistaken for a new block by the safety
        # preflight; support for reading that text belongs to the STAR parser.
        MODULE["check_particle_token_syntax"](path)

    def test_particle_numeric_coercion_cannot_silently_round_large_integers(self):
        self.rewrite(self.jobs["jobA"].star, lambda blocks: blocks["particles"].__setitem__("customValue", 9007199254740993))
        self.rewrite(self.jobs["jobB"].star, lambda blocks: blocks["particles"].__setitem__("customValue", 1.5))
        with self.assertRaisesRegex(ValueError, "(?i)numeric|precision|loss|particle|dtype"):
            self.merge()
        self.assert_no_partial_output()

    def test_no_optics_is_explicitly_validated(self):
        report = validate_star_inputs(self.cfg)
        self.assertIs(report.has_optics_block, False)
        self.assertEqual(report.optics_groups, ())
        self.assertEqual(report.jobs["jobA"].optics_groups, ())

    def test_optics_roundtrip_and_harmless_ordering_and_float_differences(self):
        self.add_optics()
        def reorder(blocks):
            optics = blocks["optics"].iloc[::-1, ::-1].copy()
            optics["rlnImagePixelSize"] += 1e-8
            blocks["optics"] = optics
        self.rewrite(self.jobs["jobB"].star, reorder)
        self.merge()
        blocks = self.read_blocks(self.pool_path("Keep", 3))
        self.assertEqual(set(blocks["optics"]["rlnOpticsGroup"]), {1, 2})
        self.assertEqual(set(blocks["particles"]["rlnOpticsGroup"]), {1})

    def test_incompatible_optics_parameters_rejected(self):
        self.add_optics()
        self.rewrite(self.jobs["jobB"].star, lambda blocks: blocks["optics"].__setitem__("rlnImagePixelSize", [9., 2.5]))
        with self.assertRaisesRegex(ValueError, "(?i)optics"):
            validate_star_inputs(self.cfg)

    def test_compatible_distinct_optics_groups_are_unioned_without_renumbering(self):
        self.add_optics()
        for name, group in (("jobA", 1), ("jobB", 2)):
            for path in (self.jobs[name].star, self.jobs[name].input_star):
                def select_group(blocks):
                    blocks["optics"] = blocks["optics"].loc[blocks["optics"]["rlnOpticsGroup"] == group]
                    blocks["particles"]["rlnOpticsGroup"] = group
                self.rewrite(path, select_group)
        self.merge()
        blocks = self.read_blocks(self.pool_path("Keep", 3))
        self.assertEqual(set(blocks["optics"]["rlnOpticsGroup"]), {1, 2})
        particles = blocks["particles"].set_index("rlnImageName")
        self.assertEqual(particles.loc["a001", "rlnOpticsGroup"], 1)
        self.assertEqual(particles.loc["b002", "rlnOpticsGroup"], 2)

    def test_undefined_and_invalid_optics_group_references_rejected(self):
        self.add_optics()
        for value in (3, 0, -1, 1.5, "?"):
            with self.subTest(value=value):
                self.rewrite(self.jobs["jobA"].star, lambda blocks: blocks["particles"].__setitem__("rlnOpticsGroup", value))
                with self.assertRaisesRegex(ValueError, "(?i)optics"):
                    validate_star_inputs(self.cfg)

    def test_duplicate_optics_definitions_rejected_even_if_identical(self):
        self.add_optics()
        self.rewrite(self.jobs["jobA"].star, lambda blocks: blocks.__setitem__("optics", self.pd.concat([blocks["optics"], blocks["optics"].iloc[:1]], ignore_index=True)))
        with self.assertRaisesRegex(ValueError, "(?i)duplicate.*optics|optics.*duplicate"):
            validate_star_inputs(self.cfg)

    def test_mixed_optics_presence_and_missing_reference_columns_rejected(self):
        self.add_optics()
        self.rewrite(self.jobs["jobA"].input_star, lambda blocks: blocks.pop("optics"))
        with self.assertRaisesRegex(ValueError, "(?i)optics"):
            validate_star_inputs(self.cfg)
        self.add_optics()
        self.rewrite(self.jobs["jobA"].star, lambda blocks: blocks.__setitem__("particles", blocks["particles"].drop(columns="rlnOpticsGroup")))
        with self.assertRaisesRegex(ValueError, "(?i)optics"):
            validate_star_inputs(self.cfg)

    def test_previous_input_and_frozen_optics_are_checked(self):
        cfg = self.add_frozen(self.cfg, "Archive", ["old001"], "copy_only")
        self.add_optics(cfg)
        for path in (self.jobs["jobA"].input_star, cfg.frozen["Archive"].star):
            with self.subTest(path=path):
                original = path.read_bytes()
                self.rewrite(path, lambda blocks: blocks["optics"].__setitem__("rlnVoltage", 200.))
                with self.assertRaisesRegex(ValueError, "(?i)optics"):
                    validate_star_inputs(cfg)
                path.write_bytes(original)

    def test_other_metadata_is_preserved_and_conflicts_rejected(self):
        for job in self.jobs.values():
            for path in (job.star, job.input_star):
                self.rewrite(path, lambda blocks: blocks.update({
                    "general": {"label": "synthetic", "pixel_size": 1.25},
                    "extra": self.pd.DataFrame({"name": ["A", "B"], "value": [1.5, 2.5]}),
                }))
        self.rewrite(self.jobs["jobB"].star, lambda blocks: blocks.__setitem__("extra", blocks["extra"].iloc[::-1, ::-1]))
        self.merge()
        blocks = self.read_blocks(self.pool_path("Keep", 3))
        self.assertEqual(blocks["general"]["label"], "synthetic")
        self.assertEqual(set(blocks["extra"]["name"]), {"A", "B"})
        self.rewrite(self.jobs["jobB"].star, lambda blocks: blocks["general"].__setitem__("pixel_size", 5.))
        with self.assertRaisesRegex(ValueError, "(?i)metadata|general|block"):
            validate_star_inputs(self.cfg)

    def test_existing_directory_and_symlink_are_never_overwritten(self):
        self.cfg.output_dir.mkdir()
        marker = self.cfg.output_dir / "keep.txt"
        marker.write_text("untouched", encoding="utf-8")
        with self.assertRaises((OSError, ValueError)):
            self.merge()
        self.assertEqual(marker.read_text(), "untouched")
        link_cfg = replace(self.cfg, output_dir=self.root / "result_link")
        link_cfg.output_dir.symlink_to(self.cfg.output_dir, target_is_directory=True)
        with self.assertRaises((OSError, ValueError)):
            self.merge(link_cfg)
        self.assertTrue(link_cfg.output_dir.is_symlink())

    def test_failed_write_and_readback_leave_no_final_or_staging_directory(self):
        execute = MODULE["execute_merge"]
        with mock.patch.object(self.starfile, "write", side_effect=OSError("injected STAR write failure")):
            with self.assertRaises(OSError):
                self.merge()
        self.assert_no_partial_output()
        for hook, error in (("verify_written_star", ValueError("injected readback corruption")),
                            ("publish_directory", OSError("injected publish failure"))):
            with self.subTest(hook=hook):
                with mock.patch.dict(execute.__globals__, {hook: mock.Mock(side_effect=error)}):
                    with self.assertRaises(type(error)):
                        self.merge()
                self.assert_no_partial_output()
        # A clean retry succeeds after either failed attempt.
        self.merge()
        self.assertTrue(self.cfg.output_dir.is_dir())

    def test_publication_race_does_not_replace_a_new_empty_destination(self):
        execute = MODULE["execute_merge"]
        publish = execute.__globals__["publish_directory"]
        created_inode = []
        def concurrent_destination(staging, destination):
            Path(destination).mkdir()
            created_inode.append(Path(destination).stat().st_ino)
            return publish(staging, destination)
        with mock.patch.dict(execute.__globals__, {"publish_directory": concurrent_destination}):
            with self.assertRaises((OSError, ValueError)):
                self.merge()
        self.assertEqual(self.cfg.output_dir.stat().st_ino, created_inode[0])
        self.assertEqual(list(self.cfg.output_dir.iterdir()), [])
        self.assertFalse([path for path in self.root.iterdir() if path.name.startswith(".") and path.is_dir()])

    def unsupported_rename_library(self):
        def unsupported(*args):
            MODULE["ctypes"].set_errno(errno.EINVAL)
            return -1
        library = mock.Mock()
        library.renameat2 = mock.Mock(side_effect=unsupported)
        return library

    def test_unsupported_atomic_rename_publishes_complete_relative_symlink(self):
        with mock.patch.object(MODULE["ctypes"], "CDLL", return_value=self.unsupported_rename_library()):
            output = Path(self.merge())
        self.assertTrue(output.is_symlink())
        target = output.readlink()
        self.assertFalse(target.is_absolute())
        self.assertTrue(target.name.startswith(f".{output.name}.data-"))
        backing = output.parent / target
        self.assertTrue(backing.is_dir())
        self.assertEqual(list(output.parent.glob(f".{output.name}.data-*")), [backing])
        particles = self.read_blocks(self.pool_path("Keep", 3))["particles"]
        self.assertEqual(particles["rlnImageName"].tolist(), ["a001", "a003", "b002"])
        manifest = json.loads((output / "manifest.json").read_text())
        self.assertIs(manifest["validation"]["input_output_conserved"], True)
        self.assertEqual(manifest["counts"]["total"], 5)
        self.assertTrue((backing / "provenance.json").is_file())

    def test_symlink_publication_fallback_preserves_raced_destination_and_cleans_backing(self):
        execute = MODULE["execute_merge"]
        publish = execute.__globals__["publish_directory"]
        created_inode = []
        def concurrent_destination(staging, destination):
            Path(destination).mkdir()
            created_inode.append(Path(destination).stat().st_ino)
            return publish(staging, destination)
        with mock.patch.object(MODULE["ctypes"], "CDLL", return_value=self.unsupported_rename_library()):
            with mock.patch.dict(execute.__globals__, {"publish_directory": concurrent_destination}):
                with self.assertRaises((OSError, ValueError)):
                    self.merge()
        self.assertFalse(self.cfg.output_dir.is_symlink())
        self.assertEqual(self.cfg.output_dir.stat().st_ino, created_inode[0])
        self.assertEqual(list(self.cfg.output_dir.iterdir()), [])
        self.assertEqual(list(self.root.glob(f".{self.cfg.output_dir.name}.data-*")), [])

    def test_interrupt_after_symlink_publication_keeps_complete_backing(self):
        execute = MODULE["execute_merge"]
        publish = execute.__globals__["publish_directory"]
        def publish_then_interrupt(staging, destination):
            publish(staging, destination)
            raise KeyboardInterrupt("injected immediately after publication")
        with mock.patch.object(MODULE["ctypes"], "CDLL", return_value=self.unsupported_rename_library()):
            with mock.patch.dict(execute.__globals__, {"publish_directory": publish_then_interrupt}):
                with self.assertRaises(KeyboardInterrupt):
                    self.merge()
        self.assertTrue(self.cfg.output_dir.is_symlink())
        backing = self.cfg.output_dir.parent / self.cfg.output_dir.readlink()
        self.assertTrue(backing.is_dir())
        manifest = json.loads((self.cfg.output_dir / "manifest.json").read_text())
        self.assertIs(manifest["validation"]["input_output_conserved"], True)
        self.assertEqual(manifest["counts"]["total"], 5)
        all_ids = []
        for item in manifest["outputs"].values():
            path = self.cfg.output_dir / item["file"]
            self.assertEqual(item["sha256"], hashlib.sha256(path.read_bytes()).hexdigest())
            if item["count"]:
                all_ids.extend(self.read_blocks(path)["particles"]["rlnImageName"])
        self.assertEqual(set(all_ids), self.input_ids)
        self.assertEqual(len(all_ids), len(self.input_ids))
        self.assertTrue((backing / "provenance.json").is_file())

    def test_serialized_particle_corruption_is_detected_before_publication(self):
        original_write = self.starfile.write
        def corrupt(*args, **kwargs):
            original_write(*args, **kwargs)
            path = Path(args[1] if len(args) > 1 else kwargs["filename"])
            text = path.read_text(encoding="utf-8")
            path.write_text(text.replace("a001", "CORRUPTED_ID"), encoding="utf-8")
        with mock.patch.object(self.starfile, "write", side_effect=corrupt):
            with self.assertRaisesRegex(ValueError, "(?i)identity|identit|particle|round.trip"):
                self.merge()
        self.assert_no_partial_output()

    def test_default_execution_refuses_missing_previous_input_and_override_is_explicit(self):
        job = replace(self.jobs["jobA"], input_star=None)
        cfg = replace(self.cfg, jobs={**self.jobs, job.name: job})
        with self.assertRaisesRegex(ValueError, "input_star"):
            self.merge(cfg)
        self.assert_no_partial_output(cfg)
        self.merge(cfg, allow_unverified_inputs=True)
        manifest = json.loads((cfg.output_dir / "manifest.json").read_text())
        self.assertIs(manifest["execution"]["allow_unverified_inputs"], True)
        self.assertIsNone(manifest["validation"]["all_previous_inputs_verified"])
        summary = (cfg.output_dir / "summary.txt").read_text()
        self.assertIn("input_star", summary)
        self.assertIn("not", summary.lower())

    def test_native_provenance_extends_history_and_preserves_recorded_details(self):
        previous = {
            "schema_version": 1, "round": "Round1", "id_columns": ["rlnImageName"],
            "particles": [
                {"identity": [identity], "pool": "InputjobA" if identity.startswith("a") else "InputjobB",
                 "history": [{"round": "Round1", "pool": "InputjobA" if identity.startswith("a") else "InputjobB",
                              "source_type": "job", "source": "previousJob", "custom_annotation": "retain me"}]}
                for identity in sorted(self.input_ids)
            ],
        }
        path = self.root / "previous_provenance.json"
        path.write_text(json.dumps(previous), encoding="utf-8")
        # Provenance supplies the known original pool when configuration omits it.
        jobs = {**self.jobs, "jobA": replace(self.jobs["jobA"], input_pool=None)}
        cfg = replace(self.cfg, previous_round="Round1", previous_provenance=path, jobs=jobs)
        self.merge(cfg)
        provenance = json.loads((cfg.output_dir / "provenance.json").read_text())
        self.assertEqual(provenance["schema_version"], 1)
        self.assertEqual(provenance["id_columns"], ["rlnImageName"])
        records = {row["identity"][0]: row for row in provenance["particles"]}
        self.assertEqual(set(records), self.input_ids)
        for identity, row in records.items():
            self.assertEqual(row["history"][0]["custom_annotation"], "retain me")
            self.assertEqual(len(row["history"]), 2)
            event = row["history"][-1]
            self.assertEqual(event["round"], "Round2")
            self.assertEqual(event["source_type"], "job")
            self.assertEqual(event["source"], "jobA" if identity.startswith("a") else "jobB")
            self.assertEqual(event["input_pool"], "InputjobA" if identity.startswith("a") else "InputjobB")
            self.assertEqual(str(event["iteration"]).lstrip("0"), "25")
            self.assertIn(event["class_number"], (1, 2))
            self.assertEqual(row["pool"], "Keep" if event["class_number"] == 1 else "Reject")

    def test_unsupported_previous_provenance_fails_without_output(self):
        path = self.root / "unsupported.json"
        path.write_text(json.dumps({"schema_version": 99, "particles": []}), encoding="utf-8")
        cfg = replace(self.cfg, previous_provenance=path)
        with self.assertRaisesRegex(ValueError, "(?i)provenance|schema|unsupported"):
            self.merge(cfg)
        self.assert_no_partial_output(cfg)

    def test_native_provenance_shape_and_json_numbers_are_checked_during_validation(self):
        base = {
            "schema_version": 1, "round": "Round1", "id_columns": ["rlnImageName"],
            "particles": [{"identity": [identity],
                           "pool": "InputjobA" if identity.startswith("a") else "InputjobB",
                           "history": [{"round": "Round1"}]}
                          for identity in sorted(self.input_ids)],
        }
        cases = []
        for value in (None, [], {"researcher_note": 42}):
            data = json.loads(json.dumps(base))
            data["particles"][0]["legacy_fields"] = value
            cases.append((f"legacy_fields={value!r}", data))
        for value in (None, "", "   "):
            data = json.loads(json.dumps(base))
            if value is None:
                data.pop("round")
            else:
                data["round"] = value
            cases.append((f"round={value!r}", data))
        for value in (float("nan"), float("inf"), float("-inf")):
            data = json.loads(json.dumps(base))
            data["particles"][0]["unknown_prior_field"] = value
            cases.append((f"nonfinite={value!r}", data))
        path = self.root / "previous_provenance.json"
        cfg = replace(self.cfg, previous_provenance=path, previous_round=None)
        for label, data in cases:
            with self.subTest(case=label):
                path.write_text(json.dumps(data), encoding="utf-8")
                with self.assertRaises(ValueError):
                    validate_star_inputs(cfg)
                self.assert_no_partial_output(cfg)

    def test_incomplete_or_duplicate_previous_lineage_is_not_fabricated(self):
        rows = [{"identity": [identity], "pool": "InputjobA" if identity.startswith("a") else "InputjobB",
                 "history": [{"round": "Round1", "pool": "InputjobA" if identity.startswith("a") else "InputjobB"}]}
                for identity in sorted(self.input_ids)]
        path = self.root / "previous_provenance.json"
        for records in (rows[:-1], rows + [rows[0]]):
            with self.subTest(records=len(records)):
                path.write_text(json.dumps({"schema_version": 1, "round": "Round1",
                                            "id_columns": ["rlnImageName"], "particles": records}), encoding="utf-8")
                cfg = replace(self.cfg, previous_round="Round1", previous_provenance=path)
                with self.assertRaisesRegex(ValueError, "(?i)provenance|identity|identit|lineage|duplicate"):
                    self.merge(cfg)
                self.assert_no_partial_output(cfg)

    def test_historical_csv_provenance_preserves_all_columns(self):
        path = self.root / "previous.csv"
        rows = [{"rlnImageName": identity, "round1_pool": "InputjobA" if identity.startswith("a") else "InputjobB",
                 "round1_job": "oldJob", "researcher_note": f"note for {identity}"}
                for identity in sorted(self.input_ids)]
        self.pd.DataFrame(rows).to_csv(path, index=False)
        cfg = replace(self.cfg, previous_round="Round1", previous_provenance=path)
        self.merge(cfg)
        provenance = json.loads((cfg.output_dir / "provenance.json").read_text())
        for record in provenance["particles"]:
            # Imported fields must remain machine-readable even if the historical
            # schema has columns unknown to the current CLI.
            self.assertIn(f"note for {record['identity'][0]}", json.dumps(record))

    def test_dry_run_with_existing_output_is_read_only(self):
        self.cfg.output_dir.mkdir()
        (self.cfg.output_dir / "keep.txt").write_text("untouched", encoding="utf-8")
        config = self.root / "dry_run.toml"
        config.write_text(textwrap.dedent(f"""
            round = "Round2"
            project_dir = "{self.root}"
            output_dir = "result"
            pools = ["Keep", "Reject", "Empty"]
            id_column = "rlnImageName"
            [jobs.jobA]
            star = "jobA/run_it025_data.star"
            input_star = "jobA/previous.star"
            k = 3
            [jobs.jobA.classes]
            Keep = [1]
            Reject = [2]
            Empty = [3]
        """), encoding="utf-8")
        before = self.snapshot_inputs()
        directories_before = sorted(str(path.relative_to(self.root)) for path in self.root.rglob("*") if path.is_dir())
        result = subprocess.run([sys.executable, str(CLI), "--config", str(config), "--dry-run"], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(before, self.snapshot_inputs())
        self.assertEqual(directories_before, sorted(str(path.relative_to(self.root)) for path in self.root.rglob("*") if path.is_dir()))


if __name__ == "__main__":
    unittest.main()
