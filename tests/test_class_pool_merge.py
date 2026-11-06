"""Synthetic tests for class_pool_merge config, STAR validation, and recovery."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import textwrap
from typing import Any
import unittest


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

        # Invoking without --dry-run must fail until write phase is implemented.
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
        self.assertIn("not yet implemented", no_dry_res.stderr)
        self.assertFalse((self.root / "Round13_split").exists())

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
        self.assertIsNone(report.frozen)
        self.assertIsNone(report.frozen_particles)
        self.assertIsNone(report.unique_particle_ids)
        summary = format_validation_summary(cfg, report)
        self.assertIn("active counts + configured frozen expectations", summary)
        self.assertIn("Frozen STAR files have not been validated.", summary)
        self.assertNotIn("Validated pool totals:", summary)
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
    """Unfinished Phase 2B-1B checks must not appear as verified results."""

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

    def test_identity_optics_and_previous_inputs_remain_unverified(self):
        # Duplicate identities and a mismatched previous input are outside the
        # recovered class/count checks; neither may be reported as verified.
        data_path = self.cfg.jobs["job001"].star
        data_path.write_text(
            data_path.read_text(encoding="utf-8").replace("part_00002", "part_00001"),
            encoding="utf-8",
        )
        input_path = self.root / "previous.star"
        write_synthetic_data_star(input_path, [1], image_prefix="different/particle")
        job = replace(self.cfg.jobs["job001"], input_star=input_path)
        for identity in ("auto", "rlnImageName", ("rlnTomoName", "rlnTomoParticleId")):
            with self.subTest(identity=identity):
                cfg = replace(self.cfg, id_column=identity, jobs={job.name: job})
                report = validate_star_inputs(cfg)
                self.assertEqual(report.active_particles, 2)
                self.assertEqual(report.total_particles, 2)
                self.assertEqual(report.frozen, {})
                self.assertEqual(report.frozen_particles, 0)
                self.assertIsNone(report.id_columns)
                self.assertIsNone(report.unique_particle_ids)
                self.assertIsNone(report.optics_groups)
                self.assertIsNone(report.has_optics_block)
                self.assertIsNone(report.all_previous_inputs_verified)
                self.assertIsNone(report.jobs[job.name].input_star_verified)
                self.assertIsNone(report.jobs[job.name].optics_groups)
                summary = format_validation_summary(cfg, report)
                self.assertIn("Validated pool totals:", summary)
                self.assertIn("particle identity/uniqueness", summary)
                self.assertIn("optics compatibility", summary)
                self.assertIn("previous inputs", summary)

    def test_frozen_expectations_are_separate_from_unmeasured_counts(self):
        frozen_path = self.root / "frozen.star"
        write_synthetic_data_star(frozen_path, [1, 1, 2, 2])
        for expected in (10, None):
            with self.subTest(expected=expected):
                frozen = MODULE["FrozenPoolConfig"](
                    pool="Junk", star=frozen_path, mode="append", expected_particles=expected,
                )
                cfg = replace(self.cfg, frozen={"Junk": frozen})
                report = validate_star_inputs(cfg)
                self.assertEqual(report.active_particles, 2)
                self.assertEqual(report.frozen_expected_particles, expected or 0)
                self.assertEqual(report.total_particles, 2 + (expected or 0))
                self.assertEqual(report.pool_totals, {"Good": 1, "Junk": 1 + (expected or 0)})
                self.assertIsNone(report.frozen)
                self.assertIsNone(report.frozen_particles)
                self.assertIsNone(report.unique_particle_ids)
                summary = format_validation_summary(cfg, report)
                self.assertIn("active counts + configured frozen expectations", summary)
                self.assertIn("Frozen STAR files have not been validated.", summary)
                self.assertNotIn("Validated pool totals:", summary)
                if expected is None:
                    self.assertIn("without expected_particles are excluded from totals", summary)


if __name__ == "__main__":
    unittest.main()
