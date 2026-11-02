"""Synthetic integration checks; external refinement commands are stubbed."""

import hashlib
import json
import os
from pathlib import Path
import runpy
import struct
import subprocess
import sys
import tempfile
import unittest


CLI = Path(__file__).resolve().parents[1] / "m" / "m_progressive"
MODULE = runpy.run_path(str(CLI))


def write_mrc(path, values):
    header = bytearray(1024)
    struct.pack_into("<4i", header, 0, 2, 2, 2, 2)
    struct.pack_into("<3i", header, 28, 2, 2, 2)
    struct.pack_into("<3f", header, 40, 4, 4, 4)
    struct.pack_into("<3i", header, 64, 1, 2, 3)
    path.write_bytes(header + struct.pack("<8f", *values))


class ProgressiveTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="m_progressive_test_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.inputs = self.root / "inputs with spaces"
        self.job = self.inputs / "job"
        self.archive = self.inputs / "snapshot"
        self.processing = self.inputs / "processing"
        self.data = self.inputs / "raw"
        for folder in (self.job, self.archive, self.processing / "tiltstack" / "tomo001", self.data):
            folder.mkdir(parents=True)
        self.settings = self.inputs / "warp.settings"
        self.settings.write_text('''<Settings>
<Param Name="ProcessingFolder" Value="processing" />
<Param Name="DataFolder" Value="raw" />
<Param Name="UseRelativePaths" Value="True" />
</Settings>''')
        (self.archive / "tomo001.xml").write_text("<Synthetic />")
        (self.job / "run_data.star").write_text('''data_particles
loop_
_rlnImageName #1
_rlnTomoName #2
1@particles.mrcs tomo001.tomostar
''')
        (self.job / "RELION_JOB_EXIT_SUCCESS").touch()
        (self.job / "note.txt").write_text("--sym C1 --particle_diameter 100")
        for name in MODULE["REQUIRED_RELION"][1:]:
            write_mrc(self.job / name, [0.] * 8)
        self.mask = self.inputs / "binary.mrc"
        write_mrc(self.mask, [1.] + [0.] * 7)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        stub = '''import json, os, pathlib, sys
args = sys.argv[1:]
with open(os.environ['STUB_LOG'], 'a') as log:
    log.write(json.dumps([pathlib.Path(sys.argv[0]).name] + args) + '\\n')
if args and args[0] == 'create_population':
    directory = pathlib.Path(args[args.index('--directory') + 1])
    directory.mkdir()
    (directory / (args[args.index('--name') + 1] + '.population')).write_text('stub')
if pathlib.Path(sys.argv[0]).name == 'MCore' and os.environ.get('STUB_FAIL'):
    sys.exit(7)
'''
        for name in ("MTools", "MCore"):
            executable = self.bin / name
            executable.write_text("#!" + sys.executable + "\n" + stub)
            executable.chmod(0o755)

    def command(self, out, *extra):
        return [sys.executable, str(CLI), str(self.job), "--settings", str(self.settings),
                "--pre-m-xml", str(self.archive), "--mask", str(self.mask),
                "--out", str(out), *extra]

    def prepare(self, name="run"):
        out = self.root / name
        result = subprocess.run(self.command(out), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return out

    def run_script(self, out, **overrides):
        env = {k: v for k, v in os.environ.items()
               if not k.startswith("SLURM") and k not in
               ("GPU_ID", "START_STEP", "END_STEP", "CPU_MEMORY", "CTF_BATCH", "STUB_FAIL")}
        env.update(PATH=str(self.bin) + os.pathsep + os.environ["PATH"],
                   STUB_LOG=str(out / "stub.log"))
        env.update(overrides)
        return subprocess.run(["bash", str(out / "run_m_progressive.sh")],
                              env=env, capture_output=True, text=True)

    def test_preparation_and_gpu_selection_without_slurm(self):
        before = {p: hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in self.inputs.rglob("*") if p.is_file()}
        for gpu in ("0", "1"):
            out = self.prepare("run ' gpu" + gpu)
            summary = json.loads((out / "input_summary.json").read_text())
            self.assertEqual(summary["star_counted_rows"], 1)
            self.assertEqual(summary["mask_qc"]["foreground_voxels"], 1)
            copied = out / "warp_processing/tomo001.xml"
            self.assertFalse(copied.is_symlink())
            self.assertEqual(copied.read_bytes(), (self.archive / copied.name).read_bytes())
            self.assertEqual((out / "warp_processing/tiltstack").resolve(),
                             self.processing / "tiltstack")
            settings = (out / "warp.settings").read_text()
            self.assertEqual(MODULE["get_setting"](settings, "ProcessingFolder"),
                             str(out / "warp_processing"))
            self.assertEqual(MODULE["get_setting"](settings, "UseRelativePaths"), "False")
            result = self.run_script(out, GPU_ID=gpu, CPU_MEMORY="1", CTF_BATCH="4")
            self.assertEqual(result.returncode, 0, result.stderr)
            calls = [json.loads(line) for line in (out / "stub.log").read_text().splitlines()]
            self.assertEqual([call[1] for call in calls[:3]],
                             ["create_population", "create_source", "create_species"])
            cores = [call for call in calls if call[0] == "MCore"]
            self.assertEqual(len(cores), 6)
            for call in cores:
                self.assertEqual(call[call.index("--devicelist") + 1], gpu)
                self.assertIn("--cpu_memory", call)
                if "--ctf_batch" in call:
                    self.assertEqual(call[call.index("--ctf_batch") + 1], "4")
            self.assertEqual(len(list(out.glob(".completed_step_*"))), 9)
            self.assertNotEqual(self.run_script(out, GPU_ID=gpu).returncode, 0)
        self.assertEqual(before, {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before})

    def test_resume_and_xml_integrity(self):
        out = self.prepare()
        copied = out / "warp_processing/tomo001.xml"
        copied.write_text("<Changed />")
        self.assertIn("Working XML modified", self.run_script(out).stderr)
        self.assertFalse((out / "stub.log").exists())
        copied.write_bytes((self.archive / copied.name).read_bytes())
        self.assertIn("M population missing", self.run_script(out, START_STEP="5").stderr)
        result = self.run_script(out, END_STEP="4")
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.run_script(out, START_STEP="5")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(len(list(out.glob(".completed_step_*"))), 9)

    def test_failed_stage_is_not_marked_complete(self):
        out = self.prepare()
        result = self.run_script(out, STUB_FAIL="1")
        self.assertEqual(result.returncode, 7)
        self.assertTrue((out / ".completed_step_3").exists())
        self.assertFalse((out / ".completed_step_4").exists())

    def test_invalid_runtime_values(self):
        out = self.prepare()
        for overrides in ({"GPU_ID": "bad"}, {"START_STEP": "0"},
                          {"START_STEP": "9", "END_STEP": "4"},
                          {"CPU_MEMORY": "2"}, {"CTF_BATCH": "0"}):
            self.assertNotEqual(self.run_script(out, **overrides).returncode, 0)
        self.assertFalse((out / "stub.log").exists())

    def test_validation_and_no_overwrite(self):
        out = self.prepare()
        self.assertNotEqual(subprocess.run(self.command(out), capture_output=True).returncode, 0)
        for extra in (["--sym", "D7"], ["--image-warp", "3x0"],
                      ["--mask-threshold", "nan"], ["--ctf-batch", "0"]):
            failed_out = self.root / "invalid"
            result = subprocess.run(self.command(failed_out, *extra), capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(failed_out.exists())
        write_mrc(self.mask, [0.5] + [0.] * 7)
        result = subprocess.run(self.command(self.root / "bad_mask"), capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / "bad_mask").exists())
        self.assertFalse(list(self.root.glob(".m_progressive_*")))


if __name__ == "__main__":
    unittest.main()
