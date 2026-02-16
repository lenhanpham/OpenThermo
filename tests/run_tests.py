#!/usr/bin/env python3
"""OpenThermo regression test runner.

Usage:
    python run_tests.py                          # Run all tests
    python run_tests.py --generate               # Generate reference values
    python run_tests.py --test <id>              # Run a single test
    python run_tests.py --executable <path>      # Override binary path
    python run_tests.py --verbose                # Show details on failure
"""

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_REF_FILE = SCRIPT_DIR / "reference_values.json"

DEFAULT_TOLERANCES = {
    "energy_au": 1e-6,
    "energy_kJ": 0.01,
    "energy_kcal": 0.01,
    "entropy_J": 0.01,
    "heat_capacity_J": 0.01,
}

PATTERNS = {
    "point_group": re.compile(
        r"Point group:\s+(\S+)"
    ),
    "rotational_symmetry_number": re.compile(
        r"Rotational symmetry number:\s+(\d+)"
    ),
    "electronic_energy_au": re.compile(
        r"Electronic energy:\s+([-\d.Ee+]+)\s+a\.u\."
    ),
    "zpe_kJ": re.compile(
        r"Zero point energy \(ZPE\):\s+([\d.Ee+]+)\s+kJ/mol"
    ),
    "zpe_kcal": re.compile(
        r"Zero point energy \(ZPE\):\s+[\d.Ee+]+\s+kJ/mol\s+([\d.Ee+]+)\s+kcal/mol"
    ),
    "zpe_au": re.compile(
        r"Zero point energy \(ZPE\):.*?([\d.Ee+]+)\s+a\.u\."
    ),
    "thermal_correction_U_kJ": re.compile(
        r"Thermal correction to U:\s+([\d.Ee+]+)\s+kJ/mol"
    ),
    "thermal_correction_U_au": re.compile(
        r"Thermal correction to U:.*?([\d.Ee+]+)\s+a\.u\."
    ),
    "thermal_correction_H_kJ": re.compile(
        r"Thermal correction to H:\s+([\d.Ee+]+)\s+kJ/mol"
    ),
    "thermal_correction_H_au": re.compile(
        r"Thermal correction to H:.*?([\d.Ee+]+)\s+a\.u\."
    ),
    "thermal_correction_G_kJ": re.compile(
        r"Thermal correction to G:\s+([\d.Ee+]+)\s+kJ/mol"
    ),
    "thermal_correction_G_au": re.compile(
        r"Thermal correction to G:.*?([\d.Ee+]+)\s+a\.u\."
    ),
    "total_CV_J": re.compile(
        r"Total CV:\s+([\d.Ee+]+)\s+J/mol/K"
    ),
    "total_CP_J": re.compile(
        r"Total CP:\s+([\d.Ee+]+)\s+J/mol/K"
    ),
    "total_S_J": re.compile(
        r"Total S:\s+([\d.Ee+]+)\s+J/mol/K"
    ),
    "E_plus_ZPE_au": re.compile(
        r"Sum of electronic energy and ZPE.*?:\s+([-\d.Ee+]+)\s+a\.u\."
    ),
    "E_plus_U_au": re.compile(
        r"Sum of electronic energy and thermal correction to U:\s+([-\d.Ee+]+)\s+a\.u\."
    ),
    "E_plus_H_au": re.compile(
        r"Sum of electronic energy and thermal correction to H:\s+([-\d.Ee+]+)\s+a\.u\."
    ),
    "E_plus_G_au": re.compile(
        r"Sum of electronic energy and thermal correction to G:\s+([-\d.Ee+]+)\s+a\.u\."
    ),
}

BATCH_PATTERNS = {
    "weighted_U_au": re.compile(
        r"^\s+U:\s+([-\d.Ee+]+)\s+a\.u\.", re.MULTILINE
    ),
    "weighted_H_au": re.compile(
        r"^\s+H:\s+([-\d.Ee+]+)\s+a\.u\.", re.MULTILINE
    ),
    "weighted_G_au": re.compile(
        r"^\s+G:\s+([-\d.Ee+]+)\s+a\.u\.", re.MULTILINE
    ),
    "weighted_S_J": re.compile(
        r"^\s+S:\s+([\d.Ee+]+)\s+J/mol/K", re.MULTILINE
    ),
    "weighted_CV_J": re.compile(
        r"^\s+CV:\s+([\d.Ee+]+)\s+J/mol/K", re.MULTILINE
    ),
    "weighted_CP_J": re.compile(
        r"^\s+CP:\s+([\d.Ee+]+)\s+J/mol/K", re.MULTILINE
    ),
    "conf_entropy_J": re.compile(
        r"Conformation entropy:\s+([\d.Ee+]+)\s+J/mol/K"
    ),
}

STRING_FIELDS = {"point_group"}
INT_FIELDS = {"rotational_symmetry_number"}

DEFAULT_TESTS = [
    {
        "id": "gaussian_h2co_freq",
        "description": "Gaussian 16 H2CO frequency calculation",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": [],
    },
    {
        "id": "gaussian_large",
        "description": "Gaussian large molecule",
        "input_file": "gaussian.log",
        "extra_args": [],
    },
    {
        "id": "gaussian_c2h5",
        "description": "Gaussian 16 C2H5 opt+freq",
        "input_file": "gaussian_16_C2H5_optfreq.out",
        "extra_args": [],
    },
    {
        "id": "orca_optfreq",
        "description": "ORCA opt+freq calculation",
        "input_file": "ORCA-optFreq.log",
        "extra_args": [],
    },
    {
        "id": "orca_ethanol",
        "description": "ORCA ethanol opt+freq",
        "input_file": "ORCA_ethanol_optfreq.out",
        "extra_args": [],
    },
    {
        "id": "gamess_ch4",
        "description": "GAMESS CH4 methane",
        "input_file": "gamess-CH4.log",
        "extra_args": [],
    },
    {
        "id": "gamess_h2co",
        "description": "GAMESS H2CO opt+freq",
        "input_file": "GMS_H2CO_optfreq.out",
        "extra_args": [],
    },
    {
        "id": "nwchem_nh2coh",
        "description": "NWChem NH2COH opt+freq",
        "input_file": "NWChem_NH2COH_optfreq.out",
        "extra_args": [],
    },
    {
        "id": "cp2k_co2_freq",
        "description": "CP2K CO2 frequency calculation",
        "input_file": "CP2K_CO2_freq.out",
        "extra_args": [],
    },
    {
        "id": "vasp_co_freq",
        "description": "VASP CO molecule frequency",
        "input_file": "vasp_CO_freq/OUTCAR",
        "extra_args": [],
    },
    {
        "id": "gaussian_h2co_T500",
        "description": "H2CO at 500 K",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-T", "500"],
    },
    {
        "id": "gaussian_h2co_P2",
        "description": "H2CO at 2 atm",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-P", "2.0"],
    },
    {
        "id": "gaussian_h2co_truhlar",
        "description": "H2CO with Truhlar quasi-RRHO",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-lowvibmeth", "1"],
    },
    {
        "id": "gaussian_h2co_grimme",
        "description": "H2CO with Grimme quasi-RRHO",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-lowvibmeth", "2"],
    },
    {
        "id": "gaussian_h2co_minenkov",
        "description": "H2CO with Minenkov quasi-RRHO",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-lowvibmeth", "3"],
    },
    {
        "id": "gaussian_h2co_headgordon",
        "description": "H2CO with Head-Gordon quasi-RRHO (energy + entropy interpolation)",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-lowvibmeth", "4"],
    },
    {
        "id": "gaussian_h2co_headgordon_no_entropy",
        "description": "H2CO with Head-Gordon quasi-RRHO (energy only, no entropy interpolation)",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-lowvibmeth", "4", "-hg_entropy", "false"],
    },
    {
        "id": "gaussian_h2co_scaled",
        "description": "H2CO with ZPE scaling 0.96",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-sclZPE", "0.96"],
    },
    {
        "id": "gaussian_h2co_conc",
        "description": "H2CO with concentration 1.0 M",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-conc", "1.0"],
    },
    {
        "id": "nwchem_nh2coh_imagreal",
        "description": "NH2COH with imagreal threshold 20",
        "input_file": "NWChem_NH2COH_optfreq.out",
        "extra_args": ["-imagreal", "20"],
    },
    {
        "id": "vasp_co_ni111_ipmode1",
        "description": "VASP CO/Ni(111) condensed phase",
        "input_file": "vasp_CO_Ni111_freq/OUTCAR",
        "extra_args": ["-ipmode", "1"],
    },
    {
        "id": "gaussian_h2co_otm",
        "description": "H2CO .otm file generation",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-outotm", "1"],
        "check_files": [".otm"],
    },
    {
        "id": "gaussian_h2co_vibcon",
        "description": "H2CO .vibcon file generation",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-prtvib", "-1"],
        "check_files": [".vibcon"],
    },
    {
        "id": "gaussian_h2co_Tscan",
        "description": "H2CO temperature scan 200-400K",
        "input_file": "Gaussian_16_H2CO_freq.out",
        "extra_args": ["-T", "200", "400", "50"],
        "check_files": [".UHG", ".SCq"],
    },
    {
        "id": "error_nonexistent_file",
        "description": "Non-existent input file",
        "input_file": "NONEXISTENT_FILE.out",
        "extra_args": [],
        "expect_error": True,
    },
    {
        "id": "error_invalid_input",
        "description": "Invalid input file",
        "input_file": "invalid_input.txt",
        "extra_args": [],
        "expect_error": True,
    },
    {
        "id": "batch_conformers",
        "description": "Conformer ensemble (3 BIH)",
        "input_file": "conformer_ensemble.list",
        "extra_args": [],
        "test_mode": "batch",
    },
]


def find_executable():
    """Locate the OpenThermo binary in standard build locations."""
    candidates = [
        PROJECT_ROOT / "build" / "OpenThermo.exe",
        PROJECT_ROOT / "build" / "OpenThermo",
        PROJECT_ROOT / "build" / "Release" / "OpenThermo.exe",
        PROJECT_ROOT / "build" / "Debug" / "OpenThermo.exe",
        PROJECT_ROOT / "build" / "Debug" / "OpenThermo",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def parse_output(stdout_text, patterns=None):
    """Extract thermodynamic values from OpenThermo stdout."""
    if patterns is None:
        patterns = PATTERNS
    result = {}
    for key, pattern in patterns.items():
        match = pattern.search(stdout_text)
        if match:
            val = match.group(1)
            if key in STRING_FIELDS:
                result[key] = val
            elif key in INT_FIELDS:
                result[key] = int(val)
            else:
                result[key] = float(val)
        else:
            result[key] = None
    return result


def run_single_test(executable, test_def, tests_dir):
    """Run OpenThermo on one test case.

    Returns (parsed_output, raw_stdout, raw_stderr, return_code).
    """
    input_path = tests_dir / test_def["input_file"]
    expect_error = test_def.get("expect_error", False)

    if not expect_error and not input_path.exists():
        return None, "", f"Input file not found: {input_path}", -1

    cmd = [str(executable), str(input_path), "-noset"]
    cmd += test_def.get("extra_args", [])

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=120,
            cwd=str(tests_dir),
        )
        if test_def.get("test_mode") == "batch":
            parsed = parse_output(result.stdout, BATCH_PATTERNS)
        else:
            parsed = parse_output(result.stdout)
        return parsed, result.stdout, result.stderr, result.returncode
    except subprocess.TimeoutExpired:
        return None, "", "TIMEOUT (120s)", -1


def get_tolerance_for_key(key, tolerances):
    """Map a value key to its appropriate tolerance."""
    if "_au" in key:
        return tolerances.get("energy_au", 1e-6)
    if "_kJ" in key:
        return tolerances.get("energy_kJ", 0.01)
    if "_kcal" in key:
        return tolerances.get("energy_kcal", 0.01)
    if key.startswith("total_S"):
        return tolerances.get("entropy_J", 0.01)
    if key.startswith("total_C"):
        return tolerances.get("heat_capacity_J", 0.01)
    return 1e-6


def compare_values(actual, expected, tolerances):
    """Compare actual vs expected values.

    Returns list of (key, actual_val, expected_val, passed, detail).
    """
    results = []
    for key, exp_val in expected.items():
        if exp_val is None:
            continue
        act_val = actual.get(key)
        if act_val is None:
            results.append((key, None, exp_val, False, "value not found in output"))
            continue

        if key in STRING_FIELDS:
            passed = act_val == exp_val
            detail = "" if passed else f"got '{act_val}', expected '{exp_val}'"
        elif key in INT_FIELDS:
            passed = act_val == exp_val
            detail = "" if passed else f"got {act_val}, expected {exp_val}"
        else:
            tol = get_tolerance_for_key(key, tolerances)
            delta = abs(act_val - exp_val)
            passed = delta <= tol
            detail = "" if passed else f"expected {exp_val}, got {act_val} (delta={delta:.2e}, tol={tol:.0e})"

        results.append((key, act_val, exp_val, passed, detail))
    return results


def find_output_file(input_path, ext):
    """Find an output file created by OpenThermo.

    Tests run with cwd=tests_dir, so output files land next to the input.
    Returns the path if found, else None.
    """
    fpath = input_path.parent / (input_path.stem + ext)
    return fpath if fpath.exists() else None


def load_reference(ref_file):
    """Load reference_values.json."""
    with open(ref_file, "r") as f:
        return json.load(f)


def save_reference(ref_data, ref_file):
    """Write reference_values.json."""
    with open(ref_file, "w") as f:
        json.dump(ref_data, f, indent=2)
    print(f"Reference values written to {ref_file}")


def generate_references(executable, ref_file):
    """Run all tests and write/update reference_values.json."""
    tests_dir = SCRIPT_DIR

    # Load existing or create skeleton
    if ref_file.exists():
        ref_data = load_reference(ref_file)
        test_defs = ref_data.get("tests", [])
    else:
        test_defs = []
        ref_data = {}

    # If no tests defined, use defaults
    if not test_defs:
        test_defs = DEFAULT_TESTS

    # Get git commit hash
    commit_hash = "unknown"
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, cwd=str(PROJECT_ROOT)
        )
        if result.returncode == 0:
            commit_hash = result.stdout.strip()
    except FileNotFoundError:
        pass

    ref_data["_meta"] = {
        "version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generated_by": f"OpenThermo ({commit_hash})",
        "tolerances": DEFAULT_TOLERANCES,
    }

    print(f"Generating reference values using: {executable}")
    print(f"Commit: {commit_hash}")
    print()

    updated_tests = []
    for test_def in test_defs:
        test_id = test_def["id"]
        print(f"  Running {test_id}...", end=" ", flush=True)

        if test_def.get("expect_error", False):
            test_def["expected"] = {}
            print("SKIP (error test)")
            updated_tests.append(test_def)
            continue

        parsed, stdout, stderr, rc = run_single_test(executable, test_def, tests_dir)

        # Clean up generated output files
        input_path = tests_dir / test_def["input_file"]
        for ext in test_def.get("check_files", []):
            fpath = find_output_file(input_path, ext)
            if fpath:
                fpath.unlink()

        if rc != 0 or parsed is None:
            print(f"FAILED (exit code {rc})")
            if stderr:
                print(f"    stderr: {stderr[:200]}")
            test_def["expected"] = {}
        else:
            # Filter out None values
            expected = {k: v for k, v in parsed.items() if v is not None}
            test_def["expected"] = expected
            print(f"OK ({len(expected)} values captured)")

        updated_tests.append(test_def)

    ref_data["tests"] = updated_tests
    save_reference(ref_data, ref_file)


def run_tests(executable, ref_file, test_filter=None, verbose=False):
    """Run regression tests and report results."""
    tests_dir = SCRIPT_DIR

    if not ref_file.exists():
        print(f"Error: Reference file not found: {ref_file}")
        print("Run with --generate first to create reference values.")
        return 2

    ref_data = load_reference(ref_file)
    tolerances = ref_data.get("_meta", {}).get("tolerances", DEFAULT_TOLERANCES)
    test_defs = ref_data.get("tests", [])

    if not test_defs:
        print("Error: No tests defined in reference file.")
        return 2

    meta = ref_data.get("_meta", {})
    print("OpenThermo Regression Test Suite")
    print("=================================")
    print(f"Executable: {executable}")
    print(f"Reference:  {ref_file}")
    if meta.get("generated_at"):
        print(f"Generated:  {meta['generated_at']} ({meta.get('generated_by', '?')})")
    print()

    passed_count = 0
    failed_count = 0
    error_count = 0
    total_time = 0.0

    for test_def in test_defs:
        test_id = test_def["id"]
        description = test_def.get("description", "")
        expect_error = test_def.get("expect_error", False)
        check_files = test_def.get("check_files", [])

        if test_filter and test_id != test_filter:
            continue

        if not expect_error:
            expected = test_def.get("expected", {})
            if not expected and not check_files:
                print(f"[SKIP] {test_id:<25} - {description} (no reference values)")
                continue

        start = time.time()
        parsed, stdout, stderr, rc = run_single_test(executable, test_def, tests_dir)
        elapsed = time.time() - start
        total_time += elapsed

        # Compute input path for file checks
        input_path = tests_dir / test_def["input_file"]

        # Handle error tests (expect non-zero exit code)
        if expect_error:
            if rc != 0:
                passed_count += 1
                print(f"[PASS] {test_id:<25} - {description} ({elapsed:.1f}s)")
            else:
                failed_count += 1
                print(f"[FAIL] {test_id:<25} - {description} ({elapsed:.1f}s)")
                print(f"         Expected non-zero exit code, got 0")
            continue

        if rc != 0 or parsed is None:
            error_count += 1
            print(f"[ERR ] {test_id:<25} - {description} ({elapsed:.1f}s)")
            print(f"         Exit code: {rc}")
            if stderr:
                for line in stderr.strip().split("\n")[:3]:
                    print(f"         {line}")
            if verbose and stdout:
                print(f"         --- stdout (last 500 chars) ---")
                print(f"         {stdout[-500:]}")
            # Clean up any generated files even on error
            for ext in check_files:
                fpath = find_output_file(input_path, ext)
                if fpath:
                    fpath.unlink()
            continue

        # Check output file generation
        files_ok = True
        missing_files = []
        for ext in check_files:
            fpath = find_output_file(input_path, ext)
            if fpath is None:
                files_ok = False
                missing_files.append(ext)
            else:
                fpath.unlink()

        # Compare numerical values
        expected = test_def.get("expected", {})
        if expected:
            comparisons = compare_values(parsed, expected, tolerances)
            values_ok = all(c[3] for c in comparisons)
        else:
            comparisons = []
            values_ok = True

        all_passed = values_ok and files_ok

        if all_passed:
            passed_count += 1
            print(f"[PASS] {test_id:<25} - {description} ({elapsed:.1f}s)")
        else:
            failed_count += 1
            print(f"[FAIL] {test_id:<25} - {description} ({elapsed:.1f}s)")
            for ext in missing_files:
                print(f"         Output file not created: {ext}")
            for key, act, exp, ok, detail in comparisons:
                if not ok:
                    print(f"         {key}: {detail}")

            if verbose and stdout:
                print(f"         --- Full stdout ---")
                for line in stdout.strip().split("\n"):
                    print(f"         | {line}")

    print()
    print("=================================")
    total = passed_count + failed_count + error_count
    print(f"Results: {passed_count} passed, {failed_count} failed, {error_count} errors out of {total} tests ({total_time:.1f}s)")

    if failed_count > 0 or error_count > 0:
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description="OpenThermo regression test runner")
    parser.add_argument(
        "--generate", action="store_true",
        help="Generate reference values from current binary"
    )
    parser.add_argument(
        "--executable", type=str, default=None,
        help="Path to OpenThermo executable"
    )
    parser.add_argument(
        "--test", type=str, default=None,
        help="Run only the named test (by id)"
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Show detailed output on failure"
    )
    parser.add_argument(
        "--ref-file", type=str, default=None,
        help="Path to reference_values.json"
    )
    args = parser.parse_args()

    # Locate executable
    if args.executable:
        executable = Path(args.executable).resolve()
    else:
        executable = find_executable()

    if executable is None or not executable.exists():
        print(f"Error: OpenThermo executable not found.")
        print("Build first with 'make' or specify with --executable.")
        sys.exit(2)

    ref_file = Path(args.ref_file) if args.ref_file else DEFAULT_REF_FILE

    if args.generate:
        generate_references(executable, ref_file)
    else:
        rc = run_tests(executable, ref_file, test_filter=args.test, verbose=args.verbose)
        sys.exit(rc)


if __name__ == "__main__":
    main()
