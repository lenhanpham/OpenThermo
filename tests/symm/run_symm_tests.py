#!/usr/bin/env python3
"""Symmetry module test driver: runs test_symm binary over all tests/symm/ geometries.
Usage:
  python3 tests/symm/run_symm_tests.py [--delta 0.01] [--sweep] [--timeout 120]
Outputs JSON + human-readable summary. Crash-safe (subprocess per file).
"""
import argparse, json, re, subprocess, sys, time
from pathlib import Path

SYMM_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SYMM_DIR.parent.parent
BINARY = Path("/tmp/test_symm")

# Filename -> expected point group (strict). None = TBD (report only).
MANUAL = {
    "C2v.xyz": "C2v", "D2h.xyz": "D2h", "D2h-2nd.xyz": "D2h",
    "D4d.xyz": "D4d", "Ih.xyz": "Ih", "c60.xyz": "Ih",
    "C4v": "C4v", "400-atoms.xyz": None,  # unknown a priori; report only
}
# G_Cinfv file content label is "Cinfv" but code uses "Civ"; same for Dinfh->Dih.
LABEL_FIX = {"Cinfv": "Civ", "Dinfh": "Dih"}

def expected_for(fname: str):
    if fname in MANUAL:
        return MANUAL[fname]
    if fname.startswith("G_"):
        base = fname[2:]
        if base.endswith(".xyz"):
            base = base[:-4]
        return LABEL_FIX.get(base, base)
    if fname.startswith("M_"):
        core = fname[2:]
        if core.endswith(".xyz"):
            core = core[:-4]
        core = re.sub(r"_\d+$", "", core)  # strip trailing _N instance index
        return core or None
    return None

def geom_files():
    out = []
    for p in sorted(SYMM_DIR.iterdir()):
        if not p.is_file():
            continue
        if p.suffix in (".cpp", ".h") or p.name in ("run_symm_tests.py", "test_symm.cpp"):
            continue
        out.append(p)
    return out

def run_one(path: Path, delta: float, timeout: int, mode: str = "pg_determ"):
    t0 = time.time()
    try:
        r = subprocess.run([str(BINARY), str(path), "--delta", str(delta), "--mode", mode],
                           capture_output=True, text=True, timeout=timeout)
        dt = time.time() - t0
        txt = (r.stdout or "").strip()
        try:
            data = json.loads(txt.splitlines()[-1]) if txt else {}
        except Exception:
            data = {"raw": txt[-300:]}
        data["_rc"] = r.returncode
        data["_stderr_tail"] = (r.stderr or "")[-300:]
        data["_wall_s"] = round(dt, 2)
        if r.returncode not in (0, 1):
            data["_crash"] = True
        return data
    except subprocess.TimeoutExpired:
        return {"file": str(path), "_rc": -9, "_timeout": True, "_wall_s": timeout}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.01)
    ap.add_argument("--sweep", action="store_true",
                    help="sweep deltas [0.001,0.005,0.01,0.05,0.1] (400-atoms only at default)")
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--mode", type=str, default="pg_determ",
                    help="harness mode: pg_determ (single delta) or detect (full detectPG retry grid)")
    ap.add_argument("--out", type=str, default="/tmp/symm_results.json")
    args = ap.parse_args()

    deltas = [0.001, 0.005, 0.01, 0.05, 0.1] if args.sweep else [args.delta]
    files = geom_files()
    print(f"Testing {len(files)} geometries, deltas={deltas}, mode={args.mode}", flush=True)
    results = []
    for path in files:
        exp = expected_for(path.name)
        entry = {"file": path.name, "expected": exp, "runs": {}}
        for d in deltas:
            if path.name == "400-atoms.xyz" and args.sweep and d != 0.01:
                entry["runs"][str(d)] = {"skipped": "heavy file, default delta only"}
                continue
            data = run_one(path, d, args.timeout, args.mode)
            got = data.get("pg")
            if exp is None:
                verdict = "REPORT"
            elif data.get("_crash") or data.get("_timeout") or data.get("error"):
                verdict = "FAIL"
            else:
                verdict = "PASS" if got == exp else "FAIL"
            data["_verdict"] = verdict
            entry["runs"][str(d)] = data
            print(f"  {path.name:22s} d={d:<6} exp={str(exp):6s} got={str(got):6s} {verdict} "
                  f"({data.get('_wall_s', '?')}s rc={data.get('_rc')})", flush=True)
        results.append(entry)

    with open(args.out, "w") as f:
        json.dump({"deltas": deltas, "results": results}, f, indent=1)
    # summary at default delta
    d0 = str(0.01 if args.sweep else args.delta)
    n_pass = sum(1 for e in results if e["runs"].get(d0, {}).get("_verdict") == "PASS")
    n_fail = sum(1 for e in results if e["runs"].get(d0, {}).get("_verdict") == "FAIL")
    n_rep = sum(1 for e in results if e["runs"].get(d0, {}).get("_verdict") == "REPORT")
    print(f"\nSummary @ delta={d0}: {n_pass} PASS, {n_fail} FAIL, {n_rep} REPORT-only (of {len(results)})")
    print(f"Full results: {args.out}")

if __name__ == "__main__":
    main()
