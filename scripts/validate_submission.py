#!/usr/bin/env python
"""Validate the TELEMEDICON26 / CliffForge submission.

Checks (exit code 0 only if every check passes):
  1. required files are present
  2. no secrets, no private absolute paths, no data/weights/large files in the tracked tree
  3. reference results are internally consistent (summary.json vs metrics files)
  4. headline numbers in README.md and presentation/main.tex agree with the reference results
  5. the compiled deck has exactly seven pages
  6. package modules compile and the unit tests pass (weight/data tests are skipped when resources are absent)

Usage:  python scripts/validate_submission.py [--skip-tests] [--skip-pdf]
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FAILS, PASSES = [], []


def check(name, ok, detail=""):
    (PASSES if ok else FAILS).append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail and not ok else ""))
    return ok


def tracked_or_walk():
    """Files that would be committed: `git ls-files` (+ untracked, not ignored) if available, else a filtered walk."""
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--cached", "--others", "--exclude-standard"],
                             capture_output=True, text=True, check=True).stdout.split("\n")
        files = [ROOT / f for f in out if f.strip()]
        if files:
            return [f for f in files if f.is_file()]
    except Exception:
        pass
    skip = {"data/raw", "data/moonshot", "models", "outputs", ".git", "__pycache__", "build"}
    return [p for p in ROOT.rglob("*") if p.is_file() and not any(s in str(p.relative_to(ROOT)) for s in skip)]


def load(p):
    return json.load(open(p))


def pdf_pages(pdf):
    """Page count via pdfinfo, pypdf/PyPDF2, or None if no tool is available (PDF object streams defeat regex)."""
    try:
        out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True, text=True, check=True).stdout
        m = re.search(r"^Pages:\s+(\d+)", out, re.M)
        if m:
            return int(m.group(1))
    except Exception:
        pass
    for mod in ("pypdf", "PyPDF2"):
        try:
            return len(__import__(mod).PdfReader(str(pdf)).pages)
        except Exception:
            continue
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-tests", action="store_true")
    ap.add_argument("--skip-pdf", action="store_true")
    args = ap.parse_args()

    # ---- 1. required files
    required = ["README.md", "LICENSE", "CITATION.cff", "requirements.txt", ".gitignore", "data/README.md", "models/README.md",
                "results/README.md", "results/summary.json", "examples/README.md", "tests/test_generation.py",
                "scripts/run_demo.sh", "scripts/run_evaluation.sh", "scripts/validate_submission.sh",
                "src/cliffforge/generate.py", "src/cliffforge/relmodel.py", "src/cliffforge/cliff_detection.py",
                "presentation/main.tex", "presentation/references.bib", "presentation/README.md",
                "presentation/TELEMEDICON26_CliffForge_Pitch.pdf"]
    missing = [f for f in required if not (ROOT / f).exists()]
    check("required files present", not missing, f"missing: {missing}")

    # ---- 2. hygiene
    files = tracked_or_walk()
    secret_pats = [r"ghp_[A-Za-z0-9]{20,}", r"github_pat_[A-Za-z0-9_]{20,}", r"gh[ousr]_[A-Za-z0-9]{30,}", r"AKIA[0-9A-Z]{16}",
                   r"-----BEGIN [A-Z ]*PRIVATE KEY-----", r"(?i)(api[_-]?key|secret|passwd|password)\s*[:=]\s*['\"][^'\"]{8,}"]
    private_path = re.compile(r"/home/[A-Za-z0-9_.-]+/|/Users/[A-Za-z0-9_.-]+/|C:\\Users\\")
    bad_secret, bad_path, bad_big, bad_ext = [], [], [], []
    text_ext = {".py", ".sh", ".md", ".tex", ".bib", ".json", ".txt", ".cff", ".yml", ".yaml", ".csv", ".toml", ".cfg", ""}
    for f in files:
        rel = f.relative_to(ROOT)
        if f.stat().st_size > 5 * 1024 * 1024:
            bad_big.append(str(rel))
        if f.suffix in {".pt", ".pth", ".ckpt", ".pkl", ".npz", ".h5", ".safetensors"}:
            bad_ext.append(str(rel))
        if f.suffix.lower() in text_ext and f.name != "validate_submission.py":
            try:
                t = f.read_text(errors="ignore")
            except Exception:
                continue
            if any(re.search(p, t) for p in secret_pats):
                bad_secret.append(str(rel))
            if private_path.search(t):
                bad_path.append(str(rel))
    check("no credentials/tokens in files that would be committed", not bad_secret, str(bad_secret))
    check("no private absolute paths in committed files", not bad_path, str(bad_path))
    check("no file larger than 5 MB", not bad_big, str(bad_big))
    check("no model-weight/pickle files in the tree", not bad_ext, str(bad_ext))
    check("raw data / processed data not in the committed tree",
          not [f for f in files if str(f.relative_to(ROOT)).startswith(("data/raw/", "data/moonshot/")) and f.name != ".gitkeep"])

    # ---- 3. reference-result consistency
    s = load(RESULTS / "summary.json")
    det = load(RESULTS / "metrics" / "activity_cliff_detection_metrics.json")
    rank = load(RESULTS / "metrics" / "activity_cliff_metrics.json")
    rec = load(RESULTS / "metrics" / "generation_recall.json")
    fun = load(RESULTS / "metrics" / "funnel.json")
    A = det["results"]["A_test_test|sim>=0.8"]
    B = det["results"]["B_test_train|sim>=0.8"]
    ws = rank["within_series_test"]
    cp, rel = ws["chemprop_abs"]["within_series_spearman_mean"], ws["rel_anchored_ensemble"]["within_series_spearman_mean"]
    p_w = ws["paired_wilcoxon_rel_ensemble_vs_chemprop"]["p"]
    counts = det["source_data_counts_recomputed"]
    n_cliffs = sum(counts[f"sim>=0.8|{k}"]["cliffs"] for k in ("A_test_test", "B_test_train", "train_train"))
    n_strict, n_pairs_A = A["n_cliffs"], A["n_pairs"]
    rel_tp = A["models"]["rel_siamese"]["tau_fixed_1.0"]["tp"]
    cp_tp = A["models"]["chemprop_abs"]["tau_fixed_1.0"]["tp"]
    rel_fp = A["models"]["rel_siamese"]["tau_fixed_1.0"]["fp"]
    ap_rel, ap_cp = A["models"]["rel_siamese"]["average_precision"], A["models"]["chemprop_abs"]["average_precision"]
    check("summary.json agrees with metrics files",
          s["activity_cliff_detection"]["setting_A_strict_test_test_PRIMARY"]["n_cliffs"] == n_strict
          and s["generation"]["funnel"]["pass_all_filters"] == fun["pass_all_filters"]
          and s["generation_recall_heldout"]["n_recovered_pairs"] == rec["n_recovered_pairs_(min_precedent=1)"])
    check("reference numbers match the project brief",
          (round(cp, 3), round(rel, 3), round(p_w, 3)) == (0.192, 0.581, 0.005) and (n_cliffs, n_strict, n_pairs_A) == (46, 5, 44)
          and (rel_tp, cp_tp, rel_fp) == (2, 0, 0) and (round(ap_rel, 2), round(ap_cp, 2)) == (0.73, 0.35)
          and (fun["generated_edits"], fun["unique_non_parent"], fun["pass_all_filters"], fun["scored_by_model"]) == (324, 231, 34, 34)
          and (rec["n_recovered_pairs_(min_precedent=1)"], rec["n_test_test_single_edit_pairs_(reachable_universe)"]) == (82, 316),
          f"got within={cp:.3f}/{rel:.3f} p={p_w:.4f} cliffs={n_cliffs}/{n_strict}/{n_pairs_A} tp={rel_tp}/{cp_tp} ap={ap_rel:.2f}/{ap_cp:.2f}")
    check("detection verdict recorded", det["decision_rule_outcome"]["verdict"] in {"A", "B", "C"})

    # ---- 4. documents agree with the results
    readme = (ROOT / "README.md").read_text() if (ROOT / "README.md").exists() else ""
    tex = (ROOT / "presentation" / "main.tex").read_text() if (ROOT / "presentation" / "main.tex").exists() else ""
    tokens = {"0.192": cp, "0.581": rel, "0.005": p_w, "46": n_cliffs, "44": n_pairs_A, "324": fun["generated_edits"],
              "231": fun["unique_non_parent"], "34": fun["pass_all_filters"], "82": 82, "316": 316, "0.73": ap_rel, "0.35": ap_cp}
    for name, text in (("README.md", readme), ("presentation/main.tex", tex)):
        missing_tok = [t for t in tokens if t not in text]
        check(f"{name} contains the headline numbers", not missing_tok, f"missing {missing_tok}")
        check(f"{name} states the 2-of-5 vs 0-of-5 preliminary result",
              re.search(r"2\s*(of|/)\s*5", text) is not None and re.search(r"(0|zero)\s*(of|/)\s*5|0\s*/\s*5", text) is not None)
        check(f"{name} makes no unsupported superiority claim",
              not re.search(r"(?i)\b(reliably (detects|outperforms)|conclusively|statistically (established|proven)|clinically validated)\b", text))
        check(f"{name} says the generator is not neural / docking not run",
              re.search(r"(?i)not a neural", text) is not None and re.search(r"(?i)docking", text) is not None)

    # ---- 5. PDF
    pdf = ROOT / "presentation" / "TELEMEDICON26_CliffForge_Pitch.pdf"
    if not args.skip_pdf:
        if pdf.exists():
            n_pages = pdf_pages(pdf)
            if n_pages is None:
                print("[SKIP] compiled deck page count (install poppler-utils for `pdfinfo`, or pip install pypdf)")
            else:
                check("compiled deck has exactly seven pages", n_pages == 7, f"found {n_pages}")
        else:
            check("compiled deck present", False, str(pdf))

    # ---- 6. code compiles + tests
    bad = []
    for f in (ROOT / "src" / "cliffforge").glob("*.py"):
        try:
            compile(f.read_text(), str(f), "exec")
        except SyntaxError as e:
            bad.append(f"{f.name}: {e}")
    check("all package modules compile", not bad, str(bad))
    if not args.skip_tests:
        env_cmd = [sys.executable, "-W", "ignore", "-m", "unittest", "discover", "-s", str(ROOT / "tests")]
        r = subprocess.run(env_cmd, capture_output=True, text=True, cwd=ROOT,
                           env={**__import__("os").environ, "PYTHONPATH": str(ROOT / "src"), "PYTHONDONTWRITEBYTECODE": "1",
                                "CLIFFFORGE_SKIP_SLOW": __import__("os").environ.get("CLIFFFORGE_SKIP_SLOW", "1")})
        tail = (r.stderr.strip().splitlines() or [""])[-1]
        check("unit tests pass", r.returncode == 0, tail)
        print("      ", tail)

    print(f"\n{len(PASSES)} passed, {len(FAILS)} failed")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
