#!/usr/bin/env python
"""Validate the TELEMEDICON26 / CliffForge submission.

Checks (exit code 0 only if every check passes):
  1. required files are present
  2. no secrets, no private absolute paths, no data/weights/large files in tracked files
  3. Moonshot and Dablander result bundles are internally consistent
  4. the eleven-slide deck reports both datasets' results and limitations accurately
  5. the compiled deck has exactly eleven pages
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
MOONSHOT = RESULTS / "moonshot_mpro_tanimoto" / "metrics"
DABLANDER = RESULTS / "dablander_mmp_5fold"
FAILS, PASSES = [], []


def check(name, ok, detail=""):
    (PASSES if ok else FAILS).append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  -- {detail}" if detail and not ok else ""))
    return ok


def tracked_files():
    """Scan committed files; large local result archives remain outside this hygiene check."""
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "ls-files", "--cached"],
                             capture_output=True, text=True, check=True).stdout.split("\n")
        files = [ROOT / f for f in out if f.strip()]
        if files:
            return [f for f in files if f.is_file()]
    except Exception:
        pass
    skip = {"data/raw", "data/moonshot", "models", "outputs", ".git", "__pycache__", "build",
            "results/dablander_mmp_5fold/raw", "results/dablander_mmp_5fold/fold_"}
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
                "presentation/TELEMEDICON26_CliffForge_Pitch.pdf",
                "results/moonshot_mpro_tanimoto/metrics/activity_cliff_detection_metrics.json",
                "results/moonshot_mpro_tanimoto/metrics/activity_cliff_metrics.json",
                "results/moonshot_mpro_tanimoto/metrics/generation_recall.json",
                "results/moonshot_mpro_tanimoto/metrics/funnel.json",
                "results/dablander_mmp_5fold/protocol.json",
                "results/dablander_mmp_5fold/performance_summary.json"]
    missing = [f for f in required if not (ROOT / f).exists()]
    check("required files present", not missing, f"missing: {missing}")
    if missing:
        return 1

    # ---- 2. hygiene
    files = tracked_files()
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
            # Public links can contain path segments such as /news/home/2026/.
            text_without_urls = re.sub(r"https?://[^\s<>]+", "", t)
            if private_path.search(text_without_urls):
                bad_path.append(str(rel))
    check("no credentials/tokens in files that would be committed", not bad_secret, str(bad_secret))
    check("no private absolute paths in committed files", not bad_path, str(bad_path))
    check("no tracked file larger than 5 MB", not bad_big, str(bad_big))
    check("no tracked model-weight/pickle files", not bad_ext, str(bad_ext))
    check("raw data / processed data not in the committed tree",
          not [f for f in files if str(f.relative_to(ROOT)).startswith(("data/raw/", "data/moonshot/")) and f.name != ".gitkeep"])

    # ---- 3. reference-result consistency
    s = load(RESULTS / "summary.json")
    det = load(MOONSHOT / "activity_cliff_detection_metrics.json")
    rank = load(MOONSHOT / "activity_cliff_metrics.json")
    rec = load(MOONSHOT / "generation_recall.json")
    fun = load(MOONSHOT / "funnel.json")
    dab_protocol = load(DABLANDER / "protocol.json")
    dab = load(DABLANDER / "performance_summary.json")
    A = det["results"]["A_test_test|sim>=0.8"]
    B = det["results"]["B_test_train|sim>=0.8"]
    ws = rank["within_series_test"]
    cp, rel = ws["chemprop_abs"]["within_series_spearman_mean"], ws["rel_anchored_ensemble"]["within_series_spearman_mean"]
    counts = det["source_data_counts_recomputed"]
    n_cliffs = sum(counts[f"sim>=0.8|{k}"]["cliffs"] for k in ("A_test_test", "B_test_train", "train_train"))
    n_strict, n_pairs_A = A["n_cliffs"], A["n_pairs"]
    rel_tp = A["models"]["rel_siamese"]["tau_fixed_1.0"]["tp"]
    cp_tp = A["models"]["chemprop_abs"]["tau_fixed_1.0"]["tp"]
    rel_fp = A["models"]["rel_siamese"]["tau_fixed_1.0"]["fp"]
    ap_rel, ap_cp = A["models"]["rel_siamese"]["average_precision"], A["models"]["chemprop_abs"]["average_precision"]
    moonshot_consistent = (
          s["dataset"]["n_molecules"] == rank["n_molecules"]
          and s["activity_cliff_detection"]["setting_A_strict_test_test_PRIMARY"]["n_cliffs"] == n_strict
          and s["activity_cliff_detection"]["setting_A_strict_test_test_PRIMARY"]["n_pairs"] == n_pairs_A
          and s["generation"]["funnel"]["pass_all_filters"] == fun["pass_all_filters"]
          and s["generation_recall_heldout"]["n_recovered_pairs"] == rec["n_recovered_pairs_(min_precedent=1)"]
          and abs(s["activity_cliff"]["within_series_test"]["chemprop_abs"]["mean_spearman"] - cp) < 1e-12
          and abs(s["activity_cliff"]["within_series_test"]["rel_anchored_ensemble"]["mean_spearman"] - rel) < 1e-12)
    check("Moonshot summary agrees with relocated metrics", moonshot_consistent)
    check("Moonshot benchmark numbers match frozen reference",
          (round(cp, 3), round(rel, 3)) == (0.192, 0.581) and (n_cliffs, n_strict, n_pairs_A) == (46, 5, 44)
          and (rel_tp, cp_tp, rel_fp) == (2, 0, 0) and (round(ap_rel, 2), round(ap_cp, 2)) == (0.73, 0.35)
          and (fun["generated_edits"], fun["unique_non_parent"], fun["pass_all_filters"], fun["scored_by_model"]) == (324, 231, 34, 34)
          and (rec["n_recovered_pairs_(min_precedent=1)"], rec["n_test_test_single_edit_pairs_(reachable_universe)"]) == (82, 316),
          f"got within={cp:.3f}/{rel:.3f} cliffs={n_cliffs}/{n_strict}/{n_pairs_A} tp={rel_tp}/{cp_tp} ap={ap_rel:.2f}/{ap_cp:.2f}")
    check("Moonshot lead-setting limitation recorded",
          B["n_cliffs"] == 14
          and B["models"]["rel_siamese"]["average_precision"] < B["models"]["chemprop_abs"]["average_precision"]
          and det["decision_rule_outcome"]["verdict"] == "B")

    dc = dab["counts"]
    dtest = dab["cliff_detection"]["test_test"]
    dpoints = dab["primary_metric_points"]
    dcliff = dtest["models"]
    check("Dablander protocol and counts agree with evaluation",
          dab["protocol"] == dab_protocol
          and (dc["source_molecules"], dc["source_mmps"]) == (1924, 12594)
          and dc["strict_binary_cliffs"] == dtest["n_cliffs"] == 106
          and dc["strict_binary_noncliffs"] == dtest["n_noncliffs"] == 2026
          and dc["strict_test_test_pairs"] == dtest["n_binary_pairs"] + dtest["n_halfAC_excluded_only_here"]
          and dab_protocol["cliff"]["positive_abs_delta_gte"] == 2.0
          and dab_protocol["cliff"]["negative_abs_delta_lte"] == 1.0)
    check("Dablander three headline metrics agree with detailed results",
          all(abs(dpoints["strict_AP"][m] - dcliff[m]["average_precision"]) < 1e-12
              for m in ("cliffforge", "chemprop"))
          and all(abs(dpoints["lead_top10_gain"][m] -
                      sum(fold["budget"] * fold["mean_measured_gain"]
                          for fold in dab["lead_secondary_by_fold"][m]) /
                      sum(fold["budget"] for fold in dab["lead_secondary_by_fold"][m])) < 1e-12
                  for m in ("cliffforge", "chemprop"))
          and (dcliff["cliffforge"]["tp"], dcliff["chemprop"]["tp"]) == (20, 0)
          and (round(dpoints["within_series_rho"]["cliffforge"], 3),
               round(dpoints["within_series_rho"]["chemprop"], 3)) == (0.509, 0.343))
    check("Dablander adjusted inference is labelled inconclusive where appropriate",
          dab["primary_paired_contrasts"]["within_series_rho__vs__chemprop"]["superiority_supported"]
          and not dab["primary_paired_contrasts"]["strict_AP__vs__chemprop"]["superiority_supported"]
          and not dab["primary_paired_contrasts"]["lead_top10_gain__vs__chemprop"]["superiority_supported"])

    # ---- 4. deck agrees with the results
    readme = (ROOT / "README.md").read_text() if (ROOT / "README.md").exists() else ""
    tex = (ROOT / "presentation" / "main.tex").read_text() if (ROOT / "presentation" / "main.tex").exists() else ""
    check("source has exactly eleven slides", len(re.findall(r"(?m)^\s*\\begin\{frame\}", tex)) == 11)
    moonshot_tokens = ["1,102", "961", "141", str(n_cliffs), str(n_pairs_A),
                       f"{ap_rel:.2f}", f"{ap_cp:.2f}", f"{rel:.3f}", f"{cp:.3f}",
                       "53/82", "47/82", "82/316"]
    missing_moonshot = [t for t in moonshot_tokens if t not in tex]
    check("deck reports Moonshot cliff, ranking and proposal metrics", not missing_moonshot,
          f"missing {missing_moonshot}")
    check("deck reports Moonshot's 2-of-5 vs 0-of-5 cliff flags",
          re.search(r"2\s*(of|/)\s*5", tex) is not None and re.search(r"0\s*(of|/)\s*5", tex) is not None)
    dab_tokens = ["1,924", "12,594", "2,132", "106", "343", "20/106", "0/106",
                  f"{dpoints['strict_AP']['cliffforge']:.3f}", f"{dpoints['strict_AP']['chemprop']:.3f}",
                  f"{dpoints['within_series_rho']['cliffforge']:.3f}", f"{dpoints['within_series_rho']['chemprop']:.3f}",
                  f"+{dpoints['lead_top10_gain']['cliffforge']:.3f}", f"+{dpoints['lead_top10_gain']['chemprop']:.3f}",
                  "2,026"]
    missing_dab = [t for t in dab_tokens if t not in tex]
    check("deck reports Dablander cliff, ranking and lead metrics", not missing_dab,
          f"missing {missing_dab}")
    check("deck labels the key interpretation limits",
          re.search(r"(?i)no measured gain for generated candidates", tex) is not None
          and re.search(r"(?i)different cliff definitions", tex) is not None
          and re.search(r"(?i)inconclusive", tex) is not None)
    for name, text in (("README.md", readme), ("presentation/main.tex", tex)):
        check(f"{name} makes no unsupported superiority claim",
              not re.search(r"(?i)\b(reliably (detects|outperforms)|conclusively|statistically (established|proven)|clinically validated)\b", text))
        check(f"{name} distinguishes prototype generation and docking",
              re.search(r"(?i)not a neural", text) is not None and re.search(r"(?i)docking", text) is not None)

    # ---- 5. PDF
    pdf = ROOT / "presentation" / "TELEMEDICON26_CliffForge_Pitch.pdf"
    if not args.skip_pdf:
        if pdf.exists():
            n_pages = pdf_pages(pdf)
            if n_pages is None:
                print("[SKIP] compiled deck page count (install poppler-utils for `pdfinfo`, or pip install pypdf)")
            else:
                check("compiled deck has exactly eleven pages", n_pages == 11, f"found {n_pages}")
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
