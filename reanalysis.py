"""
Re-analysis of the DeNovoSeer test-set predictions (Qiu et al., Brief Bioinform 2026, bbag428).

Usage:
    python reanalysis.py --s1 path/to/Supplementary_File_S1.csv

The Gene4Denovo annotation files are downloaded from the authors' GitHub
repository on the first run and kept in ./data. Results go to ./results.
"""

import argparse
import os
import urllib.request

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = "https://media.githubusercontent.com/media/dadaxiaoxiao2333/DeNovoSeer/main/data/Gene4Denovo"
FILES = {
    "annovar": "annovar_coding_Gene4Denovo.csv",
    "clinical": "Clinical_coding_Gene4Denovo.csv",
}
EVIDENCE = ["pvs1", "ps1", "ps2", "ps3", "ps4", "pm1", "pm2", "pm3", "pm4", "pm5", "pm6",
            "pp1", "pp2", "pp3", "pp4", "ba1", "bs1", "bs2", "bs3", "bs4",
            "bp1", "bp2", "bp3", "bp4", "bp5", "bp7"]
SPLICE_DS = ["DS_AG", "DS_AL", "DS_DG", "DS_DL"]
N_BOOT = 2000
SEED = 2026


def fetch(data_dir):
    os.makedirs(data_dir, exist_ok=True)
    paths = {}
    for key, name in FILES.items():
        path = os.path.join(data_dir, name)
        if not os.path.exists(path):
            print(f"downloading {name} ...")
            urllib.request.urlretrieve(f"{REPO}/{name}", path)
        paths[key] = path
    return paths


def variant_id(df):
    return ("chr" + df["Chr"].astype(str) + "_start" + df["Start"].astype(str)
            + "_end" + df["End"].astype(str) + "_ref" + df["Ref"] + "_alt" + df["Alt"])


def load(paths):
    ann = pd.read_csv(paths["annovar"], na_values=".", low_memory=False)
    clin = pd.read_csv(paths["clinical"], dtype=str, na_values=".")
    ann["variant_index"] = variant_id(ann)
    clin["variant_index"] = variant_id(clin)
    # the two files are row-aligned in the authors' preprocessing; check anyway
    assert (ann["variant_index"].values == clin["variant_index"].values).all()
    for col in EVIDENCE:
        clin[col] = pd.to_numeric(clin[col]).fillna(0)
    ann["consequence"] = ann["ExonicFunc.refGeneWithVer"].fillna(ann["Func.refGeneWithVer"])
    ann["label_source"] = np.where(ann["ClinVar_label"].notna(), "ClinVar",
                                   np.where(ann["label"] != -1, "rule", "unlabeled"))
    return ann, clin


def boot_auc(y, x, rng, n=N_BOOT):
    y, x = np.asarray(y), np.asarray(x)
    vals = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            vals.append(roc_auc_score(y[i], x[i]))
    return np.percentile(vals, [2.5, 97.5])


def paired_diff(y, a, b, rng, n=N_BOOT):
    y, a, b = np.asarray(y), np.asarray(a), np.asarray(b)
    diffs = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            diffs.append(roc_auc_score(y[i], a[i]) - roc_auc_score(y[i], b[i]))
    return np.median(diffs), *np.percentile(diffs, [2.5, 97.5])


def score_row(group, name, y, x, rng):
    ok = x.notna()
    lo, hi = boot_auc(y[ok], x[ok], rng)
    return {"group": group, "score": name, "n": int(ok.sum()), "n_pathogenic": int(y[ok].sum()),
            "AUC": roc_auc_score(y[ok], x[ok]), "AUC_lo": lo, "AUC_hi": hi,
            "AP": average_precision_score(y[ok], x[ok])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--s1", required=True, help="Supplementary File S1 (variant-level predictions)")
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(SEED)

    ann, clin = load(fetch(args.data))
    s1 = pd.read_csv(args.s1)

    feats = ann[["variant_index", "label_source", "consequence", "Clinical_evidence_score",
                 "BayesDel_addAF_score", "AlphaMissense_score", "CADD_phred"] + SPLICE_DS]
    test = s1.merge(feats, on="variant_index", how="left", validate="1:1")
    test = test.merge(clin[["variant_index", "pm2"]], on="variant_index", how="left")
    assert test["label_source"].notna().all(), "some S1 variants not found in the annotation file"
    y = test["true_label"]

    lines = []
    def say(msg=""):
        print(msg)
        lines.append(msg)

    # 1. performance by label source
    rows = []
    groups = {"all": test, "ClinVar": test[test.label_source == "ClinVar"],
              "rule": test[test.label_source == "rule"]}
    scores = {"DeNovoSeer": "mean_pred_prob", "clinical score": "Clinical_evidence_score",
              "PM2": "pm2", "BayesDel_addAF": "BayesDel_addAF_score",
              "AlphaMissense": "AlphaMissense_score", "CADD_phred": "CADD_phred"}
    for g, d in groups.items():
        for name, col in scores.items():
            rows.append(score_row(g, name, d["true_label"], d[col], rng))
    perf = pd.DataFrame(rows)
    perf.to_csv(os.path.join(args.out, "auc_by_label_source.csv"), index=False)

    say("1) AUC by label source (S1 test predictions, each variant counted once)")
    for g, d in groups.items():
        say(f"   {g}: n={len(d)}, pathogenic={int(d.true_label.sum())}")
    say(perf.round(3).to_string(index=False))

    # 2. missense only, ClinVar group
    cv_mis = groups["ClinVar"][groups["ClinVar"].consequence == "nonsynonymous SNV"]
    say("\n2) ClinVar group, missense only")
    for name, col in [("DeNovoSeer", "mean_pred_prob"), ("AlphaMissense", "AlphaMissense_score")]:
        r = score_row("ClinVar missense", name, cv_mis.true_label, cv_mis[col], rng)
        say(f"   {name}: n={r['n']}, AUC={r['AUC']:.3f} ({r['AUC_lo']:.3f}-{r['AUC_hi']:.3f})")

    cvb = groups["ClinVar"][groups["ClinVar"].true_label == 0]
    say(f"   ClinVar benign that are synonymous: {(cvb.consequence == 'synonymous SNV').sum()}"
        f"/{len(cvb)}")

    # 3. rule group, same variants for every tool
    rule = groups["rule"].dropna(subset=["AlphaMissense_score", "CADD_phred"])
    say(f"\n3) Rule group, variants scored by all three tools (n={len(rule)}, "
        f"pathogenic={int(rule.true_label.sum())})")
    for name, col in [("DeNovoSeer", "mean_pred_prob"), ("AlphaMissense", "AlphaMissense_score"),
                      ("CADD_phred", "CADD_phred")]:
        r = score_row("rule common", name, rule.true_label, rule[col], rng)
        say(f"   {name}: AUC={r['AUC']:.3f} ({r['AUC_lo']:.3f}-{r['AUC_hi']:.3f})")
    for name, col in [("AlphaMissense", "AlphaMissense_score"), ("CADD_phred", "CADD_phred")]:
        med, lo, hi = paired_diff(rule.true_label, rule.mean_pred_prob, rule[col], rng)
        say(f"   DeNovoSeer minus {name}: {med:.3f} ({lo:.3f}-{hi:.3f})")

    # 4. synonymous variants with a rule-based pathogenic label
    syn = groups["rule"][groups["rule"].consequence == "synonymous SNV"]
    syn_p, syn_b = syn[syn.true_label == 1], test[(test.consequence == "synonymous SNV") & (y == 0)]
    called = syn_p[syn_p.mean_pred_prob >= 0.5]
    say("\n4) Synonymous variants")
    say(f"   rule-pathogenic synonymous in test: {len(syn_p)}")
    say(f"   called pathogenic (p >= 0.5): {len(called)} ({len(called) / len(syn_p):.1%})")
    say(f"   benign synonymous called pathogenic: {(syn_b.mean_pred_prob >= 0.5).mean():.1%}")
    say(f"   called ones with no SpliceAI delta score: {called[SPLICE_DS].isna().all(axis=1).sum()}")
    syn_p.to_csv(os.path.join(args.out, "rule_pathogenic_synonymous.csv"), index=False)

    # 5. PS2 / PM6 in the released data
    say("\n5) PS2 / PM6 in the released clinical evidence file")
    say(f"   variants: {len(clin)}")
    for col in ["ps2", "pm6"]:
        counts = clin[col].astype(int).value_counts().to_dict()
        say(f"   {col.upper()} value counts: {counts}")

    # 6. repeated sample IDs behind rule-based pathogenic labels
    rp = ann[(ann.label_source == "rule") & (ann.label == 1)].copy()
    ids = rp["Otherinfo"].astype(str).str.split(";")
    rp["n_records"] = ids.apply(len)
    rp["n_unique_ids"] = ids.apply(lambda x: len(set(x)))
    dup = rp[rp.n_unique_ids == 1]
    say("\n6) Rule-based pathogenic labels, full Gene4Denovo file")
    say(f"   labels: {len(rp)}; all records share one sample ID: {len(dup)}")
    dup[["variant_index", "Otherinfo", "Phenotype", "Explanation"]].to_csv(
        os.path.join(args.out, "repeated_sample_ids.csv"), index=False)

    with open(os.path.join(args.out, "summary.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
