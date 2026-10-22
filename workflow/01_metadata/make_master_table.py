#!/usr/bin/env python3
"""
make_master_table.py — the sample metadata table: one row per sample with
everything known about it.

Columns: dataset, bioproject, sample_id, biosample, subject_id, sex, body_site,
runs, n_runs, raw_gbases, sequencer, layout, center; visit matching
(reference_sample = the person's tonsil sample, accession_distance,
same_visit_likely); KneadData read counts (raw_pairs, clean_pairs,
clean_unpaired_reads, pct_pairs_kept); then every BioSample attribute (bs_*).
Nothing is dropped, even if it is the same in all samples.

Make the table for one dataset:
  make_master_table.py --dataset hmp_oral --bioproject PRJNA48479 \
      --runs-full RUNS_FULL.tsv --biosample BIOSAMPLE_ATTRIBUTES.tsv \
      [--kneaddata kneaddata_read_counts.tsv] [--reference tonsil_runs_full.tsv] \
      --out OUT.tsv
Combine several datasets into one table (union of columns):
  make_master_table.py --combine OUT.tsv IN1.tsv IN2.tsv ...
"""
import argparse, collections, csv, re, sys

ADMIN = re.compile(r"^(ena-|insdc |ena_checklist|ena-checklist)", re.I)
SAME_VISIT_MAX = 1000          # accession distance treated as "same visit"


def clean(name):
    return re.sub(r"[^0-9a-z]+", "_", name.strip().lower()).strip("_")


def acc(s):
    d = re.sub(r"\D", "", s or "")
    return int(d) if d else None


def write(path, rows, cols):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def combine(out, inputs):
    rows, cols = [], []
    for f in inputs:
        for r in csv.DictReader(open(f), delimiter="\t"):
            rows.append(r)
            for c in r:
                if c not in cols:
                    cols.append(c)
    write(out, rows, cols)
    print(f"{len(rows)} samples x {len(cols)} columns -> {out}")


def build(a):
    runs = list(csv.DictReader(open(a.runs_full), delimiter="\t"))
    attrs = collections.defaultdict(dict)
    for r in csv.DictReader(open(a.biosample), delimiter="\t"):
        name = r["harmonized_name"] or r["attribute_name"]
        if not ADMIN.match(name):
            attrs[r["biosample"]][clean(name)] = r["value"]
    ref = {}
    if a.reference:
        ref = {r["subject_id"]: r["sample_id"] for r in csv.DictReader(open(a.reference), delimiter="\t")
               if r.get("subject_id")}
    kd = {}
    if a.kneaddata:
        for r in csv.DictReader(open(a.kneaddata), delimiter="\t"):
            kd[list(r.values())[0]] = r

    rows = collections.OrderedDict()
    for r in runs:
        s = r["sample_id"]
        x = rows.setdefault(s, {
            "dataset": a.dataset, "bioproject": a.bioproject, "sample_id": s,
            "biosample": r.get("biosample", ""), "subject_id": r.get("subject_id", ""),
            "body_site": re.sub(r"^G_DNA_", "", r.get("body_site", "")),
            "runs": [], "raw_bases": 0, "sequencer": set(), "layout": r.get("library_layout", ""),
            "center": r.get("center", "")})
        x["runs"].append(r["run"])
        x["raw_bases"] += int(r.get("bases") or 0)
        x["sequencer"].add(r.get("model", ""))

    bs_cols = []
    for x in rows.values():
        b = attrs.get(x["biosample"], {})
        x["sex"] = b.get("host_sex") or b.get("sex") or ""
        for k, v in b.items():
            x["bs_" + k] = v
            if "bs_" + k not in bs_cols:
                bs_cols.append("bs_" + k)
        x["n_runs"] = len(x["runs"])
        x["runs"] = ";".join(x["runs"])
        x["raw_gbases"] = f"{x.pop('raw_bases') / 1e9:.3f}"
        x["sequencer"] = ";".join(sorted(x["sequencer"]))
        rs = ref.get(x["subject_id"], "")
        x["reference_sample"] = rs
        if rs and acc(rs) is not None and acc(x["sample_id"]) is not None:
            d = abs(acc(x["sample_id"]) - acc(rs))
            x["accession_distance"] = d
            x["same_visit_likely"] = "yes" if d < SAME_VISIT_MAX else "no"
        k = kd.get(x["sample_id"])
        if k:
            raw = float(k.get("raw pair1") or 0)
            fin = float(k.get("final pair1") or 0)
            orph = float(k.get("final orphan1") or 0) + float(k.get("final orphan2") or 0)
            x["raw_pairs"] = int(raw)
            x["clean_pairs"] = int(fin)
            x["clean_unpaired_reads"] = int(orph)
            x["pct_pairs_kept"] = f"{100 * fin / raw:.1f}" if raw else ""

    cols = ["dataset", "bioproject", "sample_id", "biosample", "subject_id", "sex", "body_site",
            "runs", "n_runs", "raw_gbases", "sequencer", "layout", "center"]
    if ref:
        cols += ["reference_sample", "accession_distance", "same_visit_likely"]
    if kd:
        cols += ["raw_pairs", "clean_pairs", "clean_unpaired_reads", "pct_pairs_kept"]
    cols += bs_cols
    write(a.out, list(rows.values()), cols)
    missing_kd = [s for s in rows if kd and s not in kd]
    print(f"{len(rows)} samples x {len(cols)} columns -> {a.out}")
    if missing_kd:
        print(f"  note: no KneadData counts yet for {len(missing_kd)} samples (e.g. {missing_kd[:3]})")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--combine":
        combine(sys.argv[2], sys.argv[3:]); sys.exit()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", required=True)
    p.add_argument("--bioproject", required=True)
    p.add_argument("--runs-full", required=True)
    p.add_argument("--biosample", required=True)
    p.add_argument("--kneaddata", default="")
    p.add_argument("--reference", default="")
    p.add_argument("--out", required=True)
    build(p.parse_args())
