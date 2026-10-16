#!/usr/bin/env python3
"""
make_sample_table.py — one row per selected sample, with every characteristic
available in the public metadata (sex, body site, subject ID, tissue, ...).

Combines the SRA RunInfo table (sequencing details, and for HMP Subject_ID and
Sex) with all BioSample attributes (turned from one-row-per-attribute into one
column per attribute). Runs of the same sample are merged into one row.

Writes, in --outdir:
  LABEL_sample_table.tsv          the table (use it as phyloseq sample_data)
  LABEL_sample_table_dropped.tsv  attributes left out of the table: the same
                                  value in every sample, or database bookkeeping
"""
import argparse, collections, csv, os, re, sys

RUN_COLS = ["Subject_ID", "Sex", "Body_Site", "LibraryLayout", "Platform",
            "Model", "CenterName", "ScientificName"]
ADMIN = re.compile(r"^(ena-|insdc |ena_checklist|ena-checklist)", re.I)


def clean(name):
    return re.sub(r"[^0-9a-z]+", "_", name.strip().lower()).strip("_") or "attribute"


def main(a):
    keep = [l.split()[0] for l in open(a.samples) if l.strip() and not l.startswith("#")]
    keep_set = set(keep)
    runs = [r for r in csv.DictReader(open(a.runinfo, newline=""))
            if (r.get("Sample") or r.get("BioSample")) in keep_set]
    if not runs:
        sys.exit("No runs of the selected samples found in RunInfo")

    attrs = collections.defaultdict(dict)          # biosample -> {column: value}
    if a.biosample and os.path.exists(a.biosample):
        for r in csv.DictReader(open(a.biosample, newline=""), delimiter="\t"):
            name = r["harmonized_name"] or r["attribute_name"]
            attrs[r["biosample"]][name] = r["value"]

    rows = collections.OrderedDict()
    for r in runs:
        s = r.get("Sample") or r.get("BioSample")
        row = rows.setdefault(s, {"sample_id": s, "biosample": r.get("BioSample", ""),
                                  "runs": [], "n_runs": 0, "raw_bases": 0})
        row["runs"].append(r["Run"])
        row["n_runs"] += 1
        row["raw_bases"] += int(r.get("bases") or 0)
        for c in RUN_COLS:
            v = (r.get(c) or "").strip()
            if v and not row.get("sra_" + clean(c)):
                row["sra_" + clean(c)] = v
        for k, v in attrs.get(r.get("BioSample", ""), {}).items():
            row.setdefault("bs:" + k, v)

    cols = []
    for r in rows.values():
        for k in r:
            if k not in cols:
                cols.append(k)
    fixed = ["sample_id", "biosample", "runs", "n_runs", "raw_bases"]
    dropped, table_cols = [], list(fixed)
    for k in cols:
        if k in fixed:
            continue
        values = [str(r.get(k, "")) for r in rows.values()]
        name = k[3:] if k.startswith("bs:") else k
        if k.startswith("bs:") and ADMIN.match(name):
            dropped.append((name, "database bookkeeping", values[0] if len(set(values)) == 1 else "varies"))
        elif len(set(values)) == 1:
            dropped.append((name, "same in all samples", values[0]))
        else:
            table_cols.append(k)

    header = [c if not c.startswith("bs:") else clean(c[3:]) for c in table_cols]
    seen = collections.Counter()
    for i, h in enumerate(header):
        seen[h] += 1
        if seen[h] > 1:
            header[i] = f"{h}_{seen[h]}"

    os.makedirs(a.outdir, exist_ok=True)
    out = os.path.join(a.outdir, f"{a.label}_sample_table.tsv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        for s in keep:
            if s not in rows:
                continue
            r = rows[s]
            w.writerow([";".join(r[c]) if c == "runs" else r.get(c, "") for c in table_cols])
    with open(os.path.join(a.outdir, f"{a.label}_sample_table_dropped.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["attribute", "why_left_out", "value"])
        w.writerows(dropped)

    missing = [s for s in keep if s not in rows]
    print(f"{len(rows)} samples x {len(header)} columns -> {out}")
    print("columns: " + ", ".join(header))
    if missing:
        print(f"WARNING: {len(missing)} selected samples not found in RunInfo: {missing[:5]}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runinfo", required=True)
    p.add_argument("--biosample", default="")
    p.add_argument("--samples", required=True, help="one selected sample ID per line")
    p.add_argument("--label", required=True)
    p.add_argument("--outdir", required=True)
    main(p.parse_args())
