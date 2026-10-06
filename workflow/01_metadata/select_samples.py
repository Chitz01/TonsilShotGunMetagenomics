#!/usr/bin/env python3
"""
select_samples.py — pick runs from one body site (default: tonsil) using the
metadata downloaded by fetch_metadata.py.

A run is selected when the regex matches its SRA Body_Site value or any of its
BioSample attribute values (case-insensitive), and its LibraryStrategy is in
--strategy. Every match records WHICH field matched, so you can check that the
regex did not pick up anything unintended.

Outputs in --outdir (LABEL defaults to "tonsil"):
  LABEL_runs_full.tsv          selected runs with metadata and ENA FASTQ links
  LABEL_run_accessions.txt     one run accession per line   (-> 02_download)
  LABEL_sample_ids.txt         one sample (SRS) per line     (-> 03 onwards)
  LABEL_sample_runs.tsv        sample_id, run, layout, bases (runs per sample)
  LABEL_selection_report.txt   counts, bases, every distinct matched value
  LABEL_excluded_samples.tsv   samples removed by --one-per-subject, with reason
  body_site_summary.tsv        runs/samples per body site in the WHOLE project
                               (use it to choose oral comparison sites)
"""
import argparse
import collections
import csv
import os
import re
import sys

BODY_SITE_KEYS = ("body_site", "body site", "body_habitat", "body_product",
                  "env_material", "isolation_source", "tissue", "host_tissue_sampled")


def read_runinfo(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_biosample(path):
    attrs = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return attrs
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            name = r["harmonized_name"] or r["attribute_name"]
            attrs[r["biosample"]].append((name, r["value"]))
    return attrs


def read_ena(path):
    ena = {}
    if not path or not os.path.exists(path):
        return ena
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ena[r.get("run_accession", "")] = r
    return ena


def srs_number(sample_id):
    digits = re.sub(r"\D", "", sample_id)
    return int(digits) if digits else float("inf")


def visit_number(bs_attrs):
    """Integer visit number from BioSample attributes, if present."""
    for name, value in bs_attrs:
        if "visit" in name.lower():
            m = re.search(r"\d+", value or "")
            if m:
                return int(m.group())
    return None


def one_per_subject(selected, bs):
    """Keep one sample per subject. Rule (decided before seeing any results):
    1) lowest visit number, if every sample of that subject has one;
    2) otherwise the earliest-registered sample (lowest SRS accession number).
    Returns (kept_rows, excluded_records)."""
    by_subject = collections.defaultdict(set)
    for d in selected:
        if d["subject_id"]:
            by_subject[d["subject_id"]].add((d["sample_id"], d["biosample"]))
    drop, excluded = set(), []
    for subj, samples in sorted(by_subject.items()):
        if len(samples) < 2:
            continue
        visits = {s: visit_number(bs.get(b, [])) for s, b in samples}
        if all(v is not None for v in visits.values()):
            keep = min(visits, key=lambda s: (visits[s], srs_number(s)))
            rule = "kept lowest visit number"
        else:
            keep = min(visits, key=srs_number)
            rule = "no visit number in metadata; kept earliest-registered sample (lowest SRS)"
        for s in sorted(visits):
            if s != keep:
                drop.add(s)
                excluded.append({"sample_id": s, "subject_id": subj, "kept_sample": keep,
                                 "visit_removed": visits[s], "visit_kept": visits[keep],
                                 "reason": f"repeat sample of subject {subj}; {rule}"})
    return [d for d in selected if d["sample_id"] not in drop], excluded


def body_site_label(run, bs_attrs):
    """Best available body-site description for the summary table."""
    v = (run.get("Body_Site") or "").strip()
    if v:
        return v
    for name, value in bs_attrs:
        if any(k in name.lower() for k in BODY_SITE_KEYS) and value:
            return value
    return "unknown"


def main(a):
    rx = re.compile(a.regex, re.IGNORECASE)
    keep = {s.strip().upper() for s in a.strategy.split(",") if s.strip()}
    runs = read_runinfo(a.runinfo)
    bs = read_biosample(a.biosample)
    ena = read_ena(a.ena)
    if not runs:
        sys.exit("RunInfo table is empty")
    if not bs:
        print("WARNING: no BioSample attributes; matching on RunInfo only", file=sys.stderr)

    summary = collections.defaultdict(lambda: [0, set(), 0])
    selected, matched_values = [], collections.Counter()

    for r in runs:
        attrs = bs.get(r.get("BioSample", ""), [])
        site = body_site_label(r, attrs)
        strategy = (r.get("LibraryStrategy") or "").upper()
        sample = r.get("Sample") or r.get("BioSample") or r["Run"]
        bases = int(r.get("bases") or 0)

        s = summary[(site, strategy)]
        s[0] += 1
        s[1].add(sample)
        s[2] += bases

        if keep and strategy not in keep:
            continue
        if bases < a.min_bases:
            continue

        hits = []
        if rx.search(r.get("Body_Site") or ""):
            hits.append(("RunInfo:Body_Site", r["Body_Site"]))
        for name, value in attrs:
            if rx.search(value or ""):
                hits.append((f"BioSample:{name}", value))
        if not hits:
            continue
        for field, value in hits:
            matched_values[(field, value)] += 1

        e = ena.get(r["Run"], {})
        selected.append({
            "sample_id": sample,
            "run": r["Run"],
            "biosample": r.get("BioSample", ""),
            "body_site": site,
            "library_layout": r.get("LibraryLayout", ""),
            "library_strategy": r.get("LibraryStrategy", ""),
            "platform": r.get("Platform", ""),
            "model": r.get("Model", ""),
            "spots": r.get("spots", ""),
            "bases": bases,
            "avg_length": r.get("avgLength", ""),
            "size_MB": r.get("size_MB", ""),
            "center": r.get("CenterName", ""),
            "subject_id": r.get("Subject_ID", ""),
            "matched_fields": "; ".join(f"{f}={v}" for f, v in hits),
            "ena_fastq_ftp": e.get("fastq_ftp", ""),
            "ena_fastq_md5": e.get("fastq_md5", ""),
            "ena_fastq_bytes": e.get("fastq_bytes", ""),
        })

    os.makedirs(a.outdir, exist_ok=True)
    L = a.label

    # whole-project body-site summary
    with open(os.path.join(a.outdir, "body_site_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["body_site", "library_strategy", "n_runs", "n_samples", "total_Gbases"])
        for (site, strat), (nr, ss, b) in sorted(summary.items(), key=lambda x: -x[1][0]):
            w.writerow([site, strat, nr, len(ss), f"{b / 1e9:.2f}"])

    if not selected:
        sys.exit(f"No runs matched /{a.regex}/ — check body_site_summary.tsv for the exact wording")

    n_before = len({d["sample_id"] for d in selected})
    excluded = []
    if a.one_per_subject:
        selected, excluded = one_per_subject(selected, bs)
    with open(os.path.join(a.outdir, f"{a.label}_excluded_samples.tsv"), "w", newline="") as fh:
        cols = ["sample_id", "subject_id", "kept_sample", "visit_removed", "visit_kept", "reason"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(excluded)

    selected.sort(key=lambda d: (d["sample_id"], d["run"]))
    with open(os.path.join(a.outdir, f"{L}_runs_full.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(selected[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(selected)

    with open(os.path.join(a.outdir, f"{L}_run_accessions.txt"), "w") as fh:
        fh.writelines(f"{d['run']}\n" for d in selected)

    samples = sorted({d["sample_id"] for d in selected})
    with open(os.path.join(a.outdir, f"{L}_sample_ids.txt"), "w") as fh:
        fh.writelines(f"{s}\n" for s in samples)

    with open(os.path.join(a.outdir, f"{L}_sample_runs.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "run", "library_layout", "bases"])
        for d in selected:
            w.writerow([d["sample_id"], d["run"], d["library_layout"], d["bases"]])

    # human-readable report
    per_sample = collections.defaultdict(lambda: [0, 0])
    for d in selected:
        per_sample[d["sample_id"]][0] += 1
        per_sample[d["sample_id"]][1] += d["bases"]
    layouts = collections.Counter(d["library_layout"] for d in selected)
    subjects = {d["subject_id"] for d in selected if d["subject_id"]}
    with open(os.path.join(a.outdir, f"{L}_selection_report.txt"), "w") as fh:
        fh.write(f"Selection label : {L}\nRegex           : /{a.regex}/i\n")
        fh.write(f"Strategies kept : {', '.join(sorted(keep)) or 'all'}\n")
        fh.write(f"Min bases/run   : {a.min_bases}\n\n")
        fh.write(f"One per subject : {'yes' if a.one_per_subject else 'no'}\n")
        fh.write(f"Samples matched : {n_before}\n")
        fh.write(f"Samples removed : {len(excluded)} (see {L}_excluded_samples.tsv)\n")
        for e in excluded:
            fh.write(f"    - {e['sample_id']} (subject {e['subject_id']}; kept {e['kept_sample']}): {e['reason']}\n")
        fh.write(f"Runs selected   : {len(selected)}\n")
        fh.write(f"Samples         : {len(samples)}\n")
        fh.write(f"Subjects        : {len(subjects) if subjects else 'not reported'}\n")
        fh.write(f"Total size      : {sum(d['bases'] for d in selected) / 1e9:.2f} Gbases\n")
        fh.write(f"Library layouts : {dict(layouts)}\n\n")
        fh.write("Distinct matched values (CHECK THESE ARE ALL TONSIL):\n")
        for (field, value), n in matched_values.most_common():
            fh.write(f"  {n:5d}  {field} = {value}\n")
        fh.write("\nPer sample (runs, Gbases):\n")
        for s in samples:
            n, b = per_sample[s]
            fh.write(f"  {s}\t{n}\t{b / 1e9:.2f}\n")

    print(f"Selected {len(selected)} runs from {len(samples)} samples -> {a.outdir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runinfo", required=True)
    p.add_argument("--biosample", default="")
    p.add_argument("--ena", default="")
    p.add_argument("--regex", default="tonsil")
    p.add_argument("--label", default="tonsil")
    p.add_argument("--strategy", default="WGS", help="comma list; empty = keep all")
    p.add_argument("--min-bases", type=int, default=0)
    p.add_argument("--one-per-subject", action="store_true",
                   help="keep one sample per subject (see one_per_subject())")
    p.add_argument("--outdir", required=True)
    main(p.parse_args())
