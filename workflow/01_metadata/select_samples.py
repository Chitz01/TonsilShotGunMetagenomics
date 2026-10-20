#!/usr/bin/env python3
"""
select_samples.py — pick runs from the body site(s) of interest using the
metadata downloaded by fetch_metadata.py.

A run is selected when the regex matches its SRA Body_Site value or any of its
BioSample attribute values (case-insensitive), its LibraryStrategy is in
--strategy, and it passes the optional filters below. Every match records WHICH
field matched, so the selection can be checked by eye.

Optional filters:
  --subjects FILE     keep only these subject IDs (one per line)
  --platform REGEX    keep only runs whose Platform/Model match (e.g. Illumina)
  --min-bases N       drop runs with fewer bases
  --one-per-subject {none,subject,site}
        subject: one sample per person in total (tonsil dataset)
        site:    one sample per person AND body site (matched oral sites)
        Rule, fixed in advance: lowest visit number if every candidate has one,
        otherwise the earliest-registered sample (lowest SRS/ERS number).

Outputs in --outdir (LABEL = --label):
  LABEL_runs_full.tsv          selected runs with metadata and ENA FASTQ links
  LABEL_run_accessions.txt     one run accession per line   (-> 02_download)
  LABEL_sample_ids.txt         one sample per line           (-> 03 onwards)
  LABEL_sample_runs.tsv        sample_id, run, layout, bases
  LABEL_excluded_samples.tsv   samples removed by --one-per-subject, with reason
  LABEL_selection_report.txt   counts, filters, matched values, per-site table
  body_site_summary.tsv        runs/samples per body site in the WHOLE project
"""
import argparse, collections, csv, os, re, sys

BODY_SITE_KEYS = ("body_site", "body site", "body_habitat", "body_product",
                  "env_material", "isolation_source", "tissue", "host_tissue_sampled")


def read_runinfo(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_biosample(path):
    attrs = collections.defaultdict(list)
    if path and os.path.exists(path):
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                attrs[r["biosample"]].append((r["harmonized_name"] or r["attribute_name"], r["value"]))
    return attrs


def read_ena(path):
    ena = {}
    if path and os.path.exists(path):
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                ena[r.get("run_accession", "")] = r
    return ena


def body_site_label(run, bs_attrs):
    v = (run.get("Body_Site") or "").strip()
    if v:
        return v
    for name, value in bs_attrs:
        if any(k in name.lower() for k in BODY_SITE_KEYS) and value:
            return value
    return "unknown"


def acc_number(sample_id):
    digits = re.sub(r"\D", "", sample_id)
    return int(digits) if digits else float("inf")


def visit_number(bs_attrs):
    for name, value in bs_attrs:
        if "visit" in name.lower():
            m = re.search(r"\d+", value or "")
            if m:
                return int(m.group())
    return None


def one_per(selected, bs, mode):
    """Keep one sample per subject (mode 'subject') or per subject+site ('site')."""
    groups = collections.defaultdict(set)
    for d in selected:
        if d["subject_id"]:
            key = (d["subject_id"], d["body_site"]) if mode == "site" else (d["subject_id"],)
            groups[key].add((d["sample_id"], d["biosample"]))
    drop, excluded = set(), []
    for key, samples in sorted(groups.items()):
        if len(samples) < 2:
            continue
        visits = {s: visit_number(bs.get(b, [])) for s, b in samples}
        if all(v is not None for v in visits.values()):
            keep = min(visits, key=lambda s: (visits[s], acc_number(s)))
            rule = "kept lowest visit number"
        else:
            keep = min(visits, key=acc_number)
            rule = "no visit number in metadata; kept earliest-registered sample (lowest accession)"
        what = f"subject {key[0]}" + (f" at {key[1]}" if mode == "site" else "")
        for s in sorted(visits):
            if s != keep:
                drop.add(s)
                excluded.append({"sample_id": s, "subject_id": key[0],
                                 "body_site": key[1] if mode == "site" else "",
                                 "kept_sample": keep, "visit_removed": visits[s],
                                 "visit_kept": visits[keep],
                                 "reason": f"repeat sample of {what}; {rule}"})
    return [d for d in selected if d["sample_id"] not in drop], excluded


def main(a):
    rx = re.compile(a.regex, re.IGNORECASE)
    prx = re.compile(a.platform, re.IGNORECASE) if a.platform else None
    keep_strat = {s.strip().upper() for s in a.strategy.split(",") if s.strip()}
    subjects = None
    if a.subjects:
        subjects = {l.strip() for l in open(a.subjects) if l.strip() and not l.startswith("#")}
    runs, bs, ena = read_runinfo(a.runinfo), read_biosample(a.biosample), read_ena(a.ena)
    if not runs:
        sys.exit("RunInfo table is empty")

    summary = collections.defaultdict(lambda: [0, set(), 0])
    selected, matched_values = [], collections.Counter()
    filtered = collections.Counter()
    for r in runs:
        attrs = bs.get(r.get("BioSample", ""), [])
        site = body_site_label(r, attrs)
        strategy = (r.get("LibraryStrategy") or "").upper()
        sample = r.get("Sample") or r.get("BioSample") or r["Run"]
        bases = int(r.get("bases") or 0)
        subj = (r.get("Subject_ID") or "").strip()
        s = summary[(site, strategy)]
        s[0] += 1; s[1].add(sample); s[2] += bases

        hits = []
        if rx.search(r.get("Body_Site") or ""):
            hits.append(("RunInfo:Body_Site", r["Body_Site"]))
        for name, value in attrs:
            if rx.search(value or ""):
                hits.append((f"BioSample:{name}", value))
        if not hits:
            continue
        if keep_strat and strategy not in keep_strat:
            filtered["library strategy not " + "/".join(sorted(keep_strat))] += 1; continue
        if subjects is not None and subj not in subjects:
            filtered["subject not in --subjects list"] += 1; continue
        if prx and not prx.search(f"{r.get('Platform', '')} {r.get('Model', '')}"):
            filtered[f"platform not /{a.platform}/"] += 1; continue
        if bases < a.min_bases:
            filtered[f"fewer than {a.min_bases} bases"] += 1; continue
        for field, value in hits:
            matched_values[(field, value)] += 1
        e = ena.get(r["Run"], {})
        selected.append({
            "sample_id": sample, "run": r["Run"], "biosample": r.get("BioSample", ""),
            "body_site": site, "library_layout": r.get("LibraryLayout", ""),
            "library_strategy": r.get("LibraryStrategy", ""), "platform": r.get("Platform", ""),
            "model": r.get("Model", ""), "spots": r.get("spots", ""), "bases": bases,
            "avg_length": r.get("avgLength", ""), "size_MB": r.get("size_MB", ""),
            "center": r.get("CenterName", ""), "subject_id": subj,
            "matched_fields": "; ".join(f"{f}={v}" for f, v in hits),
            "ena_fastq_ftp": e.get("fastq_ftp", ""), "ena_fastq_md5": e.get("fastq_md5", ""),
            "ena_fastq_bytes": e.get("fastq_bytes", ""),
        })

    os.makedirs(a.outdir, exist_ok=True)
    L = a.label
    with open(os.path.join(a.outdir, "body_site_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["body_site", "library_strategy", "n_runs", "n_samples", "total_Gbases"])
        for (site, strat), (nr, ss, b) in sorted(summary.items(), key=lambda x: -x[1][0]):
            w.writerow([site, strat, nr, len(ss), f"{b / 1e9:.2f}"])
    if not selected:
        sys.exit(f"No runs left after matching /{a.regex}/ and filters {dict(filtered)}; check body_site_summary.tsv")

    n_before = len({d["sample_id"] for d in selected})
    excluded = []
    if a.one_per_subject in ("subject", "site"):
        selected, excluded = one_per(selected, bs, a.one_per_subject)

    with open(os.path.join(a.outdir, f"{L}_excluded_samples.tsv"), "w", newline="") as fh:
        cols = ["sample_id", "subject_id", "body_site", "kept_sample", "visit_removed", "visit_kept", "reason"]
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(excluded)

    selected.sort(key=lambda d: (d["sample_id"], d["run"]))
    with open(os.path.join(a.outdir, f"{L}_runs_full.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(selected[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(selected)
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

    per_sample = collections.defaultdict(lambda: [0, 0, "", ""])
    for d in selected:
        p = per_sample[d["sample_id"]]
        p[0] += 1; p[1] += d["bases"]; p[2] = d["body_site"]; p[3] = d["subject_id"]
    per_site = collections.defaultdict(lambda: [set(), set(), 0])
    for s, (n, b, site, subj) in per_sample.items():
        x = per_site[site]; x[0].add(s); x[2] += b
        if subj: x[1].add(subj)
    layouts = collections.Counter(d["library_layout"] for d in selected)
    models = collections.Counter(d["model"] for d in selected)
    subj_all = {d["subject_id"] for d in selected if d["subject_id"]}
    with open(os.path.join(a.outdir, f"{L}_selection_report.txt"), "w") as fh:
        fh.write(f"Selection label : {L}\nRegex           : /{a.regex}/i\n")
        fh.write(f"Strategies kept : {', '.join(sorted(keep_strat)) or 'all'}\n")
        fh.write(f"Subjects filter : {a.subjects + f' ({len(subjects)} IDs)' if subjects is not None else 'none'}\n")
        fh.write(f"Platform filter : {'/' + a.platform + '/i' if a.platform else 'none'}\n")
        fh.write(f"Min bases/run   : {a.min_bases}\n")
        fh.write(f"One sample per  : {a.one_per_subject}\n\n")
        for why, n in filtered.items():
            fh.write(f"Runs dropped    : {n} ({why})\n")
        fh.write(f"Samples matched : {n_before}\n")
        fh.write(f"Samples removed : {len(excluded)} (see {L}_excluded_samples.tsv)\n")
        for e in excluded:
            fh.write(f"    - {e['sample_id']} (kept {e['kept_sample']}): {e['reason']}\n")
        fh.write(f"Runs selected   : {len(selected)}\nSamples         : {len(samples)}\n")
        fh.write(f"Subjects        : {len(subj_all) if subj_all else 'not reported'}\n")
        fh.write(f"Total size      : {sum(d['bases'] for d in selected) / 1e9:.2f} Gbases\n")
        fh.write(f"Library layouts : {dict(layouts)}\nSequencers      : {dict(models)}\n\n")
        fh.write("Per body site (samples, people, Gbases):\n")
        for site, (ss, pp, b) in sorted(per_site.items(), key=lambda x: -len(x[1][0])):
            fh.write(f"  {len(ss):4d} {len(pp):4d} {b / 1e9:8.2f}  {site}\n")
        fh.write("\nDistinct matched values (CHECK THESE ARE THE INTENDED BODY SITES):\n")
        for (field, value), n in matched_values.most_common():
            fh.write(f"  {n:5d}  {field} = {value}\n")
        fh.write("\nPer sample (runs, Gbases, subject, body site):\n")
        for s in samples:
            n, b, site, subj = per_sample[s]
            fh.write(f"  {s}\t{n}\t{b / 1e9:.2f}\t{subj}\t{site}\n")
    print(f"Selected {len(selected)} runs from {len(samples)} samples -> {a.outdir}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runinfo", required=True)
    p.add_argument("--biosample", default="")
    p.add_argument("--ena", default="")
    p.add_argument("--regex", default="tonsil")
    p.add_argument("--label", default="tonsil")
    p.add_argument("--strategy", default="WGS", help="comma list; empty = keep all")
    p.add_argument("--subjects", default="", help="file with subject IDs to keep")
    p.add_argument("--platform", default="", help="regex on Platform + Model, e.g. Illumina")
    p.add_argument("--min-bases", type=int, default=0)
    p.add_argument("--one-per-subject", choices=["none", "subject", "site"], default="none")
    p.add_argument("--outdir", required=True)
    main(p.parse_args())
