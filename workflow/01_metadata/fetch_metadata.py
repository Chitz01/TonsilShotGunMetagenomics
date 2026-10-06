#!/usr/bin/env python3
"""
fetch_metadata.py — download run- and sample-level metadata for a BioProject.

Downloads metadata only (a few MB), never sequencing data. Produces:

  <outdir>/<BIOPROJECT>_sra_runinfo.csv        NCBI SRA RunInfo (one row per run;
                                               includes Body_Site when available)
  <outdir>/<BIOPROJECT>_biosample_attributes.tsv  every BioSample attribute, long format
  <outdir>/<BIOPROJECT>_ena_read_run.tsv       ENA file report (FASTQ URLs + MD5)
  <outdir>/fetch_log.txt                       source URLs and timestamp

Only the Python standard library is used. NCBI asks for at most 3 requests per
second without an API key; set NCBI_API_KEY in the environment to go faster.

Usage:
  python3 fetch_metadata.py --bioproject PRJNA48479 --outdir <dir> [--email you@x]
"""
import argparse
import csv
import datetime as dt
import io
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
API_KEY = os.environ.get("NCBI_API_KEY", "")
PAUSE = 0.12 if API_KEY else 0.4

ENA_FIELDS = [
    "run_accession", "sample_accession", "secondary_sample_accession",
    "experiment_accession", "study_accession", "sample_alias",
    "library_layout", "library_strategy", "library_source",
    "instrument_model", "read_count", "base_count",
    "fastq_ftp", "fastq_md5", "fastq_bytes",
]
ENA_FIELDS_MINIMAL = ["run_accession", "fastq_ftp", "fastq_md5", "fastq_bytes"]


def log(msg):
    print(f"[{dt.datetime.now():%Y-%m-%d %H:%M:%S}] {msg}", file=sys.stderr)


def request(url, params=None, post=False, retries=5):
    """GET/POST with retries and polite pausing. Returns decoded text."""
    params = dict(params or {})
    if url.startswith(EUTILS):
        params.setdefault("tool", "binp52_pipeline")
        if ARGS.email:
            params.setdefault("email", ARGS.email)
        if API_KEY:
            params.setdefault("api_key", API_KEY)
    data = None
    if post:
        data = urllib.parse.urlencode(params).encode()
        full = url
    else:
        full = f"{url}?{urllib.parse.urlencode(params)}" if params else url
    last = None
    for attempt in range(1, retries + 1):
        try:
            time.sleep(PAUSE)
            with urllib.request.urlopen(full, data=data, timeout=180) as r:
                return r.read().decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            last = e
            if e.code == 400:   # bad request: retrying will not help
                raise
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = e
        wait = 2 ** attempt
        log(f"  request failed ({last}); retry {attempt}/{retries} in {wait}s")
        time.sleep(wait)
    raise RuntimeError(f"Giving up on {url}: {last}")


def esearch(db, term):
    txt = request(f"{EUTILS}/esearch.fcgi",
                  {"db": db, "term": term, "usehistory": "y",
                   "retmax": 0, "retmode": "json"}, post=True)
    res = json.loads(txt)["esearchresult"]
    return int(res.get("count", 0)), res.get("webenv"), res.get("querykey")


# ---------------------------------------------------------------------------
# 1. SRA RunInfo
# ---------------------------------------------------------------------------
def fetch_runinfo(bioproject, path):
    count, webenv, qkey = esearch("sra", f"{bioproject}[BioProject]")
    if count == 0:
        log("  no hits with [BioProject] field; retrying as free text")
        count, webenv, qkey = esearch("sra", bioproject)
    if count == 0:
        raise SystemExit(f"No SRA records found for {bioproject}")
    log(f"SRA: {count} experiment records for {bioproject}")

    header, rows, step = None, [], 2000
    for start in range(0, count, step):
        txt = request(f"{EUTILS}/efetch.fcgi",
                      {"db": "sra", "rettype": "runinfo", "retmode": "text",
                       "WebEnv": webenv, "query_key": qkey,
                       "retstart": start, "retmax": step})
        reader = csv.reader(io.StringIO(txt))
        for rec in reader:
            if not rec or not any(rec):
                continue
            if rec[0] == "Run":           # header (repeated per chunk)
                header = header or rec
                continue
            rows.append(rec)
        log(f"  runinfo: {min(start + step, count)}/{count} records fetched")
    if header is None:
        raise SystemExit("RunInfo download returned no header; NCBI may be unavailable")

    seen, unique = set(), []
    for r in rows:                         # de-duplicate runs
        if r[0] not in seen:
            seen.add(r[0])
            unique.append(r)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(unique)
    log(f"Wrote {len(unique)} runs -> {path}")
    return header, unique


# ---------------------------------------------------------------------------
# 2. BioSample attributes
# ---------------------------------------------------------------------------
def fetch_biosamples(accessions, path, batch=100):
    out = open(path, "w", newline="")
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(["biosample", "sra_sample", "attribute_name",
                "harmonized_name", "value"])
    accessions = sorted(a for a in set(accessions) if a)
    n_written = 0
    for i in range(0, len(accessions), batch):
        chunk = accessions[i:i + batch]
        term = " OR ".join(f"{a}[accn]" for a in chunk)
        count, webenv, qkey = esearch("biosample", term)
        if count == 0:
            log(f"  biosample batch {i // batch + 1}: no hits")
            continue
        xml_txt = request(f"{EUTILS}/efetch.fcgi",
                          {"db": "biosample", "retmode": "xml",
                           "WebEnv": webenv, "query_key": qkey,
                           "retmax": count})
        root = ET.fromstring(xml_txt)
        for bs in root.iter("BioSample"):
            acc = bs.get("accession", "")
            srs = ""
            for id_el in bs.iter("Id"):
                if id_el.get("db") == "SRA":
                    srs = (id_el.text or "").strip()
            title = bs.findtext("Description/Title", default="") or ""
            if title:
                w.writerow([acc, srs, "title", "title", title.strip()])
            for att in bs.iter("Attribute"):
                w.writerow([acc, srs,
                            att.get("attribute_name", ""),
                            att.get("harmonized_name", ""),
                            (att.text or "").strip()])
            n_written += 1
        log(f"  biosamples: {min(i + batch, len(accessions))}/{len(accessions)}")
    out.close()
    log(f"Wrote attributes for {n_written} BioSamples -> {path}")


# ---------------------------------------------------------------------------
# 3. ENA file report (FASTQ URLs and MD5 sums; used by DOWNLOAD_METHOD=ena)
# ---------------------------------------------------------------------------
def fetch_ena(bioproject, path):
    for fields in (ENA_FIELDS, ENA_FIELDS_MINIMAL):
        try:
            txt = request(ENA_FILEREPORT,
                          {"accession": bioproject, "result": "read_run",
                           "fields": ",".join(fields), "format": "tsv",
                           "download": "true", "limit": 0})
            with open(path, "w") as fh:
                fh.write(txt)
            n = max(0, len(txt.strip().splitlines()) - 1)
            log(f"Wrote ENA report ({n} runs) -> {path}")
            return
        except urllib.error.HTTPError as e:
            log(f"  ENA request with {len(fields)} fields failed ({e}); trying fewer fields")
    log("WARNING: ENA report could not be downloaded (only needed for DOWNLOAD_METHOD=ena)")


def main():
    os.makedirs(ARGS.outdir, exist_ok=True)
    bp = ARGS.bioproject
    runinfo = os.path.join(ARGS.outdir, f"{bp}_sra_runinfo.csv")
    biosamp = os.path.join(ARGS.outdir, f"{bp}_biosample_attributes.tsv")
    ena = os.path.join(ARGS.outdir, f"{bp}_ena_read_run.tsv")

    header, rows = fetch_runinfo(bp, runinfo)
    bs_col = header.index("BioSample") if "BioSample" in header else None
    if bs_col is not None and not ARGS.skip_biosample:
        fetch_biosamples([r[bs_col] for r in rows if len(r) > bs_col], biosamp)
    if not ARGS.skip_ena:
        fetch_ena(bp, ena)

    with open(os.path.join(ARGS.outdir, "fetch_log.txt"), "w") as fh:
        fh.write(f"date\t{dt.datetime.now().isoformat()}\n")
        fh.write(f"bioproject\t{bp}\n")
        fh.write(f"ncbi_eutils\t{EUTILS}\nena_portal\t{ENA_FILEREPORT}\n")
    log("Metadata download finished.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bioproject", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--email", default="", help="sent to NCBI as contact (optional)")
    p.add_argument("--skip-biosample", action="store_true")
    p.add_argument("--skip-ena", action="store_true")
    ARGS = p.parse_args()
    main()
