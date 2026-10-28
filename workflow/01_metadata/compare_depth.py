import csv, collections, re, statistics as st
B = '/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/'
ORAL = re.compile(r'gingiva|buccal|palate|tonsil|saliva|plaque|throat|tongue|mucosa|rinse', re.I)

def gb(v): return v / 1e9
def summ(xs): xs = sorted(xs); return f"{st.median(xs):6.2f} {xs[0]:6.2f} {xs[-1]:6.2f}"

def per_sample(runinfo, site_of, keep_sample=None):
    d = {}
    for r in csv.DictReader(open(runinfo)):
        if r['LibraryStrategy'] != 'WGS':
            continue
        s = r['Sample']
        if keep_sample and s not in keep_sample:
            continue
        site = site_of.get(r['BioSample'], 'unknown')
        x = d.setdefault(s, [site, r.get('Subject_ID', ''), 0, r['LibraryLayout'], r['Model'], []])
        x[2] += int(r['bases'] or 0)
        x[5].append(float(r['avgLength'] or 0))
    return d

def attr(path, name):
    m = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        if name in (r['attribute_name'], r['harmonized_name']):
            m[r['biosample']] = r['value']
    return m

hdr = f"{'dataset / site':40s} {'samples':>7s} {'people':>6s} {'+tonsil':>7s}  {'raw Gb/sample: median   min    max':34s} {'layout':7s} sequencer"
print(hdr); print('-' * len(hdr))

hmp_site = attr(B + '01_metadata/raw/PRJNA48479_biosample_attributes.tsv', 'isolation_source')
hmp = per_sample(B + '01_metadata/raw/PRJNA48479_sra_runinfo.csv', hmp_site)
tonsil_ids = set(open(B + '01_metadata/selected/tonsil_sample_ids.txt').read().split())
tonsil_people = {v[1] for k, v in hmp.items() if k in tonsil_ids}
by_site = collections.defaultdict(list)
for s, v in hmp.items():
    if ORAL.search(v[0]):
        by_site[v[0]].append(v)
for site, vs in sorted(by_site.items(), key=lambda x: -len(x[1])):
    people = {v[1] for v in vs if v[1]}
    models = collections.Counter(v[4] for v in vs)
    print(f"{'HMP  ' + site.replace('G_DNA_', ''):40s} {len(vs):7d} {len(people):6d} {len(people & tonsil_people):7d}  "
          f"{summ([gb(v[2]) for v in vs]):34s} {vs[0][3]:7s} {dict(models)}")

sel = [v for k, v in hmp.items() if k in tonsil_ids]
print(f"{'YOUR 23 tonsil samples (raw)':40s} {len(sel):7d} {len({v[1] for v in sel}):6d} {'':7s}  {summ([gb(v[2]) for v in sel]):34s}")
kd = B + '03_preprocessing/kneaddata/kneaddata_read_counts.tsv'
try:
    rows = list(csv.DictReader(open(kd), delimiter='\t'))
    clean = [2 * float(r['final pair1']) * 100 / 1e9 for r in rows]
    print(f"{'YOUR 23 tonsil samples (CLEAN, <=)':40s} {len(clean):7d} {'':6s} {'':7s}  {summ(clean):34s} (clean pairs x 2 x 100 bp; human reads removed)")
except FileNotFoundError:
    print('  (KneadData table not found)')

fe_site = attr(B + 'oral_PRJEB36291/01_metadata/raw/PRJEB36291_biosample_attributes.tsv', 'tissue')
fe = per_sample(B + 'oral_PRJEB36291/01_metadata/raw/PRJEB36291_sra_runinfo.csv', fe_site)
by_site = collections.defaultdict(list)
for s, v in fe.items():
    by_site[v[0]].append(v)
allv = []
for site, vs in sorted(by_site.items()):
    allv += vs
    rl = st.median([l for v in vs for l in v[5]])
    print(f"{'Ferrara ' + site:40s} {len(vs):7d} {'':6s} {'':7s}  {summ([gb(v[2]) for v in vs]):34s} {vs[0][3]:7s} {vs[0][4]}, median read {rl:.0f} bp")
print(f"{'Ferrara ALL (before human removal)':40s} {len(allv):7d} {'20':>6s} {'':7s}  {summ([gb(v[2]) for v in allv]):34s}")
