import csv, os, re, collections
B = '/lunarc/nobackup/projects/lu2025-12-46/CHETANA/ANALYSIS/'
SETS = [('tonsil', 'hmp_tonsil'), ('oral', 'hmp_oral')]

def kept(d):
    q = B + d + '/07_mags/quality_filtered'
    return {f[:-3] for f in os.listdir(q) if f.endswith('.fa')}

def gtdb(d, k):
    out = {}
    for f in ('gtdbtk.bac120.summary.tsv', 'gtdbtk.ar53.summary.tsv'):
        p = B + d + '/07_mags/gtdbtk/' + f
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter='\t'):
                if r['user_genome'] in k:
                    out[r['user_genome']] = r['classification']
    return out

def rank(c, r):
    m = re.search(r + r'__([^;]*)', c)
    return m.group(1) if m and m.group(1) else ''

def near_complete(d, k):
    n = 0
    for r in csv.DictReader(open(B + d + '/07_mags/mag_quality_summary.tsv'), delimiter='\t'):
        if r['bin'] in k and float(r['completeness']) >= 90:
            n += 1
    return n

def bakta_done(d, k):
    return sum(1 for m in k if os.path.getsize(B + d + f'/07_mags/bakta/{m}/{m}.gff3') > 0
               if os.path.exists(B + d + f'/07_mags/bakta/{m}/{m}.gff3'))

res, gen = {}, {}
for name, d in SETS:
    k = kept(d); g = gtdb(d, k)
    genera = collections.Counter(rank(c, 'g') or 'unclassified' for c in g.values())
    species = {rank(c, 's') for c in g.values() if rank(c, 's')}
    res[name] = dict(mags=len(k), gtdb=len(g), near=near_complete(d, k), bakta=bakta_done(d, k),
                     genera=len([x for x in genera if x != 'unclassified']), species=len(species),
                     nospecies=sum(1 for c in g.values() if not rank(c, 's')),
                     archaea=sum(1 for c in g.values() if c.startswith('d__Archaea')))
    gen[name] = genera

rows = [('MAGs passing filter', 'mags'), ('classified by GTDB-Tk', 'gtdb'), ('annotated by Bakta', 'bakta'),
        ('near-complete (>=90%)', 'near'), ('distinct genera', 'genera'), ('distinct named species', 'species'),
        ('MAGs without species match', 'nospecies'), ('archaeal MAGs', 'archaea')]
print(f"{'':30s} {'tonsil':>8s} {'oral':>8s}")
for label, key in rows:
    print(f"{label:30s} {res['tonsil'][key]:8d} {res['oral'][key]:8d}")

t, o = set(gen['tonsil']) - {'unclassified'}, set(gen['oral']) - {'unclassified'}
print(f"\ngenera: {len(t & o)} shared, {len(t - o)} only in tonsil MAGs, {len(o - t)} only in oral MAGs")
print("only in tonsil MAGs:", ', '.join(sorted(t - o)) or '-')
print(f"\n{'genus':28s} {'tonsil':>7s} {'oral':>7s}   (top 20 by total)")
for g_, n in (gen['tonsil'] + gen['oral']).most_common(20):
    print(f"{g_:28s} {gen['tonsil'][g_]:7d} {gen['oral'][g_]:7d}")
