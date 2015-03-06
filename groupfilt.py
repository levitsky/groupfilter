from pyteomics.auxiliary import filter, fdr
from csv import DictReader
from itertools import groupby
from sys import argv

filename = argv[1]
with open(filename) as f:
    psms = list(DictReader(f, delimiter='\t'))
key = lambda x: (int(float(x['Morpheus Score'])))
score = lambda x: float(x['Morpheus Score'])
is_decoy = lambda x: x['Decoy?'] == 'True' and x['Target?'] == 'False'

def groupfilter(psms, key, score, fdr, is_decoy, reverse=False, remove_decoy=False):
    out = []
    bad = []
    gd = []
    ndec = 0
    flag = 1
    psms = sorted(psms, key=key, reverse=True)
    for k, g in groupby(psms, key=key):
        g = list(g)
        if flag:
            if not sum(map(is_decoy, g)):
                gd.extend(g)
            else:
                flag = 0
        if not flag and not all(is_decoy(x) for x in g):
            with filter(g,
                    key=score, fdr=fdr, is_decoy=is_decoy, reverse=reverse, remove_decoy=remove_decoy, correction=2, full_output=False) as good:
                goodl = list(good)
                if len(goodl) == len(g):
                    gd.extend(g)
                elif len(goodl):
                    out.extend(goodl)
                elif not goodl:
                    bad.extend(g)                        
    if bad:
        with filter(bad,
            key=score, fdr=fdr, is_decoy=is_decoy, reverse=reverse, remove_decoy=remove_decoy, correction=2, full_output=False) as b:
            b = list(b)
            out.extend(b)
    if gd:
        with filter(gd,
            key=score, fdr=fdr, is_decoy=is_decoy, reverse=reverse, remove_decoy=remove_decoy, correction=2, full_output=False) as b:
            b = list(b)
            out.extend(b)
    return out

def iterate(fdr_v):
    fdr_c = -1
    fdr_u = fdr_v
    prev = fdr_u
    next = fdr_u * 2
    while abs(fdr_c - fdr_v) >= fdr_v * 0.1:
        if fdr_v > fdr_c:
            fdr_u = fdr_u + fdr_u / 2
        else:
            fdr_u = fdr_u - fdr_u / 2
        grf = groupfilter(psms, key, score, fdr_u, is_decoy, reverse=True)
        fdr_c = fdr(grf, is_decoy=is_decoy, formula=2)

    return fdr_u
fdr_c = iterate(0.01)
grf = groupfilter(psms, key, score, fdr_c, is_decoy, reverse=True)

fdr_n = fdr(grf, is_decoy=is_decoy, formula=2)
f = filter(psms, fdr=fdr_n, key=score, reverse=True, is_decoy=is_decoy, remove_decoy=False, correction=2)
print('%s\t%s\t%s\t%.3f\t%.1f' % (filename, len(grf), len(f), fdr_n, (float(len(grf)) / len(f) - 1) * 100))
