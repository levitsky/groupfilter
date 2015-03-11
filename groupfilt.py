from pyteomics.auxiliary import filter, fdr, qvalues
from csv import DictReader, DictWriter
from itertools import groupby
from sys import argv

def groupfilter(psms, key, score, fdr, is_decoy, reverse=False, remove_decoy=False):
    out, bad, gd = [], [], []
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
            goodl = filter(g, key=score, fdr=fdr, is_decoy=is_decoy,
                    reverse=reverse, remove_decoy=remove_decoy, correction=2, formula=1)
            if len(goodl) == len(g):
                gd.extend(g)
            elif len(goodl):
                out.extend(goodl)
            elif not goodl:
                bad.extend(g)
    if bad:
        b = filter(bad, key=score, fdr=fdr, is_decoy=is_decoy, reverse=reverse,
                remove_decoy=remove_decoy, correction=2, formula=1)
        out.extend(b)
    if gd:
        b = filter(gd, key=score, fdr=fdr, is_decoy=is_decoy, reverse=reverse,
                remove_decoy=remove_decoy, correction=2, formula=1)
        out.extend(b)
    return out

def iterate_fdr(fdr_v, psms, key, score, is_decoy, reverse=True):
    fdr_c = -1
    fdr_u = 0
    lf = 0.
    rf = 1.
    while fdr_c > fdr_v or abs(fdr_c - fdr_v) >= fdr_v * 0.1:
        if fdr_v > fdr_c:
            lf = fdr_u
            fdr_u += (rf - fdr_u) / 2
        else:
            rf = fdr_u
            fdr_u -= (fdr_u - lf) / 2
        grf = groupfilter(psms, key, score, fdr_u, is_decoy, reverse=True)
        if grf:
            fdr_c = fdr(grf, is_decoy=is_decoy, formula=1)
        else:
            fdr_c = 0
    return fdr_u

def write_csv_psms(out, fname, fieldnames):
    out.sort(key = lambda x: x['Q-Value (%)'])
    with open(fname, 'w') as f:
        csvwriter = DictWriter(f, fieldnames=fieldnames, delimiter='\t')#, extrasaction='ignore')
        csvwriter.writerow(dict((fn, fn) for fn in fieldnames))
        for psm in out:
            csvwriter.writerow(psm)


def filt(psms, key, score, is_decoy, fname, fieldnames):
    fdr_c = iterate_fdr(0.01, psms, key, score, is_decoy)
    grf = groupfilter(psms, key, score, fdr_c, is_decoy, reverse=True)
    q_v = qvalues(grf, key=score, is_decoy=is_decoy, reverse=True, remove_decoy=False, formula=1)
    for idx, z in enumerate(grf):
        z['Q-Value (%)'] = q_v[idx][-1] * 100
        z['Cumulative Target'] = 0 #TODO
        z['Cumulative Decoy'] = 0 #TODO
    fdr_n = fdr(grf, is_decoy=is_decoy, formula=1)
    f = filter(psms, fdr=fdr_n, key=score, reverse=True, is_decoy=is_decoy, remove_decoy=False, correction=2, formula=1)
    write_csv_psms(grf, fname, fieldnames)
    print('{}\t{}\t{}\t{:.3f}\t{:.1f}'.format(filename, len(grf), len(f), fdr_n,
        (float(len(grf)) / len(f) - 1) * 100))

if __name__ == '__main__':
    filename = argv[1]
    oname = filename.split('.')[0] + '_groupfilter'
    with open(filename) as f:
        dcr = DictReader(f, delimiter='\t')
        fieldnames = dcr.fieldnames
        psms = list(dcr)
    key = lambda x: (int(float(x['Morpheus Score'])))
    score = lambda x: float(x['Morpheus Score'])
    is_decoy = lambda x: x['Decoy?'] == 'True' and x['Target?'] == 'False'

    filt(psms, key, score, is_decoy, fname=oname + '.PSMs.tsv', fieldnames=fieldnames)
    seq_added = set()
    peptides = []
    for psm in psms:
        if psm['Base Peptide Sequence'] not in seq_added:
            seq_added.add(psm['Base Peptide Sequence'])
            peptides.append(psm)
    filt(peptides, key, score, is_decoy, fname=oname + '.unique_peptides.tsv', fieldnames=fieldnames)

