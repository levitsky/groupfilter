from pyteomics.auxiliary import filter, fdr, qvalues
from csv import DictReader, DictWriter, writer
from itertools import groupby
from collections import defaultdict
from sys import argv

def groupfilter(psms, key, score, fdr, is_decoy, reverse=False, remove_decoy=False):
    full_out, out, bad, gd = [], [], [], []
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
                    reverse=reverse, remove_decoy=remove_decoy, correction=1, formula=1)
            if len(goodl) == len(g):
                gd.extend(g)
            elif len(goodl):
                q_v = qvalues(g, key=score, is_decoy=is_decoy, reverse=True, remove_decoy=False, formula=1, correction=1)
                for idx, z in enumerate(g):
                    z['Q-Value (%)'] = q_v[idx][-1] * 100
                full_out.extend(g)
                out.extend(goodl)
            elif not goodl:
                bad.extend(g)
    if bad:
        b = filter(bad, key=score, fdr=fdr, is_decoy=is_decoy, reverse=reverse,
                remove_decoy=remove_decoy, correction=1, formula=1)
        q_v = qvalues(bad, key=score, is_decoy=is_decoy, reverse=True, remove_decoy=False, formula=1, correction=1)
        for idx, z in enumerate(bad):
            z['Q-Value (%)'] = q_v[idx][-1] * 100
        full_out.extend(bad)
        out.extend(b)
    if gd:
        q_v = qvalues(gd, key=score, is_decoy=is_decoy, reverse=True, remove_decoy=False, formula=1, correction=1)
        for idx, z in enumerate(gd):
            z['Q-Value (%)'] = q_v[idx][-1] * 100
        full_out.extend(gd)
        b = filter(gd, key=score, fdr=fdr, is_decoy=is_decoy, reverse=reverse,
                remove_decoy=remove_decoy, correction=1, formula=1)
        out.extend(b)
    return out, full_out

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
        grf, psms = groupfilter(psms, key, score, fdr_u, is_decoy, reverse=True)
        if grf:
            fdr_c = fdr(grf, is_decoy=is_decoy, formula=1)
        else:
            fdr_c = 0
    return fdr_u

def write_csv_psms(out, fname, fieldnames):
    out.sort(key=lambda x: x['Q-Value (%)'])
    with open(fname, 'w') as f:
        csvwriter = DictWriter(f, fieldnames=fieldnames, delimiter='\t')#, extrasaction='ignore')
        csvwriter.writerow(dict((fn, fn) for fn in fieldnames))
        for psm in out:
            csvwriter.writerow(psm)


def write_csv_proteins(out, fname): #TODO
    fieldnames = ['protein', 'Q-Value (%)']
    out.sort(key=lambda x: x[1])
    with open(fname, 'w') as f:
        csvwriter = writer(f, delimiter='\t')
        # csvwriter = DictWriter(f, fieldnames=fieldnames, delimiter='\t')#, extrasaction='ignore')
        # csvwriter.writerow(dict((fn, fn) for fn in fieldnames))
        csvwriter.writerow(fieldnames)
        for prot in out:
            csvwriter.writerow(prot)


def filt(psms, key, score, qscore, is_decoy, fname, fieldnames):
    fdr_c = iterate_fdr(0.01, psms, key, score, is_decoy)
    grf, psms = groupfilter(psms, key, score, fdr_c, is_decoy, reverse=True)
    psms.sort(key=qscore)
    q_v = qvalues(psms, key=qscore, is_decoy=is_decoy, reverse=False, remove_decoy=False, formula=1)
    for idx, z in enumerate(psms):
        z['Q-Value (%)'] = q_v[idx][-1] * 100
        z['Cumulative Target'] = 0 #TODO
        z['Cumulative Decoy'] = 0 #TODO
    fdr_n = fdr(grf, is_decoy=is_decoy, formula=1)
    f = filter(psms, fdr=fdr_n, key=score, reverse=True, is_decoy=is_decoy, remove_decoy=False, correction=1, formula=1)
    write_csv_psms(psms, fname, fieldnames)
    print('{}\t{}\t{}\t{:.3f}\t{:.1f}'.format(filename, len(grf), len(f), fdr_n,
        (float(len(grf)) / len(f) - 1) * 100))
    return grf

if __name__ == '__main__':
    filename = argv[1]
    oname = filename.split('.')[0] + '_groupfilter'
    with open(filename) as f:
        dcr = DictReader(f, delimiter='\t')
        fieldnames = dcr.fieldnames
        psms = list(dcr)
    key = lambda x: (int(float(x['Morpheus Score'])))
    score = lambda x: float(x['Morpheus Score'])
    qscore = lambda x: float(x['Q-Value (%)'])
    is_decoy = lambda x: x['Decoy?'] == 'True' and x['Target?'] == 'False'

    filt(psms, key, score, qscore, is_decoy, fname=oname + '.PSMs.tsv', fieldnames=fieldnames)
    seq_added = dict()
    peptides = []
    for psm in psms:
        seq_added[psm['Base Peptide Sequence']] = max(psm['Morpheus Score'], seq_added.get(psm['Base Peptide Sequence'], 0))
    for psm in psms:
        if psm['Base Peptide Sequence'] in seq_added and psm['Morpheus Score'] == seq_added[psm['Base Peptide Sequence']]:
            # seq_added.add(psm['Base Peptide Sequence'])
            peptides.append(psm)
    fpeptides = filt(peptides, key, score, qscore, is_decoy, fname=oname + '.unique_peptides.tsv', fieldnames=fieldnames)

    prot_score = lambda x: x[1]
    prot_is_decoy = lambda x: 'DECOY_' in x[0]
    prots = defaultdict(float)
    for psm in fpeptides:
        prots[psm['Protein Description']] += float(psm['Morpheus Score'])
    prots = sorted(prots.items(), key=prot_score, reverse=True)
    q_v = qvalues(prots, key=prot_score, is_decoy=prot_is_decoy, reverse=True, remove_decoy=False, formula=1)
    prots_q = []
    for idx, z in enumerate(prots):
        prots_q.append((z[0], q_v[idx][-1] * 100))
    prots_fil = filter(prots_q, fdr=0.01, key=prot_score, reverse=False, is_decoy=prot_is_decoy, remove_decoy=False, correction=0, formula=1)
    write_csv_proteins(prots_q, fname=oname + '.protein_groups.tsv')
    print len(prots_fil)
    print fdr(prots_fil, is_decoy=prot_is_decoy, formula=1)
