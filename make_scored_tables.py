#!/usr/bin/env python

import sys
import re

vcf, table, qtl = sys.argv[1:]

var = {}
with open(vcf, 'r') as vc:
    for line in vc:
        if line.startswith('#'):
            continue
#        print line.strip()
        content = line.strip().split('\t')
        keyval = 'chr' + content[0] + ':' + content[1]
        alleles = content[3] + '>' + content[4]
        rsid = content[2]
        gene = content[7].split('ANN=')[1].split('|')[4] if ('ANN=' in line and len(content[7].split('|')) > 6) else '-'
        efftype = content[7].split('ANN=')[1].split('|')[1] if 'ANN=' in line else '-'
        prov_pred = re.findall('PROVEAN_pred=([^,^;]*)', content[7])[0] if 'PROVEAN_pred' in content[7] else '-'
        pph_pred = re.findall('Polyphen2_HVAR_pred=([^,^;]*)', content[7])[0] if 'Polyphen2_HVAR' in content[7] else '-'
        sift_pred = re.findall('SIFT_pred=([^,^;]*)', content[7])[0] if 'SIFT_pred' in content[7] else '-'
        ExAc_af = re.findall('ExAC_AF=.*?(\d[\.\de\+-]*)', line)[0] if 'ExAC' in line else '0'
        our_af = re.findall('BIOBANK_AF_v20171101=(\d\.[\de\+-]+)', line)[0] if 'BIOBANK_AF_v20171101=' in line else '0'
        var[keyval] = (alleles, rsid, gene, efftype, prov_pred, pph_pred, sift_pred, ExAc_af, our_af)

with open(table, 'r') as tab:
    for line in tab:
        if qtl == 'q':
            if line.startswith('VAR'):
                print '\t'.join(['Location', 'Alleles', 'rsID', 'Gene', 'Effect', 'Predictions', 'ExAC', 'Biobank', 'P', 'Beta'])
                continue
            if not line.startswith('chr'):
                continue
            content = line.strip().split('\t')
            vardata = var[content[0]] if '.' not in content[0] else var[re.findall('chr.*?:\d+', content[0])[0]]
            outlist = [content[0], vardata[0], vardata[1], vardata[2], vardata[3], '/'.join(vardata[4:7]), vardata[7], vardata[8],
                       content[8], content[5]]
            print '\t'.join(outlist)
        else:
            if line.startswith('VAR'):
                print '\t'.join(['Location', 'Alleles', 'rsID', 'Gene', 'Effect', 'Predictions', 'ExAC', 'Biobank', 'REFA', 'HETA', 'HOMA', 'REFU', 'HETU', 'HOMU', 'P', 'OR/beta', 'S1D', 'S1R', 'S2D', 'S2R', 'PDOM', 'PREC'])
                continue
            if not line.startswith('chr'):
                continue
            content = line.strip().split('\t')
            vardata = var[content[0]] if '.' not in content[0] else var[re.findall('chr.*?:\d+', content[0])[0]]
            s1d = str(10 * (int(content[14]) + int(content[15])) - 50 * (int(content[17]) + int(content[18])))
            s1r = str(10 * (int(content[15])) - 50 * (int(content[18])))
            s2d = str(10 * (int(content[17]) + int(content[18])) - 50 * (int(content[14]) + int(content[15])))
            s2r = str(10 * (int(content[18])) - 50 * (int(content[15])))
            outlist = [content[0], vardata[0], vardata[1], vardata[2], vardata[3], '/'.join(vardata[4:7]), vardata[7], vardata[8],
                       '\t'.join(content[13:21]), s1d, s1r, s2d, s2r, content[21], content[23]]
            print '\t'.join(outlist)
