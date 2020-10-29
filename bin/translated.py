#!/usr/bin/env  python3

from pathlib import Path
import argparse
import sys
import csv
import os
import wget


def main():

    antibiotics_code = "bin/antibiotics_code.v3.tsv"
    argannot = Path(os.path.join(
                    os.path.expanduser('~'), antibiotics_code))
    try:
        argannot.resolve()
    except FileNotFoundError:
        url = 'https://git.io/fjtPO'
        wget.download(url, str(argannot))

    code = {}
    with open(str(argannot), 'r', encoding='utf-8') as code_table:
        for gene in csv.reader(code_table, delimiter='\t'):
            code[gene[1]] = gene[0]

    srst_results = {}
    with open(sys.argv[1], 'r') as argannotable:
        for k in csv.reader(argannotable, delimiter='\t'):
            srst_results[k[0]] = k[1:]

    result = []
    for id in sorted(srst_results):
        assignid = {}
        if srst_results[id]:
            for gene in srst_results[id]:
                for gene_id in sorted(code):
                    if gene == gene_id:
                        assignid[id] = code[gene_id]
                result.append(assignid.copy())
        else:
            assignid[id] = "NP"
            result.append(assignid.copy())

    super_dict = {}
    for x in result:
        for k, v in x.items():
            super_dict.setdefault(k, []).append(v)

    with open(sys.argv[2], 'a', encoding='utf-8') as output:
        for k, v in sorted(super_dict.items()):
            if args.freq:
                print(k + '\t' + ', '.join(v), file=output)
            else:
                print(k + '\t' + ', '.join(sorted(set(v))), file=output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage='translate.py srst2_argannot.tsv output_translated.tsv',
                                     description="This script translate the srst2_argannot.tsv "
                                                 "tableto an antibiotic category")
    parser.add_argument("[1] srst2_argannot.tsv",
                        help="srst2_argannot.tsv file")
    parser.add_argument("[2] output_translated.tsv",
                        help="output file translated")
    parser.add_argument("-freq", action="store_true",
                        help="Translate without uniqueness")

    args = parser.parse_args()
    main()

    if len(sys.argv) == 1:
        parser.print_help()
    sys.exit(1)
