#!/usr/bin/python3

from pathlib import Path
import argparse
import sys
import csv
import os
import wget


def main():

    parser = argparse.ArgumentParser(
        usage='translate.py srst2_argannot.tsv output_translated.tsv',
        description="This script translate the srst2_argannot.tsv table \
                     to an antibiotic category")
    parser.add_argument(
                    "[1] srst2_argannot.tsv", help="srst2_argannot.tsv file")
    parser.add_argument(
                    "[2] output_translated.tsv", help="output file translated")
    parser.add_argument("-freq", action="store_true",
                        help="Translate without uniqueness ")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    argannot = Path(os.path.join(
                    os.path.expanduser('~'), "bin/antibiotics_code.tsv"))
    try:
        arg_path = argannot.resolve()
        # print(arg_path)
    except FileNotFoundError:
        url = 'https://git.io/vpNpz'
        wget.download(url, str(argannot))

    code = {}
    with open(str(argannot), 'r', encoding='utf-8') as code_table:
        for gene in csv.reader(code_table, delimiter='\t'):
            code[gene[1]] = gene[0]
    # print(code)

    '''
    for gene in sorted(code):
       print(gene,':',code[gene])
    '''

    srst_results = {}
    with open(sys.argv[1], 'r') as argannotable:
        for k in csv.reader(argannotable, delimiter='\t'):
            srst_results[k[0]] = k[1:]
    # print(srst_results)

    '''
    for id in sorted(srst_results):
        if srst_results[id]:
            print(id,':',srst_results[id])
        else:
            print(id)
    '''

    result = []
    for id in sorted(srst_results):
        assignid = {}
        if srst_results[id]:
            for gene in srst_results[id]:
                for gene_id in sorted(code):
                    if gene == gene_id:
                        assignid[id] = code[gene_id]
                # print(assignid)
                result.append(assignid.copy())
        else:
            assignid[id] = "NP"
            # print(assignid)
            result.append(assignid.copy())

    # print(result)

    super_dict = {}
    for x in result:
        for k, v in x.items():
            super_dict.setdefault(k, []).append(v)

    # print(super_dict)

    with open(sys.argv[2], 'a', encoding='utf-8') as output:
        for k, v in sorted(super_dict.items()):
            if args.freq:
                print(k + '\t' + ', '.join(v), file=output)
            else:
                print(k + '\t' + ', '.join(sorted(set(v))), file=output)


if __name__ == '__main__':
    main()
