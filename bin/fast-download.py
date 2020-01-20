#!/usr/bin/env python3

import argparse
import logging
import subprocess
import sys


def main():

    if args.ssh_key == 'linux':
        ssh_key_file = '$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh'
    elif args.ssh_key == 'osx':
        app = r"$HOME/Applications/Aspera Connect.app"
        ssh_key_file = '%s/Contents/Resources/asperaweb_id_dsa.openssh' % (app)
    else:
        ssh_key_file = args.ssh_key
    logging.info("Using aspera ssh key file: {}".format(ssh_key_file))

    run_id = args.run_identifier
    output_directory = args.output_directory

    logging.info("Querying ENA for FTP paths for {}..".format(run_id))
    report = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession="
    option = "&result=read_run&fields=fastq_ftp&download=txt"
    # run_id = "PRJNA385215"
    text = subprocess.check_output(
                "curl --silent '{}{}{}'".format(report, run_id, option),
                shell=True)

    ftp_urls = []
    header = True
    for line in text.decode('utf8').split('\n'):
        if header:
            header = False
        else:
            for url in line.split(';'):
                if url.strip() != '':
                    ftp_urls.append(url.strip())
    if len(ftp_urls) == 0:
        logging.warn("No FTP download URLs found \
                      for run {}, cannot continue".format(run_id))
        sys.exit(1)
    else:
        logging.info("Found {} FTP URLs for \
                      download e.g. {}".format(len(ftp_urls), ftp_urls[0]))

    for url in ftp_urls:
        cmd = "ascp -QT -l 300m -P33001 {} -i {} \
               era-fasp@fasp.sra.ebi.ac.uk:{} {}".format(
            args.ascp_args,
            ssh_key_file,
            url.replace('ftp.sra.ebi.ac.uk', ''),
            output_directory)
        logging.info("Running command: {}".format(cmd))
        subprocess.check_call(cmd, shell=True)

    logging.info("DONE")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
                      description="Download FASTQ files from the \
                                   European Nucleotide Archive (ENA) \
                                   using Aspera Connect client.\n\n"
                                   "Requires curl and ascp in $PATH")
    parser.add_argument("run_identifier", help="Run number to download")
    parser.add_argument("--output-directory", "--output_directory",
                        help="Output directory [default: '.']", default=".")
    parser.add_argument("--ssh-key", "--ssh_key",
                        help="'linux' or 'osx' paths [default: 'linux']",
                        default="linux")
    parser.add_argument("--ascp-args", "--ascp_args",
                        help="extra arguments to pass to ascp [default: '']",
                        default="")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s %(levelname)s: %(message)s",
                        datefmt="%m/%d/%Y %I:%M:%S %p")
    main()
