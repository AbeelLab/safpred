#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
from utils import seq_utils


def main():
    parser = argparse.ArgumentParser(description = 'Download SwissProt database and preprocess it')
    parser.add_argument('--output_dir', '-o', type=str, required=True, 
                        help = 'Directory to write the output files')
    parser.add_argument('--release', '-r', type=str, required=True, 
                        help = 'Database release date to download -- only the "current" option is supported at the moment')

    args = parser.parse_args()
    output_dir = args.output_dir
    release = args.release
    
    if not os.path.exists(output_dir):
        print("Creating the output directory to download the files in")
        os.mkdir(output_dir)
    if not os.path.exists(os.path.join(output_dir, 'uniprot')):
        os.mkdir(os.path.join(output_dir, 'uniprot'))
    if not os.path.exists(os.path.join(output_dir, 'go')):
        os.mkdir(os.path.join(output_dir, 'go'))    
    
    print("Downloading the SwissProt datbase")
    seq_utils.download_sprot(output_dir, release=release)

if __name__ == '__main__':
    main()
