#!/usr/bin/env python3

""" dnamhash.py

Perform feature hashing on a data matrix.

""" 

import numpy as np
import pandas as pd
import mmh3 
import os

def feature_hash(arr, target_dim=1000):
    """ feature_hash
    
    Main hash function.

    """
    low_d_rep = [0 for _ in range(target_dim)]
    for i, el in enumerate(arr):
        hashed = mmh3.hash(str(i))
        if hashed > 0:
            low_d_rep[hashed % target_dim] += arr[i]
        else:
            low_d_rep[hashed % target_dim] -= arr[i]
    return low_d_rep

if __name__ == "__main__":
    """ run feature hashing from command line

    Apply the feature hashing script to a data file from command line.
    Specify the file name/path, and optionally specify the hash dimensions 
    (default 1000), separator symbol (default " "), and new file contianing the
    hashed features data (default args.filepath.split(".")[0] + "_hashed.csv").

    """
    import argparse
    parser = argparse.ArgumentParser(description='Arguments for dnamhash.py')
    parser.add_argument("--filepath", type=str, 
        required=True, default=None, 
        help='File name or path to file on which to do feature hashing.')
    parser.add_argument("--hashdim", type=int, 
        required=False, default=1000, 
        help='Target dimensions of hashed file.')
    parser.add_argument("--sepsymbol", type=str, 
        required=False, default=" ", 
        help='Separation symbol for file to read, defaults to single space (" ").')
    parser.add_argument("--newpath", type=str, 
        required=False, default=None, 
        help='File name or path to new hashed file.')
    if args.filepath:
        fn = args.filepath
        if os.path.exists(fn):
            if not args.newpath:
                hp = args.filepath.split(".")[0] + "_hashed.csv"
            else:
                hp = agrs.newpath
            datl = []
            nh = args.hashdim
            # write first line as colnames
            print("Writing colnames...")
            with open(fn, "r") as fncon:
                cn1 = fncon.readline()[0]
                with open(nh, "w") as nhcon:
                    cnnew = ["hashedfeature" + str(i) for i in range(fh+1)[1::]]
                    cnwrite = cn1 + cnnew
                    nhcon.write(cnwrite)
            print("Writing new hashed feature data...")
            with open(nh, "a") as nhcon:
                for li, line in enumerate(ff):
                    if li > 0:
                        cdat1 = line[0] # row label
                        rdati = line.split(args.sepsymbol)[1::] # row data
                        lmed = np.median([float(li) for li in ll 
                                if not li=="NA"]
                            )
                        # impute missing values using medians
                        ldat = [float(li) if not li=="NA" 
                            else lmed for li in ll]
                        # new hashed data
                        lfh = feature_hash(ldat, target_dim = nh)
                        newrowi = cdat1 + lfh
                        nhcon.write(','.join(str(x) for x in newrowi))
                        # print(str(li), end="\r")
        else:
            print("Please provide a valid filename or path")

    else:
        print("Please provide a valid filename or path")



