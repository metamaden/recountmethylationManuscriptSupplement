#!/usr/bin/env python3

""" dnamhash.py

Authors: Sean Maden, Abhi Nellore

Perform feature hashing on a dataset. The main function `feature_hash` performs
 the hashing (a.k.a. "hashing trick") on an array or row of data, into 
 `target_dim` dimensions. This script can be called from command line. For 
 example, use the following from a shell session for a csv table (with 500 rows 
 X 400000 columns) called "input.csv":

> python3 dnamhash.py --filepath "input.csv" --hashdim 1000 --sepsymbol ","

This creates a new csv of the hashed feature data called "input_hashed.csv", 
which has 500 rows X 1000 columns.

""" 

import numpy as np
import mmh3 
import os

def feature_hash(arr, target_dim=1000):
    """ feature_hash
    
    Performs feature hashing on an array data object. This converts a table of
    M X N1 (samples X features) to a table of M X N2 (samples X hashed features), 
    where N1 is the original column count/feature dimensions, and N2 is the new 
    column count/dimensions of the hashed dataset, and typically N1 >> N2.

    Arguments:
    * arr : A single row of numeric sample data as an array, where the array 
        length equals the original feature count.
    * target_dim : The target dimensions of the new row of feature hashed data 
        (integer, 1000).

    Returns: The hashed data row.

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
    (default 1000), separator symbol (default " "), and new file containing the
    hashed features data (default <INPUT TABLE NAME> + "_hashed.csv").

    Arguments:

    * --filepath: Path to input table (rows are samples, columns are features).
    * -- hashdim: Target number of columns for the hashed data (integer, 1000).
    * -- sepsymbol: Separator symbol for the input data table (string, " ").
    * -- newpath: Optional file name or path to hashed csv (string, 
    default <input_filename> + "_hashed.csv")

    Returns: Feature hashed dataset as a csv table.

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
    args = parser.parse_args()
    if args.filepath and os.path.exists(args.filepath):
        fn = args.filepath # input filename
        nh = args.hashdim # target hashed dimensions
        if not args.newpath:
            hp = args.filepath.split(".")[0] + "_hashed.csv"
        else:
            hp = args.newpath
        print("results will be written to file: '" + hp + "' ...")
        with open(fn, "r") as fncon:
            with open(hp, "w") as hpcon:
                print("writing colnames...")
                lread = fncon.readline()
                lform = lread.split(args.sepsymbol)
                cn1 = lform[0]
                cnnew = ["hashedfeature_" + str(i) for i in range(nh+1)[1::]]
                rowwrite = [cn1] + cnnew
                newrow = ','.join(str(x) for x in rowwrite) + "\n"
                hpcon.write(newrow)
            hpcon.close()
        fncon.close()
        print("writing new hashed feature data...")
        with open(fn, "r") as fncon:
            with open(hp, "a") as hpcon:
                for li, line in enumerate(fncon):
                    if li > 0:
                        lform = line.split(args.sepsymbol)
                        cdat1 = lform[0] # row label
                        rdati = lform[1::] # row data
                        rdati[-1] = rdati[-1].replace("\n", "")
                        # impute NAs then do feature hashing
                        fldat = [float(li) for li in rdati if not li=="NA"]
                        lmed = np.median(fldat)
                        ldat = [float(li) if not li=="NA" else lmed for li in rdati]
                        lfh = feature_hash(ldat, target_dim = nh)
                        rowwrite = [cdat1] + lfh
                        newrow = ','.join(str(x) for x in rowwrite) + "\n"
                        hpcon.write(newrow)
            hpcon.close()
        fncon.close()
        print("finished writing data to new table '" + hp + "'")
    else:
        print("please provide a valid argument for `--filepath`" +
            "to an existing file")