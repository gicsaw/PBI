#!/usr/bin/env python
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description='csv to pickle')
    parser.add_argument('-i', '--csv_file', type=str, required=True,
                        help='input csv file')
    parser.add_argument('-o', '--pickle_file', type=str, required=True,
                        help='output pickle file')

    args = parser.parse_args()

    csv_file = args.csv_file
    pickle_file = args.pickle_file

    df = pd.read_csv(csv_file)
    df.to_pickle(pickle_file)


if __name__ == "__main__":
    main()
