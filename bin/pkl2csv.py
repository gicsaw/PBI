#!/usr/bin/env python
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description='pickle to csv')
    parser.add_argument('-i', '--pickle_file', type=str, required=True,
                        help='input pickle file')
    parser.add_argument('-o', '--csv_file', type=str, required=True,
                        help='output csv file')

    args = parser.parse_args()

    pickle_file = args.pickle_file
    csv_file = args.csv_file

    df = pd.read_pickle(pickle_file)
    df.to_csv(csv_file)


if __name__ == "__main__":
    main()
