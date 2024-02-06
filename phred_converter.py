#!/usr/bin/env python

import argparse
from math import log10

def phred_to_quantile(phred_score):
    """
    Convert PHRED score to quantile.

    Parameters:
    - phred_score (float): The PHRED score to be converted.

    Returns:
    - float: The quantile value.
    """
    return 10 ** (-phred_score / 10)

def quantile_to_phred(quantile):
    """
    Convert quantile to PHRED score.

    Parameters:
    - quantile (float): The quantile value to be converted.

    Returns:
    - float: The PHRED score.
    """
    return -10 * (log10(quantile))

def main():
    parser = argparse.ArgumentParser(description='Convert PHRED score to quantile or vice versa.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--phred', type=float, help='PHRED score to be converted')
    group.add_argument('-q', '--quantile', type=float, help='Quantile to be converted')

    args = parser.parse_args()

    if args.phred:
        quantile = phred_to_quantile(args.phred)
        print(f"PHRED Score: {args.phred}, Quantile: {quantile}")
    elif args.quantile:
        phred = quantile_to_phred(args.quantile)
        print(f"Quantile: {args.quantile}, PHRED Score: {phred}")

if __name__ == "__main__":
    main()
