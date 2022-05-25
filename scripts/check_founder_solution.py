#!/usr/bin/env python3

#
# std import
#
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, maxsize
import csv
import re

csv.field_size_limit(maxsize)


def edge_set(data):

    res = set()

    _swap_ = lambda x: x == '-' and '+' or '-'
    for row in csv.reader(data, delimiter='\t'):
        if row[0] == 'L':
            if (row[2] == '-' and row[4] == '-') or \
                    (row[2] != row[4] and int(row[1]) > int(row[3])):
                res.add((int(row[3]), _swap_(row[4]), int(row[1]),
                    _swap_(row[2])))
            else:
                res.add((int(row[1]), row[2], int(row[3]), row[4]))

    return res

def read_founders(data):

    res = set()

    _trans_ = lambda x: x == '>' and '+' or '-'
    _swap_ = lambda x: x == '-' and '+' or '-'

    for row in csv.reader(data, delimiter='\t'):
        path = re.findall('[<>][0-9]+', row[1])

        for i in range(len(path)-1):
            if (path[i][0]  == '<' and path[i+1][0] == '<') or \
                    (path[i][0] != path[i+1][0] and \
                    int(path[i][1:]) > int(path[i+1][1:])):
                res.add((int(path[i+1][1:]), _swap_(_trans_(path[i+1][0])),
                    int(path[i][1:]), _swap_(_trans_(path[i][0]))))
            else:
                res.add((int(path[i][1:]), _trans_(path[i][0]),
                    int(path[i+1][1:]), _trans_(path[i+1][0])))
    return res

if __name__ == '__main__':
    description='''
    Calculates node degree from GFA (V1) and prints it to stdout.
    '''
    parser = ArgumentParser(formatter_class=ADHF, description=description)
    parser.add_argument('gfa_file', type=open, help='Graph in GFA version 1 format')
    parser.add_argument('founder_set', type=open, help='Founder set solution')

    args = parser.parse_args()

    out = stdout

    edges = edge_set(args.gfa_file)
    founders = read_founders(args.founder_set)

    print(edges == founders, file=out)

