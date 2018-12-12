#!/usr/bin/env python3
# coding=utf-8
import argparse
import logging
import glob
import os.path
import os
import errno

from . import RecurrenceRelationParser


def main():
    argParser = argparse.ArgumentParser(
        description=('Solve recurrent relation into closed-form solution'),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    argParser.add_argument('-i', '--inputdir', type=str,
                           dest='inputdir', required=True,
                           help='Input directory where files with recurrence relations are placed.')
    argParser.add_argument('-o', '--outputdir', type=str,
                           dest='outputdir', required=False,
                           help='Output directory where results are saved. Defaults to input directory')
    argParser.add_argument('-q', '--quiet', action='store_true',
                           dest='quiet', help='Only print warnings and errors.')

    args = argParser.parse_args()
    args.outputdir = args.outputdir if args.outputdir else args.inputdir

    loglevel = logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(format='%(message)s', level=loglevel)

    recurrenceParser = RecurrenceRelationParser()

    # Create output directory if it doesn't exist
    try:
        os.makedirs(args.outputdir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    # Read in the relations from the input directory
    relations = {}
    for path in glob.glob(os.path.join(args.inputdir, "comass[0-9][0-9].txt")):
        _, fn = os.path.split(path)
        with open(path, "r") as f:
            relations[fn] = recurrenceParser.parse_recurrence(f.read())


    # Maybe check equation here?

    # Write solved relations to output directory
    for fn, r in relations.items():
        path = os.path.join(args.outputdir, fn.replace(".txt", "-dir.txt"))
        with open(path, "w+") as f:
            f.write("sdir := n -> %s;\n" % r.solve())


if __name__ == '__main__':
    main()
