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
    # example run
    # python -m RecurrenceRelationSolver.RecurrenceRelationSolver -i ./exampleInOutput/ -o ./output -c 50 -p 100 -q

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
    argParser.add_argument('-c', '--check', type=int,
                           dest='check', required=False,
                           help='How many numbers to verify for the solved recurrences with the recurrence relation. Defaults to 0')
    argParser.add_argument('-p', '--precision', type=int,
                           dest='precision', required=False,
                           help='The amount of places after the decimal point that have to be equal between a test ' +
                                'of the solved equation vs the recurrence relation to be considered correct. Defaults to 4')

    args = argParser.parse_args()
    args.outputdir = args.outputdir if args.outputdir else args.inputdir
    args.check = args.check if args.check else 0
    args.precision = args.precision if args.precision else 4

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

    tolerance = 10**(-args.precision)
    for fn, r in relations.items(): 
        print("Solving %s" % fn)
        try:
            r.solve()
        except Exception as e:
            logging.error("Exception occured while solving recurrence: %s" % r.getRecurrence())
            logging.error(e, exc_info = True)
            continue

        # Verify the solved result 
        failed = False
        start = r.getLowerBoundDomain()
        for i in range(start, start + args.check):
            iterative_result = r.calculateValueFromRecurrence(i)
            solved_result = r.calculateValueFromSolved(i)
            if abs(iterative_result - solved_result) >= tolerance:
                failed = True
                logging.error("Verification of solved recurrence failed at n = %d for relation: %s" % (i, r.getRecurrence()))
                logging.error("Recurrence says: %s" % str(iterative_result))
                logging.error("Solved says: %s" % str(solved_result))
                logging.error("Delta: %s" % str(abs(iterative_result - solved_result)))
                break
        
        if failed:
            continue

        path = os.path.join(args.outputdir, fn.replace(".txt", "-dir.txt"))
        with open(path, "w+") as f:
            f.write("sdir := n -> %s;\n" % r.solve())


if __name__ == '__main__':
    main()
