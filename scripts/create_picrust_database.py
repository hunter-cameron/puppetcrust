
import argparse

from puppetcrust import database

parser = argparse.ArgumentParser(description="Filters a PICRUSt table to match a marker FASTA file. Can accept multiple tables or FASTA files and optionally selectively include/exclude genome names. Outputs a new FASTA file and trait table beginning with the string supplied to -prefix")
parser.add_argument("-fasta", help="one or more marker FASTA files. Header should be genome name", nargs="+", required=True)
parser.add_argument("-traits", help="one or more trait tables. First column should be genome names. Tables should all have the same traits", nargs="+", required=True)
parser.add_argument("-prefix", help="prefixx for the new table and markers file. [%(default)s]", default="filtered")
parser.add_argument("-subset", help="file of targeted genome names (default = to keep)")
parser.add_argument("-inverse", help="remove the targeted sequences", action="store_true")


args = parser.parse_args()

dbm = database.DatabaseManager()

# add the data files to the manager
for f in args.fasta:
    dbm.add_fasta(f)

for t in args.traits:
    dbm.add_trait_table(t)

# read in a list of genome names if given
if args.subset:
    subset = []
    with open(args.subset, 'r') as IN:
        for line in IN:
            subset.append(line.strip())
else:
    subset = []

dbm.generate_database(args.prefix, subset=subset, inverse=args.inverse, verbose=True)
