
import argparse
import pandas
from puppetcrust import trait_table

parser = argparse.ArgumentParser(description="Compares two trait tables using a given metric")
parser.add_argument("-tab1", help="the trait table to use for table 1. This is the primary table", required=True)
parser.add_argument("-tab2", help="the table to compare to", required=True)
parser.add_argument("-to_compare", help="a file with a list of names to compare", default=None)
parser.add_argument("-metric", help="the metric to use [%(default)s]", choices=["disimilarity", "spearman"], default="disimilarity")
parser.add_argument("-out", help="file to write the results [%(default)s]", default="compare_two_tables_output.tab")

args = parser.parse_args()

if args.to_compare:
    to_compare = []
    with open(args.to_compare, 'r') as IN:
        for line in IN:
            to_compare.append(line.strip())
else:
    to_compare = None

results = trait_table.TraitTableManager.compare_two_tables(args.tab1, args.tab2, to_compare, args.metric)

df = pandas.DataFrame.from_dict(results, orient="index")

df.to_csv(args.out, sep="\t", index_label="genomes")
    
