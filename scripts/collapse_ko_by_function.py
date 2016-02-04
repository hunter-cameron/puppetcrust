
import argparse
import pandas
import logging

logging.basicConfig()
LOG = logging.getLogger(__name__)

def map_ko_to_function(ko_metadata_f, level=2):
    """
    Makes a dict linking KO to pathway data.

    Level 1 is the top level.
    Level 2 is an intermediate level (Corresponding approx to COGs)
    Level 3 is the pathway level

    KO file should look like this:
    header
    ko\tdescription\tlvl1;lvl2;lvl3|lvl1;lvl2;lvl3|......

    """
    if level not in [1, 2, 3]:
        raise ValueError("Level must be 1, 2, or 3.")

    ko_to_functional = {}
    with open(ko_metadata_f, 'r') as IN:
        # skip header line
        IN.readline()

        for line in IN:
            ko_name, ko_description, ko_pathways = line.rstrip().split("\t")

            # multiple pathways sep by "|"
            for pathway in ko_pathways.split("|"):
                levels = pathway.split(";")

                try:
                    ko_to_functional[ko_name].append(";".join(levels[:level]))
                except KeyError:
                    ko_to_functional[ko_name] = [(";".join(levels[:level]))]
                except IndexError:
                    LOG.warning("{} did not have a pathway at the requested level.".format(ko_name))

    return ko_to_functional


def collapse_kos(table_f, ko_to_functional, orient, out_f):
    
    collapsed_counts = {}

    df = pandas.DataFrame.from_csv(table_f, sep="\t", header=0, index_col=0)


    if orient == "cols":
        df = df.transpose()
   
    new_df = df
    metadata_path = []
    extra_rows = []     # metadata for the extra rows added at the end for KOs with multiple pathways
    for indx in list(df.index):
        try:
            paths = ko_to_functional[indx]
        except KeyError:
            if indx[0] == "K":
                LOG.warning("{} was not present in KO metadata file, possibly because the pathway is unknown. Adding as 'unknown'".format(indx))
                metadata_path.append("unknown")
            else:
                LOG.warning("{} was not present in KO metadata file and doesn't appear to be a KO. It will be omitted from analysis.".format(indx))
                metadata_path.append("omit")

            continue

        # add the metadata
        if len(paths) > 1:
            # add a new duplicate row for each additional path
            for path in paths[1:]:
                new_df = new_df.append(df.loc[indx, :], ignore_index=True)
                extra_rows.append(path)
            
        metadata_path.append(paths[0])

    new_df["ko_pathway"] = metadata_path + extra_rows

    new_df = new_df.groupby("ko_pathway", axis=0).sum() 
    try:
        new_df.drop("omit", axis=0, inplace=True)
    except ValueError:
        pass
    
    if orient == "cols":
        new_df = new_df.transpose()

    new_df.to_csv(out_f, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is a fairly inefficient program because I had forgotten that each KO could be assigned to multiple pathways. The correction for that involved adding duplicate rows to the dataframe to represent each categorization. However, this is a really computationally expensive way of doing things")
    parser.add_argument("-table", help="a table with kos in either the rows or columns", required=True)
    parser.add_argument("-ko_metadata", help="a file with information about each KO", required=True)
    parser.add_argument("-orient", help="the orientation the KOs are in", choices=["rows", "cols"], default="rows")
    parser.add_argument("-level", help="the KO level to use [%(default)s]", choices=[1, 2, 3], type=int, default=2)
    parser.add_argument("-out", help="path to write the new table", default="ko_collapsed.tab")


    args = parser.parse_args()

    ko_to_function = map_ko_to_function(args.ko_metadata, args.level)
    collapse_kos(args.table, ko_to_function, args.orient, args.out)
