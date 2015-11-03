
import argparse
import pandas
import os


from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def read_table(table_f):
    """ Returns a dataframe from reading a table """
    df = pandas.read_csv(table_f, sep="\t", header=0, index_col=0)

    # convert potentially numerical row names into strings
    df.index = [str(i) for i in df.index]

    return df

def save_figure(ax, file):
    ax.get_figure().savefig(file)

def bootstrap_boxplot(df):
    """ Makes a boxplot from the iterations of the iterative kfolds program """

    ax = df.boxplot(return_type="axes")
    return ax

def get_mean(dfs, function="mean"):
    """ Applies gets row means over multiple dfs and returns a single df of results using key names as columns """

    # make a single df that has a column of a summary value of each df
    new_df = pandas.DataFrame()
    for name, df in dfs.items():
        if len(df.columns) == 1:
            summary = df[df.columns[0]]
        else:
            summary = df.mean(axis=1)
        if new_df is None:
            new_df[name] = summary
        else:
            new_df[name] = summary

    return new_df

def plot_by_genome_comparison(df):
    """ Plots by-genome comparison for a set of experiments. Assumes cols are experiments and rows are genomes """

    ax = df.plot(kind="box")
    return ax

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="capable of making a variety of plots. Made to be ran in ipython")
    parser.add_argument("-table", help="tables to be read into pandas", required=True, nargs="+")
    #parser.add_argument("-plot", help="the type of plot to make", required=True)
    #parser.add_argument("-out", help="filename for the figure", default="figure.png")
    
    args = parser.parse_args()

    plot2function = {"bootstrap_boxplot": bootstrap_boxplot
            }


    tables = {}
    for table_f in args.table:
        table_name = os.path.splitext(os.path.basename(table_f))[0]
        
        tables[table_name] = read_table(table_f)
