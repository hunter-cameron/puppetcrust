
import argparse
import pandas

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="capable of making a variety of plots")
    parser.add_argument("-table", help="table to be read into pandas", required=True)
    parser.add_argument("-plot", help="the type of plot to make", required=True)
    parser.add_argument("-out", help="filename for the figure", default="figure.png")
    
    args = parser.parse_args()

    plot2function = {"bootstrap_boxplot": bootstrap_boxplot
            }



    df = read_table(args.table)

    try:
        ax = plot2function[args.plot](df)
    except KeyError:
        raise ValueError("'{}' is not a valid plot type.".format(args.plot))

    save_figure(ax, args.out)
