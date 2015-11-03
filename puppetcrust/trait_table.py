
from __future__ import division


import pandas
import logging

logging.basicConfig()
LOG = logging.getLogger(__name__)


class TraitTableEntry(object):
    """ A single entry in a trait table """

    def __init__(self, name):
        self.name = name
        self.traits = {}
        self.metadata = {}

    def __str__(self):
        return "TraitTableEntry {}".format(self.name)

    def add_trait(self, trait, value):
        """ Checks traits to make sure it doesn't already exist in the dict and adds it """
        # see if we can convert the trait into a number
        try:
            value = float(value)
        except ValueError:
            pass

        # check if the trait is metadata
        if trait.startswith("metadata_"):
            name = trait.split("metadata_")[1]
                
            # check if metadata already present
            if name in self.metadata:
                raise ValueError("{} already has metadata called '{}'.".format(str(self), name))
            else:
                self.metadata[name] = value
        else:
            # check if trait already present
            if trait in self.traits:
                raise ValueError("{} already has a trait called '{}'.".format(str(self), trait))
            else:
                self.traits[trait] = value

    def write(self, fh, traits):
        to_write = [self.name]

        for trait in traits:
            try:
                to_write.append(str(self.traits[trait]))
            except KeyError:
                LOG.warning("Entry {} doesn't have trait {}. Writting 'NA'".format(str(self), trait))

        fh.write("\t".join(to_write) + "\n")

    def compare(self, other, metric, traits=None):
        """
        Compares self to other using the supplied metric.

        If traits is not suppiled, uses all the traits from self.

        Only uses traits that both have. I'm not sure if this is the intuitive default behavior.
        It might make sense to throw an error if a trait isn't found.
        Or to return a list of the ultimate traits used?

        It may be unintuitive because if I pass the same set of traits to multiple comparisons,
        each comparison may actually use different traits and not tell me.
        """

        if traits is None:
            traits = list(self.traits.keys())

        pandas_dict = {}
        for trait in traits:
            try:
                st = self.traits[trait]
            except KeyError:
                continue

            try:
                ot = other.traits[trait]
            except KeyError:
                continue


            pandas_dict[trait] = {"self": st, "other": ot}

        if not pandas_dict:
            raise ValueError("No traits were shared between the entries.")


        df = pandas.DataFrame.from_dict(pandas_dict, orient="index")


        if metric == "spearman":
            corr = df.corr(method="spearman")
            return corr.loc["self", "other"]

        elif metric == "disimilarity":
            # calcs squared Euclidian distance / num traits
            # I picked this distance metric becuase it weights predictions that are more wrong higher

            diff = df["self"] - df["other"]
            diff = diff.abs()
            diff = diff ** 2
            sq_euc = diff.sum()
            return sq_euc / len(df.index)

        elif metric == "positivepred":
            # calcs positive predictive value, which ends up being correct calls / total calls
            diff = df["self"] - df["other"]

            correct = df.count(0)
            return correct / len(diff)


        else:
            raise ValueError("metric '{}' is invalid.".format(metric))

class TraitTableManager(object):
    """ A class for parsing and manipulating trait tables """

    def __init__(self, trait_table_f):
        self.trait_table_f = trait_table_f

        # get headers
        with open(self.trait_table_f, 'r') as IN:
            headers = IN.readline().rstrip().split("\t")

            # set the entry header replacing a comment line if present
            self.entry_header = headers[0].replace("#", "")
            self.traits = headers[1:]

    def __iter__(self):
        """ Yields a TraitTableEntry for each line in a trait table """

        with open(self.trait_table_f, 'r') as IN:
            # skip header line
            IN.readline()

            for line in IN:
                # skip blank lines
                if not line:
                    continue

                try:
                    name, trait_values = line.rstrip().split("\t", 1)
                except ValueError:
                    print((line,))

                tte = TraitTableEntry(name)

                # add all traits
                for index, val in enumerate(trait_values.split("\t")):
                    tte.add_trait(self.traits[index], val)

                yield tte

    @classmethod
    def compare_two_tables(cls, tab1, tab2, to_compare=None, metric="disimilarity"):
        """ This is a convenience method that compares the entries in the list to_compare (all from tab1 if is None) from the two trait tables. Returns a dict indexed by the names in to_compare """

        # Open up a manager for each table.
        # should do a simple check here to allow users to sumbit ttm's as well as files
        ttm1 = cls(tab1)
        ttm2 = cls(tab2)

        # get a list of entries from the first table for quick comparisons
        if to_compare is None:
            LOG.info("Using all the entries from table 1 for the comparision.")
            comp_entries = [entry for entry in ttm1]
        else:

            comp_entries = [entry for entry in ttm1 if entry.name in to_compare]
            
            # check for names not found in table 1
            entry_names = [entry.name for entry in comp_entries]
            missing = []
            for name in to_compare:
                if name not in entry_names:
                    missing.append(name)

            if missing:
                print("Missing Names:")
                for name in missing:
                    print(name)

                raise ValueError("Some names from to_compare not found in table 1. See stdout for details. Refusing to continue.")
            else:
                LOG.info("Found all names from to_compare in table 1.")

        results = {entry.name: None for entry in comp_entries}
        for comp in comp_entries:
            for comp2 in ttm2:
                if comp.name == comp2.name:
                    results[comp.name] = {metric: comp.compare(comp2, metric=metric)}

                    # try to get a NSTI value
                    try:
                        results[comp.name]["NSTI"] = comp.metadata["NSTI"]
                    except KeyError:
                        try:
                            results[comp.name]["NSTI"] = comp2.metadata["NSTI"]
                        except KeyError:
                            pass


        # make sure all genomes were compared successfully
        for genome in results:
            if results[genome] is None:
                raise ValueError("Calculation failed for genome '{}'".format(genome))

        return results
    
    def get_ordered_traits(self, metadata_last=True):
        """ Returns an ordered list of traits by a natural sort algorithm that optionally sends metadata to the back. """

        def convert(char):
            """ Attempts to convert a character into an integer """
            try:
                return int(char)
            except ValueError:
                return char

        def nat_sort(entry):
            """ Performs a natural sort that will sort text and numbers in a way that makes sense """
            if metadata_last:
                if entry.startswith("metadata"):
                    # I append "~" to the beginning because of its high Unicode value
                    entry = "~~~" + "metadata"

            return [convert(char) for char in entry]


        return sorted(self.traits, key=nat_sort)

    def get_subset(self, subset_names, remove=False):
        """
        A filter around the iter method that only gets entries in the subset_names list (or removes them)
        """

        to_find = len(subset_names)
        found = 0
        for entry in self:


            # check if we have exhausted the list to speed up the search
            if found == to_find:
                if remove:
                    yield entry
                    continue
                else:
                    return

            if entry.name in subset_names:
                found += 1

                if not remove:
                    yield entry
            else:
                if remove:
                    yield entry

    def write_subset(self, path, subset_names, remove=False):
        """
        Write a new trait table including only a subset of the main one
        """

        written = []
        with open(path, 'w') as OUT:
            # write headers
            OUT.write("\t".join(["OTU"] + self.traits) + "\n")
            
            for entry in self.get_subset(subset_names, remove):
                written.append(entry.name)
                entry.write(OUT, self.traits)
            
            return written
