
import pandas
import logging

logging.basicConfig()
LOG = logging.getLogger(__name__)


class TraitTableEntry(object):
    """ A single entry in a trait table """

    def __init__(self, name):
        self.name = name
        self.traits = {}

    def __str__(self):
        return "TraitTableEntry {}".format(self.name)

    def add_trait(self, trait, value):
        """ Checks traits to make sure it doesn't already exist in the dict and adds it """
        if trait in self.traits:
            raise ValueError("{} already has a trait called '{}'.".format(str(self), trait))
        else:
            # see if we can convert the trait into a number
            try:
                value = float(value)
            except ValueError:
                pass

            self.traits[trait] = value

    def correlation(self, other, traits=None):
        """
        Finds the correlation between self and other for the listed traits

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

        corr = df.corr(method="spearman")

        return corr.loc["self", "other"]


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

                for index, val in enumerate(trait_values.split("\t")):
                    tte.add_trait(self.traits[index], val)

                yield tte

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


        Something is wrong with this method and it sometimes returns incorrect results!!!
        """

        to_find = len(subset_names)
        found = 0
        for entry in self:

            # check if we have exhausted the list to speed up the search
            if found == to_find:
                if remove:
                    yield entry
                else:
                    return

            if entry.name in subset_names:

                if not remove:
                    found += 1
                    yield entry
            else:
                if remove:
                    found += 1
                    yield entry

    @staticmethod
    def write_entry(entry, fh, traits):
        to_write = [entry.name]

        for trait in traits:
            try:
                to_write.append(str(entry.traits[trait]))
            except KeyError:
                LOG.warning("Entry {} doesn't have trait {}. Writting 'NA'".format(str(entry), trait))

