
import argparse
from Bio import Phylo
import random
import os
import pandas
import logging

from puppetcrust import trait_table, executer

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

logging.basicConfig()
LOG = logging.getLogger(__name__)


def _check_and_mkdir(path):
    if not os.path.isdir(path):
        os.mkdir(path)



class BootStrapper(object):
    """ Bundles the data and methods required for a bootstrapped experiment """
    
    def __init__(self, tree, traits, k, work_dir, to_test=None):
        self.tree_f = tree
        self.traits_f = traits
        self.k = k
        self.work_dir = work_dir
        
        # make work_dir if it doesn't exist
        _check_and_mkdir(self.work_dir)
        
        # read in to_test genomes if necessary
        if to_test is None:
            self.to_test = None
        else:
            self.to_test = []
            with open(to_test, 'r') as IN:
                for line in IN:
                    self.to_test.append(line.strip())
                    

        # init for later use
        self.kfolds = []

    def run(self, bootstrap, metric):

        LOG.info("Beginning bootstraping. Iterations={}".format(str(bootstrap)))

        # make all the partitions and start running PICRUSt
        for job in range(bootstrap):
            job_dir = self.work_dir + "/" + "kfolds" + str(job)
        
            kfolder = KFolder(tree=self.tree_f, traits=self.traits_f, work_dir=job_dir)
            kfolder.make_partitions(k=self.k, to_test=self.to_test)
            self.kfolds.append(kfolder)

        # wait for them all to complete and process the results
        pandas_dict = {}
        for indx, kfolder in enumerate(self.kfolds):
            LOG.info("Processing iteration {}...".format(indx))
            kfolder.analyze(metric=metric)
            pandas_dict["iter" + str(indx)] = kfolder.results
         
        """
            Removed because this should really be done in pandas before analysis because it is derrived data.

            # make a column with the best
            for genome in kfolder.results:
                try:
                    if kfolder.results[genome] > pandas_dict["best"][genome]:
                        pandas_dict["best"][genome]
                except KeyError:
                    pandas_dict["best"][genome] = kfolder.results[genome]
        """

        df = pandas.DataFrame.from_dict(pandas_dict, orient="columns")
        df.to_csv(self.work_dir + "/" + "summary.{}.tab".format(metric), sep="\t", index_label="genome")


class KFolder(object):
    """ Bundles the data and methods required for a single k-fold experiment """
    
    def __init__(self, tree, traits, work_dir):
        self.tree_f = tree
        self.ttm = trait_table.TraitTableManager(traits)
        self.work_dir = work_dir

        # simple check subject to race condition, I should rewrite using a try-except
        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        # attribs to be used later
        self.partitions = []
        self.results = {}

        LOG.info("Initialized K-Folds experiment in directory '{}'.".format(work_dir))


    def make_partitions(self, k, to_test=None, seed=0):
        """ Randomly makes partitions of the traits table """
        
        # try to load partitions
        for group in range(k):
            try:
                partition = Partition.load_partition(self, self.work_dir + "/" + "partition{}".format(group))
                self.partitions.append(partition)
            except ValueError:
                
                # if error at first this is likely a new run
                if group == 0:
                    break
                else:
                    raise ValueError("Successfully loaded some partitions but failed to load others. Refusing to run analysis in this directory. Please delete the output directory or select a new location.")
        else:
            LOG.info("Successfully loaded partitions.")
            return

        # get genomes from the external nodes of the tree
        # comments_are_confidence allows me to keep the digit only node names
        # this parser tries to conv names to confidence otherwise
        tree = Phylo.read(self.tree_f, "newick", comments_are_confidence=True)
        genomes = [node.name for node in tree.get_terminals()]
        
        # make sure all the genomes in to_test are in the tree and then use only to test genomes for the rest
        if to_test:
            for genome in to_test:
                if not genome in genomes:
                    raise ValueError("Genome in test set '{}' is not in the tree. Aborting.".format(genome))
            genomes = to_test


        # calculate the number of genomes per group
        genomes_per_group = int(len(genomes) / k)
        remainder = len(genomes) % k


        # make random groups
        for group in range(k):
            # add remainder to the first groups
            if remainder:
                remainder -= 1
                count = genomes_per_group + 1
            else:
                count = genomes_per_group

            sample = random.sample(genomes, count)

            # remove the sample from the original
            genomes = [g for g in genomes if g not in sample]

            partition = Partition(sample, self, self.work_dir + "/" + "partition{}".format(str((group))))
            partition.run()
            self.partitions.append(partition)
        
            LOG.info("Created partition {} with n={}.".format(str(group), str(count)))

    def analyze(self, metric):
        # wait for all to finish and then process 
        for p in self.partitions:
            try:
                p.parse_results(metric)
            except:
                LOG.info("Waiting for partition in '{}' to finish.".format(p.work_dir))
                executer.PicrustExecuter.wait_for_job(p.job_name)
                p.parse_results(metric)

            # add the results to the overall results
            self.results.update(p.results)


class Partition(object):
    def __init__(self, genomes, kfolder, work_dir):
        self.genomes = genomes
        self.kfolder = kfolder
        self.work_dir = work_dir

        _check_and_mkdir(work_dir)
        
        # status is one of "incomplete", "finished", or "failed"
        self.status = "incomplete"

        # some paths for convenience
        self.ref_traits_f = self.work_dir + "/" + "reference_traits.tab"
        self.test_genomes_f = self.work_dir + "/" + "test_genomes.txt"

        # attribs to be set later
        self.job_name = None
        self.pred_traits_f = None
        self.results = {g: None for g in self.genomes}

    @classmethod
    def load_partition(cls, kfolder, partition_dir):
        """ Loads a partition from a directory. """
        
        # see if the partition directory even exists
        if os.path.isdir(partition_dir):
            partition = cls(genomes=[], kfolder=kfolder, work_dir=partition_dir)

            # read in the genomes
            if os.path.isfile(partition.test_genomes_f):
                with open(partition.test_genomes_f, 'r') as IN:
                    for line in IN:
                        partition.genomes.append(line.strip())

                partition.results = {g: None for g in partition.genomes}

                # see what other work needs to be done
                if not os.path.isfile(partition.ref_traits_f):
                    partition.write_ref_traits()

                # check if PICRUSt has already run
                if os.path.isfile(partition.work_dir + "/" + "predicted_traits.tab"):
                    partition.pred_traits_f = partition.work_dir + "/" + "predicted_traits.tab"
                else:
                    partition.run_picrust()

                return partition

            else:
                raise ValueError("Test genomes file doesn't exist.")
        else:
            raise ValueError("Partition directory doesn't exist.")

    def run(self):
        self.write_test_genomes()
        self.write_ref_traits()
        self.run_picrust()

    def write_ref_traits(self):
        """ Write all the traits except the partition """

        self.kfolder.ttm.write_subset(self.ref_traits_f, self.genomes, remove=True)

    def write_test_genomes(self):
        """ Write a list of genome names being tested """
        
        with open(self.test_genomes_f, 'w') as OUT:
            for g in self.genomes:
                OUT.write(g + "\n")

    def run_picrust(self):
        self.job_name, self.pred_traits_f = executer.PicrustExecuter.predict_traits_wf(
                tree=self.kfolder.tree_f, trait_table=self.ref_traits_f, limit=self.test_genomes_f, base_dir=self.work_dir)

    def parse_results(self, metric):
        obs_ttm = self.kfolder.ttm
        pred_ttm = trait_table.TraitTableManager(self.pred_traits_f)

        obs_entries = [entry for entry in obs_ttm if entry.name in self.genomes]

        for obs in obs_entries:
            for pred in pred_ttm:
                if obs.name == pred.name:
                    # if changed, make sure to change the best column (right now it assumes > is better)
                    self.results[obs.name] = obs.compare(pred, metric=metric)

        # make sure all genomes were found
        for genome in self.results:
            if self.results[genome] is None:
                raise ValueError("Calculation failed for genome '{}'".format(genome))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs k-folds validation on PICRUSt predictions and reports a table with average accuracy per genome. Optionally includes a bootstrapping step.")

    parser.add_argument("-tree", help="tree that includes all the genomes to test", required=True)
    parser.add_argument("-test", help="an optional file of genome names to test (if not the whole tree)")
    parser.add_argument("-traits", help="a table of traits for each genome to test", required=True)
    parser.add_argument("-k", help="number of groups to use. Default=%(default)s", type=int, default=10)
    parser.add_argument("-bootstrap", help="number of iterations. Default=%(default)s", type=int, default=1)
    parser.add_argument("-metric", help="the metric to use for accuracy", choices=["spearman", "disimilarity"], default="spearman")
    parser.add_argument("-outdir", help="directory to store the output. Default=%(default)s", default=os.getcwd())

    args = parser.parse_args()

    LOG.setLevel(logging.INFO)

    bstrap = BootStrapper(args.tree, args.traits, args.k, args.outdir, to_test=args.test)
    bstrap.run(args.bootstrap, args.metric)
