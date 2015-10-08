
from Bio import SeqIO
import logging

from puppetcrust.trait_table import TraitTableManager


logging.basicConfig()
LOG = logging.getLogger(__name__)



class DatabaseManager(object):
    """ Manages files to make custom databases """

    def __init__(self):
        self.marker_fastas = []
        self.trait_tables = []

    def add_fasta(self, fasta_f):
        self.marker_fastas.append(fasta_f)

    def add_trait_table(self, table_f):
        self.trait_tables.append(table_f)

    def generate_database(self, prefix="new_database", subset=[], inverse=False, verbose=False):
        """ 
        Concatenates the files and removes duplicates and ensures trait tables and marker fastas have matching entries. Returns a tuple (output_fasta, output_traits). 
        
        Optionally, keeps only a subset of headers. 

        Set inverse to remove the subset instead.

        Subset is done when parsing the FASTA because the Trait Table is scaled to only things kept in 
        the FASTA.
        """
       
        # set up the output paths 
        output_fasta = prefix + ".markers.fasta"
        output_traits = prefix + ".traits.tab"

        # convert subset to a dict to store if it has been found
        subset = {k: 0 for k in subset}

        #
        ## Make a list of FASTA seqs to write
        #

        # store the found sequences indexed by header
        genomes = {}
        for fasta_f in self.marker_fastas:
            with open(fasta_f, 'r') as IN:
                for record in SeqIO.parse(IN, "fasta"):
                    if record.id in genomes:
                        LOG.warning("Duplicate header found in FASTA files: '{}'".format(record.id))
                    else:

                        # check if we are subsetting anything
                        # write it unless found && inverse OR !found && !inverse
                        if subset:
                            if record.id in subset:
                                if inverse:
                                    subset[record.id] = 1
                                    continue
                            else:
                                if not inverse:
                                    subset[record.id] = 1
                                    continue
                        genomes[record.id] = record

        # issue warnings if necessary
        if subset:
            for name in subset:
                if subset[name] != 1:
                    LOG.warning("Some ids in subset list were not found in marker FASTA. Continuing with analysis anyway...")
                    break


        #
        ## Write Trait Table
        #
        traits = None
        genomes_found = {g: 0 for g in genomes}
        with open(output_traits, 'w') as OUT:
            for trait_f in self.trait_tables:
                
                ttm = TraitTableManager(trait_f)
                
                # make sure the traits we are using are consistent throughout
                if traits is None:
                    traits = ttm.traits
                    OUT.write("\t".join(["genome"] + traits) + "\n")
                else:
                    # must use sets to ignore order
                    if set(traits) != set(ttm.traits): 
                        LOG.warning("Traits from '{}' don't match traits from the first table. Traits from the first table will be used to write the final one.".format(trait_f))

                # process each entry in the table
                for entry in ttm:
                    try:
                        if genomes_found[entry.name] == 0:
                            genomes_found[entry.name] += 1
                            entry.write(OUT, traits)
                        else:
                            # header already has an entry in the trait table. Skip.
                            continue
                    except:
                        # header is not present in the FASTA. Skip.
                        continue

        # issue more warnings if necessary
        not_found = []
        for genome in genomes_found:
            if genomes_found[genome] == 0:
                not_found.append(genome)
                if not verbose:
                    break
        
        if not_found:
            LOG.warning("Some genomes in the marker file(s) were not found in the trait table(s).")

        if verbose:
            print("{} were not found (these will not be in the output FASTA):".format(str(len(not_found))))
            for g in not_found:
                print(g)

        # write the fasta
        with open(output_fasta, 'w') as OUT:
            for genome in genomes:
                if genome in not_found:
                    continue
                else:
                    SeqIO.write(genomes[genome], OUT, "fasta")


        return output_fasta, output_traits

