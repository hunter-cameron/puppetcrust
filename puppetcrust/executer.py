
import subprocess
import os
import time

class PicrustExecuter(object):
    """ Runs PICRUSt """

    job_id = 0

    @classmethod
    def predict_traits_wf(cls, tree, trait_table, type="trait", limit=None, base_dir=None):
        """ Runs the predict_traits_wf. Returns a name and an output path """
        # make a directory to hold the analysis
        if base_dir is None:
            base_dir = os.getcwd() + "/" + "picrust_project"
        else:
            base_dir = os.path.abspath(base_dir)

        if not os.path.isdir(base_dir):
            os.mkdir(base_dir)

        if type == "trait":
            format_dir = base_dir + "/" + "format_trait"
            predict_out = base_dir + "/" + "predicted_traits.tab"
        elif type == "marker":
            format_dir = base_dir + "/" + "format_marker"
            predict_out = base_dir + "/" + "predicted_markers.tab"
        else:
            raise ValueError("type must be one of 'trait', 'marker'")

        # convert limit list to mock OTU table 
        if limit is not None:
            with open(limit, 'r') as IN, open(base_dir + "/" + "mock_otu.txt", 'w') as OUT:
                for line in IN:
                    line = line.strip()
                    line += "\t0\n"
                    OUT.write(line)

            limit = base_dir + "/" + "mock_otu.txt"


        format_cmd = cls._get_format_command(trait_table, tree, format_dir)

        # formatted paths
        fmt_table = format_dir + "/" + "trait_table.tab"
        fmt_tree = format_dir + "/" + "reference_tree.newick"
        prun_tree = format_dir + "/" + "pruned_tree.newick"


        asr_out = format_dir + "/" + "asr.tab"
        reconstruct_cmd = cls._get_asr_command(fmt_table, prun_tree, asr_out)

        predict_cmd = cls._get_predict_traits_command(fmt_table, asr_out, fmt_tree, predict_out, limit=limit)


        # link all the necessary commands into a single command
        super_command = "; ".join([format_cmd, reconstruct_cmd, predict_cmd])

        job_name = "picrust_cmd{}".format(cls.job_id)
        subprocess.call([   "bsub",
                            #"-q", "bigmem",
                            "-o", "{}/auto_picrust.out".format(base_dir),
                            "-e", "{}/auto_picrust.err".format(base_dir),
                            "-J", job_name,
                            super_command
                            ])
        cls.job_id += 1

        return job_name, predict_out

    @classmethod
    def predict_metagenome(cls, otu_table, copy_numbers, trait_table, base_dir=None):
        # make a directory to hold the analysis
        if base_dir is None:
            base_dir = os.getcwd() + "/" + "picrust_project"
        else:
            base_dir = os.path.abspath(base_dir)

        if not os.path.isdir(base_dir):
            os.mkdir(base_dir)

        filter_out = base_dir + "/" + "filtered_OTU_table.tab"
        cls._filter_otus(otu_table, trait_table, filter_out)

        convert_out = base_dir + "/" + "filtered_OTU_table.biom"
        convert_cmd = cls._get_biom_convert_command(filter_out, convert_out)

        norm_out = base_dir + "/" + "normalized_OTU_table.biom"
        norm_cmd = cls._get_normalize_command(convert_out, copy_numbers, norm_out)

        predict_out = base_dir + "/" + "predicted_metagenome.tab"
        predict_cmd = cls._get_predict_metagenome_command(norm_out, trait_table, out=predict_out)

        # link all the necessary commands into a single command
        super_command = "; ".join([convert_cmd, norm_cmd, predict_cmd])

        job_name = "picrust_cmd{}".format(cls.job_id)
        subprocess.call([   "bsub",
                            #"-q", "bigmem",
                            "-o", "{}/auto_picrust.out".format(base_dir),
                            "-e", "{}/auto_picrust.err".format(base_dir),
                            "-J", job_name,
                            super_command
                            ])
        cls.job_id += 1

        return job_name, predict_out

    @classmethod
    def wait_for_job(cls, job_name="picrust_cmd*"):
        """ waits for job to complete, checks every 10 seconds """
        while cls._job_running(job_name):
            time.sleep(10)

    @staticmethod
    def _job_running(job_name="picrust_cmd*"):

        output = subprocess.check_output([
                    "bjobs",
                    "-J", "picrust_cmd*"
                ])
        #print(output)
        if output:
            return True
        else:
            return False

    @staticmethod
    def _get_format_command(trait_tab, tree, out):
        exe = subprocess.check_output(["which", "format_tree_and_trait_table.py"]).strip()
        format_files = "python {exe} -t {tree} -i {trait_tab} -o {out}".format(exe=exe, tree=tree, trait_tab=trait_tab, out=out)

        return format_files

    @staticmethod
    def _get_asr_command(trait_table, tree, out):
        exe = subprocess.check_output(["which", "ancestral_state_reconstruction.py"]).strip()
        asr = "python {exe} -i {trait_table} -t {tree} -o {out}".format(exe=exe, trait_table=trait_table, tree=tree, out=out)

        return asr

    @staticmethod
    def _get_predict_traits_command(trait_table, asr_table, tree, out, limit=None):
        exe = subprocess.check_output(["which", "predict_traits.py"]).strip()
        predict = "python {exe} -i {trait_table} -t {tree} -r {asr_table} -o {out} -a".format(exe=exe, trait_table=trait_table, asr_table=asr_table, tree=tree, out=out)

        # add optional limit predictions to list or OTU table
        if limit is not None:
            predict += " -l {}".format(limit)

        return predict

    @staticmethod
    def _filter_otus(otu_f, traits_f, out_f="filtered_otu_table.tab"):
        """ Filters an otu table to only include OTUs that are in the trait table. Writes new table to out_f """

        # read in a list of OTUs in the trait table
        otus_to_keep = []
        with open(traits_f, "r") as IN:
            for line in IN:

                # skip comments
                if line.startswith("#"): 
                    continue

                name, counts = line.split("\t", 1)
                otus_to_keep.append(name)

        # filter table and write new one
        found_header = False        # this stores whether a line that could be the header has been found
        num_filtered = 0
        with open(otu_f, 'r') as IN, open(out_f, 'w') as OUT:
            for line in IN:

                # write comments
                if line.startswith("#"):
                    OUT.write(line)
                    continue
                
                name, counts = line.split("\t", 1)

                if name in otus_to_keep:
                    OUT.write(line)
                else:
                    # now we need to consider whether or not this could be the header
                    if found_header is False:

                        # we will do this by trying to convert a few of the "count" fields into
                        # numbers in the hope that sample names will be non-numeric
                        first, second, third, rest = counts.split("\t", 3)
                        try:
                            float(first)
                            float(second)
                            float(third)

                        except ValueError:      # one of the 3 had an alphabet character
                            OUT.write(line)
                            continue 

                        # if they were all numbers, lets assume we found the header
                        found_header = True
                    
                    num_filtered += 1

        print("Removed {} OTUs from the OTU table that had no predicted traits.".format(num_filtered))

    @staticmethod
    def _get_biom_convert_command(otu_f, out_f):
        # cannot call it directly with python or there is an error
        #exe = subprocess.check_output(["which", "biom"]).strip()
        convert = "biom convert -i {otu_f} --table-type 'otu table' -o {out_f}".format(otu_f=otu_f, out_f=out_f)

        return convert

    @staticmethod
    def _get_normalize_command(otu_table, copy_numbers, out):
        exe = subprocess.check_output(["which", "normalize_by_copy_number.py"]).strip()
        normalize = "python {exe} -i {otu_table} -c {copy_numbers} -o {out}".format(exe=exe, otu_table=otu_table, copy_numbers=copy_numbers, out=out)

        return normalize

    @staticmethod
    def _get_predict_metagenome_command(otu_table, trait_table, out):
        exe = subprocess.check_output(["which", "predict_metagenomes.py"]).strip()
        predict = "python {exe} -i {otu_table} -c {trait_table} -o {out} -f".format(exe=exe, otu_table=otu_table, trait_table=trait_table, out=out)

        return predict
