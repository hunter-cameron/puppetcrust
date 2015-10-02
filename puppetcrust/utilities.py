
def get_ko_by_function(ko_metadata_f, level=2):
    """
    Level 1 is the top level.
    Level 2 is an intermediate level (Corresponding approx to COGs)
    Level 3 is the pathway level
    """
    if level not in [1, 2, 3]:
        raise ValueError("Level must be 1, 2, or 3.")

    data = {}
    with open(ko_metadata_f, 'r') as IN:
        # skip header line
        IN.readline()

        for line in IN:
            fields = line.rstrip().split("\t")

            ko_name = fields[0]
            ko_pathways = fields[-1]

            # multiple pathways sep by "|"
            for pathway in ko_pathways.split("|"):
                levels = pathway.split(";")

                try:
                    data[";".join(levels[:level])].append(ko_name)
                except KeyError:
                    data[";".join(levels[:level])] = [ko_name]
                except IndexError:
                    LOG.warning("{} did not have a pathway at the requested level.".format(ko_name))

    return data

def get_plant_associated_kos(plant_associated_f):
    """ Reads in a database of plant associated kos; returns a dict of lineage: KOs """

    data = {}
    with open(plant_associated_f, 'r') as IN:
        # skip header
        IN.readline()

        for line in IN:
            lineage, kos = line[:-1].split("\t")

            if kos:
                data[lineage] = kos.split(";")
            else:
                data[lineage] = []


    return data

