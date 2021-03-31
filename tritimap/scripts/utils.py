import os
import sys
from tritimap import root_dir


def get_snakefile(file="Snakefile"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file;  tried %s" % sf)
    return sf


def get_configfile(file="config.yaml"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the config.yaml file;  tried %s" % sf)
    return sf


def get_samplefile(file="sample.csv"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the sample.csv file;  tried %s" % sf)
    return sf


def get_regionfile(file="region.csv"):
    sf = os.path.join(root_dir, file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the sample.csv file;  tried %s" % sf)
    return sf
