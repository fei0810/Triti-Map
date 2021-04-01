import os
import sys
import click
import logging
import subprocess
from shutil import copyfile

from tritimap import __version__
from .scripts import utils


logging.basicConfig(level=logging.INFO,
                    datefmt="%Y-%m-d %H:%M",
                    format="[%(asctime)s %(levelname)s] %(message)s")


def log_exception(msg):
    logging.critical(msg)
    logging.info(
        "Documentation is available at: https://github.com/fei0810/Triti-Map")
    logging.info(
        "Issues can be raised at: https://github.com/fei0810/Triti-Map/issues")
    sys.exit(1)


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """Triti-Map: A Snakemake-based pipeline for gene mapping in Triticeae. For more information, see: https://github.com/fei0810/Triti-Map
    """


def get_snakefile(file="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake pipeline file;  tried %s" %
                 sf)
    return sf


@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    help="Triti-Map main command. The pipeline supports three execute modules: all, only_mapping and only_assembly. First, you need to fill in the configuration file correctly.",
    short_help="Run Triti-Map pipeline"
)
@click.argument(
    "workflow",
    default="all",
    type=click.Choice(["all", "only_mapping", "only_assembly"]),
)
@click.option(
    "-d",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Triti-Map running directory.",
    default="."
)
@click.option(
    "-c",
    "--config-file",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Triti-Map config file, generated with tritimap init.",
)
@click.option(
    "-j",
    "--jobs",
    type=int,
    help="Use at most N CPU cores/jobs in parallel.",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not execute anything, and display what would be done.",
)
@click.option(
    "--profile",
    default=None,
    help="Name of profile to use for configuring Snakemake.",
)
@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)
def run_workflow(workflow, working_dir, config_file, jobs, profile, dryrun,
                 snakemake_args):

    if config_file is None and os.path.isfile('config.yaml'):
        config_file = 'config.yaml'
    elif config_file is None and os.path.isfile(
            os.path.join(working_dir, "config.yaml")):
        config_file = os.path.join(working_dir, "config.yaml")
    else:
        exit(
            "No config.yaml file found. Need generate config file with command: tritimap init"
        )

    cmd = (
        "snakemake --snakefile {snakefile} --directory {working_dir} "
        "{jobs} --rerun-incomplete "
        "--configfile '{config_file}' --nolock "
        " {profile} {dryrun} "
        # " {target_rule} "
        " {args} "
        "--config module={target_rule}").format(
            snakefile=get_snakefile(),
            working_dir=working_dir,
            jobs="--jobs {}".format(jobs) if jobs is not None else "",
            config_file=config_file,
            profile="" if
            (profile is None) else "--profile {}".foramt(profile),
            dryrun="--dryrun" if dryrun else "",
            args=" ".join(snakemake_args),
            target_rule=workflow if workflow != "none" else "")
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logging.critical(e)
        exit(1)


@cli.command(
    'init',
    context_settings={"ignore_unknown_options": True},
    help="Generate snakemake configuration file and other needed file. The command will generate three configuration files(config.yaml, sample.csv and region.csv) in the running directory.",
    short_help="Generate snakemake configuration file and other needed file."
)
@click.option(
    '-d',
    '--working-dir',
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Triti-Map running directory.",
    default="."
)
def init_workdir(working_dir):
    config_file = utils.get_configfile()
    sample_file = utils.get_samplefile()
    region_file = utils.get_regionfile()
    os.makedirs(working_dir)
    copyfile(config_file, os.path.join(working_dir, "config.yaml"))
    copyfile(sample_file, os.path.join(working_dir, "sample.csv"))
    copyfile(region_file, os.path.join(working_dir, "region.csv"))


if __name__ == "__main__":
    cli()
