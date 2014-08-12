import subprocess
import os
from os.path import join as pjoin
from os.path import exists as pexists
from collections import namedtuple
import json

FASTA = "/home/adr/data/testfido/raw/uniprot-mouse_taxonomy_10090_keyword_181_20140226_CON_REV.fasta"
XTANDEM_EXEC = "../tandem-linux-13-09-01-1/bin/64-bit/all_static/tandem.exe"


class OpenMS:
    def __init__(self, bin_path, ini_dir, log_dir):
        self._bin_path = bin_path
        self._tools = os.listdir(bin_path)
        self._ini_path = ini_dir
        self._log_dir = log_dir

    def __getattr__(self, name):
        if name not in self._tools:
            raise ValueError("OpenMS tool not found: %s" % name)

        def wrapper(input, output, logfile=None, extra_args=None, ini=None):
            if extra_args is None:
                extra_args = []
            if logfile is None:
                logfile = name

            command = [pjoin(self._bin_path, name)]
            command += ['-in'] + input
            command += ['-out'] + output
            if ini is not None:
                if isinstance(ini, list):
                    if not len(ini) == 1:
                        raise ValueError("Invalid params.")
                    ini = ini[0]
                command += ['-ini', json.loads(str(ini))['path']]
            command += extra_args

            log_std = pjoin(self._log_dir, logfile + '.stdout')
            log_err = pjoin(self._log_dir, logfile + '.stderr')
            with open(log_std, 'w') as out:
                with open(log_err, 'w') as err:
                    subprocess.check_call(command, stdout=out, stderr=err)

        return wrapper

openms = OpenMS('../OpenMS/bin', './inis', './logs')


def storage(*args):
    return [pjoin('/home/adr/data/testfido/raw', x) for x in args]


def workspace(*args):
    return [pjoin('/home/adr/data/testfido/work', x) for x in args]


def results(*args):
    return [pjoin('.', x) for x in args]


# store the content of the ini file, so that snakemake will run
# rules agrain if the parameters inside the file change.
def params(name):
    path = pjoin('inis', name + '.ini')
    try:
        with open(path, 'r') as f:
            # TODO replace makes sure that there are no wildcards
            # in the file content. This is not exactly clean ;-)
            data = {'path': path, 'content': f.read().replace('{', '[')}
            return json.dumps(data, sort_keys=True)
    except FileNotFoundError as e:
        raise ValueError("ini file '%s' not found" % path) from e


rule all:
    input: workspace('testfile.quantified.protein_grps.csv')

rule stage_data:
    input: storage("{name}.mzML")
    output: workspace("{name}.mzML")
    run:
        subprocess.check_call(['cp', '--', str(input), str(output)])

rule find_peaks:
    input: workspace("{name}.mzML")
    output: workspace("{name}.peaks.mzML")
    params: params('PeakPickerHiRes')
    run:
        openms.PeakPickerHiRes(input, output, ini=params)

rule find_features:
    input: workspace("{name}.peaks.mzML")
    output: workspace("{name}.featureXML")
    params: params('FeatureFinderCentroided')
    run:
        openms.FeatureFinderCentroided(input, output, ini=params)

rule find_peptides_ms2:
    input: workspace("{name}.mzML")
    output: workspace("{name}.peptides.idXML")
    params: params('XTandemAdapter')
    run:
        extra = [
            '-xtandem_executable', XTANDEM_EXEC,
            '-database', FASTA,
        ]
        openms.XTandemAdapter(input, output, extra_args=extra, ini=params)

rule scores_to_probs:
    input: workspace("{name}.idXML")
    output: workspace("{name}.probs.idXML")
    params: params('IDPosteriorErrorProbability')
    run:
        extra = ['-prob_correct']
        openms.IDPosteriorErrorProbability(
            input, output, extra_args=extra, ini=params
        )

rule index_peptides:
    input: workspace("{name}.idXML")
    output: workspace("{name}.indexed.idXML")
    params: params('PeptideIndexer')
    run:
        extra = [
            '-fasta', FASTA,
            '-annotate_proteins',
            '-enzyme:specificity', 'none',
            '-decoy_string', '_rev'
        ]
        openms.PeptideIndexer(input, output, extra_args=extra, ini=params)

rule identify_proteins:
    input: workspace("{name}.peptides.probs.indexed.idXML")
    output: workspace("{name}.protein_grps.idXML")
    params: params('FidoAdapter')
    run:
        openms.FidoAdapter(input, output, ini=params)

rule map_features_to_peptides:
    input: *workspace("{name}.featureXML", "{name}.peptides.probs.indexed.idXML")
    output: workspace("{name}.mapped.featureXML")
    params: params('IDMapper')
    run:
        extra = ['-id', input[1]]
        openms.IDMapper(input[0:1], output, extra_args=extra, ini=params)

rule quantify_protein_grps:
    input: *workspace("{name}.mapped.featureXML", "{name}.protein_grps.idXML")
    output: workspace("{name}.quantified.protein_grps.csv")
    params: params('ProteinQuantifier')
    run:
        extra = ['-protein_groups', input[1]]
        openms.ProteinQuantifier(
            input[0:1], output, extra_args=extra, ini=params
        )
