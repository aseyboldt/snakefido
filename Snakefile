import subprocess
import os
from os.path import join as pjoin
from os.path import exists as pexists
from collections import namedtuple
import json
import glob

BASEDIR = os.getcwd()

OPENMS_BIN = "../OpenMS/bin"
XTANDEM_EXEC = "../tandem-linux-13-09-01-1/bin/64-bit/all_static/tandem.exe"

FASTA = pjoin(BASEDIR, "18Protein_SoCe_Tr_detergents_trace_target_decoy.fasta")
files = glob.glob(pjoin(BASEDIR, 'data', '*.mzML'))
NAMES = [os.path.splitext(os.path.split(p)[1])[0] for p in files]

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

openms = OpenMS(OPENMS_BIN, pjoin(BASEDIR, 'inis'), pjoin(BASEDIR, 'logs'))


def storage(*args):
    return [pjoin(BASEDIR, 'data', x) for x in args]


def workspace(*args):
    return [pjoin(BASEDIR, 'work', x) for x in args]


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
    input: workspace( \
        *['{}.quantified.protein_grps.csv'.format(n) for n in NAMES] \
    )

rule stage_data:
    input: storage("{name}.mzML")
    output: workspace("{name}.mzML")
    run:
        subprocess.check_call(['cp', '--', str(input), str(output)])

rule PeakPickerHiRes:
    input: workspace("{name}.mzML")
    output: workspace("{name}.peaks.mzML")
    params: params('PeakPickerHiRes')
    run:
        openms.PeakPickerHiRes(input, output, ini=params)

rule FeatureFinderCentroided:
    input: workspace("{name}.peaks.mzML")
    output: workspace("{name}.featureXML")
    params: params('FeatureFinderCentroided')
    run:
        openms.FeatureFinderCentroided(input, output, ini=params)

rule XTandemAdapter:
    input: workspace("{name}.mzML")
    output: workspace("{name}.peptides.idXML")
    params: params('XTandemAdapter')
    run:
        extra = [
            '-xtandem_executable', XTANDEM_EXEC,
            '-database', FASTA,
        ]
        openms.XTandemAdapter(input, output, extra_args=extra, ini=params)

rule IDPosteriorErrorProbability:
    input: workspace("{name}.idXML")
    output: workspace("{name}.probs.idXML")
    params: params('IDPosteriorErrorProbability')
    run:
        extra = ['-prob_correct']
        openms.IDPosteriorErrorProbability(
            input, output, extra_args=extra, ini=params
        )

rule PeptideIndexer:
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

rule IDMerger:
    input: workspace(*["{}.peptides.probs.indexed.idXML".format(n) for n in NAMES])
    output: workspace("IDmerged.idXML")
    params: params('IDMerger')
    run:
        openms.IDMerger(input, output, ini=params)

#rule FalseDiscoveryRate:
#    input: workspace("IDmerged.idXML")
#    output: workspace('IDmerged_fdr.idXML')
#    params: params('FalseDiscoveryRate')
#    run:
#        openms.FalseDiscoveryRate(input, output, ini=params)


rule IDRipper:
    input: workspace("IDmerged.idXML")
    output: workspace(*["{}.ripped_consensus.idXML".format(n) for n in NAMES])
    params: params("IDRipper")
    run:
        openms.IDRipper(input, output, ini=params)


rule FidoAdapter:
    input: workspace("IDmerged.idXML")
    output: workspace("protein_grps.idXML")
    params: params('FidoAdapter')
    run:
        openms.FidoAdapter(input, output, ini=params)


rule IDMapper:
    input: *workspace("{name}.featureXML", "{name}.ripped_consensus.idXML")
    output: workspace("{name}.mapped.featureXML")
    params: params('IDMapper')
    run:
        extra = ['-id', input[1]]
        openms.IDMapper(input[0:1], output, extra_args=extra, ini=params)

rule ProteinQuantifier:
    input: *workspace("{name}.mapped.featureXML", "protein_grps.idXML")
    output: workspace("{name}.quantified.protein_grps.csv")
    params: params('ProteinQuantifier')
    run:
        extra = ['-protein_groups', input[1]]
        openms.ProteinQuantifier(
            input[0:1], output, extra_args=extra, ini=params
        )
