import sys, os, re
from posixpath import basename, join, splitext
from datetime import datetime, timedelta
from argparse import ArgumentParser
from collections.abc import Callable, Sequence

from joblib import dump, load

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

import hail as hl
import hailtop.fs as hlfs





### Arguments
clap = ArgumentParser(prog="Hail PCA wrapper",
                      description="Convenience wrapper for performing hail PCA")
clap.add_argument('-r', '--reference', default='GRCh38', choices=['GRCh37', 'GRCh38', 'GRCm38', 'CanFam3'],
                  help='the (hail-supported) reference genome of the samples')
clap.add_argument('-b', '--bucket', default=os.environ.get('WORKSPACE_BUCKET', None),
                  help='cloud bucket prefix to use for saving hail files')
clap.add_argument('--spark-conf', default="spark.executor.memory=2g",
                  help='spark configuration options as <option>=<value>')
clasp = clap.add_subparsers(required=True, metavar='', dest='proc',
                            description='Reference and Sample Operations')

# Arguments for building the reference projcetion
buildref_clap = clasp.add_parser('build-reference',
                                 help="build reference PC projection and ancestry Random-Forest")

file1_args = buildref_clap.add_argument_group('Files')
file1_args.add_argument('reference-vcf',
                        help='reference cohort VCF')
file1_args.add_argument('population-tsv',
                        help='reference sample population assignments, must have: header, all samples in reference-vcf')
file1_args.add_argument('-s', '--samples',
                        help='sample vcf, if supplied will be projected and ancestry-inferred')

buildref_clap.add_argument('-c', '--pop-col', required=True, type=int,
                           help='column with population class in population-tsv')

pca_args = buildref_clap.add_argument_group('PCA')
pca_args.add_argument('-k', default=10, type=int,
                      help='number of PCs to calculate')
pca_args.add_argument('--af-min', default=0.01, type=float,
                      help='minimum allele-frequency filter')
pca_args.add_argument('--hwe-p', default=1e-6, type=float,
                      help='Hardy-Weinberg p-value filter')
pca_args.add_argument('--ld-r2', default=0.1, type=float,
                      help='linkage disequilibrium correlation filter')

# Arguments for projecting and inferring a sample set
infer_clap = clasp.add_parser('infer-samples',
                              help='project samples and infer ancestries using premade reference')

file2_args = infer_clap.add_argument_group('Files')
file2_args.add_argument('sample-vcf',
                        help='sample cohort VCF')
file2_args.add_argument('refloadings',
                        help='Hail table with reference pc loadings and afs, (note: cannot read from local file system)')
file2_args.add_argument('refRFmodel',
                        help='joblib dump of a sklearn RandomForestClassifier trained on reference PCs -> population class')

# Defaults for interactive usage
config = {
    'reference' : 'GRCh38',
    'k': 10, # number of pcs
    'af_min': 0.01,
    'hwe_p': 1e-6,
    'ld_r2': 0.1
}

# Helpful abbreviations
type hlo = hl.MatrixTable | hl.Table
file_readers = {
    '.mt': hl.read_matrix_table,
    '.ht': hl.read_table
}

datetimeformat = '%A %B %#d, %Y, %H:%M:%S' if sys.platform == 'win32' else '%A %B %-d, %Y, %H:%M:%S'
timeformat = '%H:%M:%S'

### Entry
def run(args):
    global config

    config = vars(clap.parse_args(args))
    if config['bucket'] == None:
        print("!! No cloud bucket specified, hail files may be lost (including reference pc variant loadings)")
        config['bucket'] = "."
    config['datadir'] = join(config['bucket'], 'hail', 'data')
    config['start'] = datetime.now()

    sconf = dict(pair.split('=') for pair in config['spark_conf'].split())
    hl.init(app_name='PCA-RF',
            quiet=True,
            tmp_dir=join(config['bucket'], 'hail', 'tmp'),
            spark_conf=sconf,
            log='hail.log')

    if config['proc'] == 'build-reference':
        ref_loadings_ht, rf = build_reference(config['reference-vcf'], config['population-tsv'], config['pop_col'])
        if config['samples'] != None:
            infer_samples(config['samples'], ref_loadings_ht, rf, need_load=False)
    elif config['proc'] == 'infer-samples':
        infer_samples(config['sample-vcf'], config['refloadings'], config['refRFmodel'])

### Pipelines
def build_reference(refvcf, refpoptsv, pop_col):
    stamp('Beginning reference building', True)

    # hail can't import '.bcf' files, but let's keep it portable
    config['filebase'] = re.sub(r"\.[bv]cf\.gz", "", basename(refvcf))
    display_spark_config(mkfname('.config.txt'))
    file_stages = ['unprocessed.mt', 'filtered.mt', ['pcs.ht', 'loadings.ht']]
    progress = resume(file_stages)
    
    mt = None
    if progress < 1:
        stamp('Importing variants')
        mt = import_variants(refvcf, config['reference'],
                             cpn='unprocessed.mt')

    if progress < 2:
        stamp('Filtering')
        mt = prep_mt_pca(mt, config['af_min'], config['hwe_p'], config['ld_r2'],
                         cpn='filtered.mt')
    
    pcs_ht, loadings_ht = None, None
    if progress < 3:
        stamp('PCA')
        pcs_ht, loadings_ht = do_pca(mt, config['k'],
                                     cpn=['pcs.ht', 'loadings.ht'])

    stamp('Preparing random forest data')
    df = prep_df_rf(pcs_ht, refpoptsv, pop_col, config['k'])

    stamp('Planting forest')
    rf = make_rf_model(df, config['k'])
    
    stamp('Reference model complete')
    return loadings_ht, rf


def infer_samples(samplevcf, refloadings, refRF, need_load=True):
    samplebase = re.sub("\.[bv]cf\.gz", "", basename(samplevcf))

    if need_load:
        ref_loadings_ht = hl.read_table(refloadings)
        rf = load(refRF)
    else:
        ref_loadings_ht = refloadings
        rf = refRF

    if rf.n_features_in_ != (nL := len(ref_loadings_ht.loadings.take(1)[0])):
        print(f"Number of Random Forest PC features {rf.n_features_in_} disagrees with number in reference loadings {nL}.", file=sys.stderr)
        sys.exit(1)

    sample_mt = hl.import_vcf(samplevcf, force_bgz=samplevcf.endswith('.gz'), reference_genome=config['reference'])

    # project samples using reference loadings
    sample_pcs_ht = hl.experimental.pc_project(sample_mt.GT, ref_loadings_ht.loadings, ref_loadings_ht.af)

    PCs = PCcols(nL)
    # prepare model data
    sample_data = sample_pcs_ht.to_pandas()
    sample_data.set_index('s', inplace=True)
    sample_data.index.name = 'Sample'
    sample_data = sample_data['scores'].apply(lambda x: pd.Series(x, index=PCs))
    sample_data['Population'] = rf.predict(sample_data[PCs])
    
    ### Output results
    sample_data.to_csv(f'{samplebase}.pca_pop.tsv', sep='\t')


### Atomization
def write_to(rfn, oh: hlo):
    stamp('Writing result to {rfn}')
    oh.write(rfn)

def read_from(rfn) -> hlo:
    stamp('Reading stored result from {rfn}')
    _, ext = splitext(rfn)
    return file_readers[ext](rfn)

def stage(f: Callable[..., Sequence[hlo] | hlo]):
    def checkpoint(*args, cpn: str | Sequence[str], **kwargs):
        cpns = [cpn] if isinstance(cpn, str) else cpn
        base = kwargs.pop("base", None)
        cpfns = [mkfname(n, base) for n in cpns]

        if all(hlfs.exists(fn) for fn in cpfns):
            outs = [read_from(fn) for fn in cpfns]
            out = outs[0] if len(outs) == 1 else outs
        else:
            out = f(*args, **kwargs)
            outs = [out] if not isinstance(out, Sequence) else out
            for fn, o in zip(cpfns, outs): write_to(fn, o)
        
        return out
    return checkpoint

# Stages ('atoms')
@stage
def import_variants(vcf, reference='GRCh38'):
    stamp('Converting VCF to MatrixTable')
    hl.import_vcf(vcf,
                  force_bgz=vcf.endswith('.gz'),
                  reference_genome=reference,
                  array_elements_required=False)\
        .write('cvt.mt', overwrite=True)
    stamp('Reading variant matrix')
    return hl.read_matrix_table('cvt.mt')

@stage
def prep_mt_pca(mt: hl.MatrixTable, aft=0.01, hwe_pt=1e-6, ld_r2=0.1):
    stamp(f'Filtering {mt.count()} variants')

    # missing calls > ref
    filled_gt = hl.if_else(hl.is_defined(mt.GT), mt.GT, hl.Call([0, 0]))
    mt = mt.annotate_entries(GT=filled_gt)
    mt = hl.variant_qc(mt)

    # remove rare variants
    mt = mt.filter_rows(mt.variant_qc.AF[1] > aft)   
    
    # remove variants out of hardy weinberg equilibrium
    mt = mt.filter_rows(mt.variant_qc.p_value_hwe > hwe_pt)

    stamp(f'LD pruning on {mt.count()} variants')
    # remove variants in linkage disequilibrium
    ld_keep = hl.ld_prune(mt.GT, r2=ld_r2)
    mt = mt.filter_rows(hl.is_defined(ld_keep[mt.row_key]))
     
    stamp('Calculating AFs')
    mt = mt.annotate_rows(af = hl.agg.mean(mt.GT.n_alt_alleles()) / 2)

    stamp(mt.count(), 'variants remaining')
    return mt

@stage
def do_pca(mt: hl.MatrixTable, k=5):
    stamp('Calculating PCs')
    _, pcs_ht, loadings_ht = hl.hwe_normalized_pca(mt.GT, k=k, compute_loadings=True)
    stamp('Annotating loadings with AFs')
    loadings_ht = loadings_ht.annotate(af=mt.rows()[loadings_ht.key].af)  
    return pcs_ht, loadings_ht

def prep_df_rf(pcs_ht: hl.Table, refpoptsv, pop_col, k=5):
    PCs = PCcols(k)

    stamp('Converting PCs hail.Table to pandas.DataFrame')
    df = pcs_ht.to_pandas()
    df.set_index('s', inplace=True)
    df.index.name = 'Sample'
    df = df['scores'].apply(lambda x: pd.Series(x, index = PCs)) # Expands nested list

    stamp('Integrating sample populations')
    ref_pops = pd.read_csv(refpoptsv, sep='\t', index_col=[0], usecols=[0, pop_col])
    df['Population'] = ref_pops.loc[df.index]

    stamp('Writing to csv')
    df.to_csv(f'{config['filebase']}.pca_pop.tsv')
    return df

def make_rf_model(df: pd.DataFrame, k=5):
    # train random forest
    PCs = PCcols(k)
    tX = df[PCs]
    ty = df['Population']

    stamp('Training model')
    rf = RandomForestClassifier()
    rf_cvscores = cross_val_score(rf, tX, ty, cv=5)

    stamp('Saving trained model; cv scores: ' + ' '.join(['{:.3f}']*5).format(*rf_cvscores))
    dump(rf, f'{config['filebase']}.pop_rf.sklearn.joblib')

    return rf


### Utilities
def PCcols(n):
    return [f'PC{i}' for i in range(1, n+1)]

# Book keeping
def mkfname(fn, base=None):
    if base is None: base = config['filebase']
    return join(config['datadir'], base + '.' + fn)

def resume(stages):
    for i, stagef in enumerate(stages):
        if isinstance(stagef, str): stagef = [stagef]
        rfns = [mkfname(s) for s in stagef]
        if not all(hlfs.exists(fn) for fn in rfns):
            break
    return i-1

def stamp(message, fulltime=False):
    ct = datetime.now()
    td = format_td(config['start'] - td)

    if not fulltime:
        print(ct.strftime(timeformat),
              f'({td}) {message}')
    else:
        print(ct.strftime(datetimeformat))
        print('  ' + message)

def format_td(td: timedelta):
    s = int(td.total_seconds())
    return f'{s//3600:0>3}:{s//60%60:0>2}:{s%60:0>2}'

# Spark Introspection
def display_spark_config(to=None):
    spark_conf = hl.spark_context().getConf().getAll()
    sc_tree = branchy_tree(spark_conf)
    if to is None:
        display_tree(sc_tree)
        return
    with open(to, 'w') as out:
        display_tree(sc_tree, to)

def branchy_tree(paths, prefix=""):
    if len(paths) == 1: return {paths[0][0]: paths[0][1]}
    
    # find greatest common prefix
    ps = paths[0][0].split('.')
    gcp = ""
    i = 0

    while i < len(ps):
        ncp = gcp + ps[i]

        if all(p.startswith(ncp) for p, q in paths):
            gcp = ncp + "."
            i += 1
        else:
            break

    if i == len(ps):
        raise ValueError(f'Non-unique path encountered: {ncp}')

    tree = dict()
    paths = [(p.removeprefix(gcp), q) for p, q in paths]
    while len(paths) > 0:
        apre = paths[0][0].split('.')[0]
        agrp = [p for p in paths if p[0].startswith(apre)]
        for a in agrp: paths.remove(a)

        tree |= branchy_tree(agrp, gcp)
        
    return {gcp[:-1]: tree}

def display_tree(tree, i=0, to=None):
    for k, v in tree.items():
        print(" |"*i + "-" + str(k), end='', file=to)
        if type(v) == dict:
            print()
            display_tree(v, i+1, to)
        else:
            print(" = " + str(v), to)


if __name__ == "__main__":
    run(sys.argv[1:])