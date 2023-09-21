import yaml
import os
from argparse import ArgumentParser, Namespace
import re

def get_args() -> Namespace:
    """Parses given arguments
    Returns:
        Namespace: parsed arguments
    """

    parser = ArgumentParser(description="parameters")
    parser.add_argument("--insertionFile", default='Input/BEDInsertionTesting/DLBCLPB.BED', type=str,)
    parser.add_argument("--transposonSystem", default='PB', type=str,)
    parser.add_argument("--mutagenesisMethod", default='predefinedFeatures', type=str,)
    parser.add_argument("--annotation", default='genes', type=str,)
    parser.add_argument("--customAnnotation", default='Input/testAnnot.BED', type=str,)
    parser.add_argument("--multestCorrection", default='bonferroni', type=str,)
    parser.add_argument("--snakemakeRule", default="--cores 8", type=str)
    parser.add_argument("--customFeatures", default="Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed,Input/FeatureBeds/DNASEmesc_SRX1452763.05.bed", type=str) 
    parser.add_argument("--outputDir", default="Output", type=str,)
    parser.add_argument("--downloadDir", default="PrecomputedData", type=str,)
    parser.add_argument("--usePrecomputedFeatures", default="True", type=str,)
    return parser.parse_args()




def get_annotation(annot):
    if annot[0].endswith(".BED") is True:
        annotation = 'custom'
        customAnnotation = annot[0]
    elif bool(re.search("\d*kb|genes", annot[0])) is True: # annot can only be in the format of "genes" or "10kb", "20kb", ...
        annotation = annot[0]
        customAnnotation = False
    else:
        raise ValueError('Annotation can only be in the format of "genes" or "10kb", "20kb", ... or the path to a .BED file')
    return annotation, customAnnotation




args = get_args()
args = vars(args)

args['insertionFile'] = args['insertionFile'].split(',')
args['datasets'] = list(map(lambda x: x.split('/')[-1][:-4], args['insertionFile']))
args['transposonSystem'] = args['transposonSystem'].split(',')
args['multest_correction'] = args['multestCorrection'].split(',')
del(args['multestCorrection'])
args['mutagenesis_method'] = args['mutagenesisMethod'].split(',')
args['mutagenesis_method'] = list(map(lambda mu_model: mu_model.replace('null', 'noMutagenesis'), args['mutagenesis_method']))
del(args['mutagenesisMethod'])
args['annotation'], args['customAnnotation'] = get_annotation(args['annotation'].split(','))
args['customFeatures'] = args['customFeatures'].split(',')

with open('./user_config.yaml', 'w') as file:
    documents = yaml.dump(args, file, default_flow_style=False)

#os.system("./run_slurm.sh" + " " + args['snakemake_rule'])
os.system("snakemake " + args['snakemakeRule'] + " --latency-wait 60")
