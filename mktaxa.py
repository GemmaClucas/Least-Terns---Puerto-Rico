#! /usr/bin/env python
import sys
from pandas import Series, DataFrame
import numpy as np
from multiprocessing import Pool

if len(sys.argv) != 4:
    print("Usage: ./mktaxa.py ref_reads ref_taxa rep_seqs")
    sys.exit()
    
from qiime2 import Artifact
from qiime2.plugins import feature_classifier
reference_reads = Artifact.load(sys.argv[1])
reference_taxonomy = Artifact.load(sys.argv[2])
repseqs = Artifact.load(sys.argv[3])

def run_blast(pid):
    return feature_classifier.methods.classify_consensus_blast(query = repseqs, reference_reads = reference_reads, reference_taxonomy = reference_taxonomy, maxaccepts = 1, perc_identity = pid, query_cov = .4).classification.view(DataFrame)

confidence = np.linspace(1.0, .7, 80)

if __name__ == '__main__':
    with Pool(8) as pool:
        taxonomy = pool.map(run_blast, confidence)

combined = taxonomy[0]
combined.Confidence = combined.Confidence.astype(float)
for i in range(1, len(confidence)):
    conf = confidence[i]
    result = taxonomy[i]
    result.Confidence = result.Confidence.astype(float) * conf
    for ind in combined.index:
        if combined.loc[ind, "Taxon"] == "Unassigned":
            combined.loc[ind, "Taxon"] = result.Taxon[ind]
            combined.loc[ind, "Confidence"] = result.Confidence[ind]
            
Artifact.import_data("FeatureData[Taxonomy]", combined).save("superblast_taxonomy")
