import os
import time
import matplotlib.pyplot as plt
import numpy as np
import copy
import scipy.cluster.hierarchy as sch
from scipy.stats import scoreatpercentile
import matplotlib.image as mpimg
from IPython.display import display
from IPython.display import Image
from Bio.Seq import Seq
from Bio import motifs
import colorsys
from pysca import scaTools as sca
import pickle as pickle
from optparse import OptionParser
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process an input file, an input PDB, and produce a PML output."
    )
    parser.add_argument(
        '-i', '--infile',
        required=True,
        help='Path to the input file'
    )
    parser.add_argument(
        '-o', '--output_pml',
        required=True,
        help='Path for the output PML file'
    )
    parser.add_argument(
        '-p', '--pdb_file',
        required=True,
        help='Path to the input PDB file'
    )
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    infile = args.infile
    output_pml = args.output_pml
    pdb_file = args.pdb_file
    
    print(f"Reading input file: {infile}")
    print(f"Reading PDB file: {pdb_file}")
    print(f"Writing PML output to: {output_pml}")
    
    db = pickle.load(open(infile,"rb"))
    Dseq = db['sequence']
    Dsca = db['sca']
    Dsect = db['sector']
    
    print("After processing, the alignment size is %i sequences and %i positions" % (Dseq['Nseq'], Dseq['Npos']))
    print("With sequence weights, there are %i effective sequences" % (Dseq['effseqs']))
    
    Dseq["ats"]
    
    listS = [Dsca['simMat'][i,j] for i in range(Dsca['simMat'].shape[0]) for j in range(i+1, Dsca['simMat'].shape[1])]
    Z = sch.linkage(Dsca['simMat'],method = 'complete', metric = 'cityblock')
    R = sch.dendrogram(Z,no_plot = True)
    ind = R['leaves']
    
    for n,ipos in enumerate(Dsect['ics']):
        sort_ipos = sorted(ipos.items)
        ats_ipos = ([Dseq['ats'][s] for s in sort_ipos])
        ic_pymol = ('+'.join(map(str, ats_ipos)))
        print('IC %i is composed of %i positions:' % (n+1,len(ats_ipos)))
        print(ic_pymol + "\n")
    
    sectors = list()
    colors = [0.0, 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88]
    
    for n, ic in enumerate(Dsect['ics']):
        s = sca.Unit()
        s.items = sorted(ic.items)
        s.col = colors[n] if n < len(colors) else 0.66
        sectors.append(s)
    
    for i,k in enumerate(sectors):
        sort_ipos = sorted(k.items)
        ats_ipos = ([Dseq['ats'][s] for s in sort_ipos])
        print(len(ats_ipos))
        ic_pymol = ('+'.join(map(str, ats_ipos)))
        print('Sector %i is composed of %i positions:' % (i+1,len(ats_ipos)))
        print(ic_pymol + "\n")
    
    print(Dseq["ats"])
    
    basename = os.path.basename(pdb_file)
    pdb = os.path.splitext(basename)[0]
    
    with open("ats_" + os.path.basename(infile) +".txt", "w") as f:
        f.write(','.join(map(str, Dseq['ats'])))
    
    ats_str = [str(x) for x in Dseq['ats']]
    sca.writePymol(pdb, sectors, Dsect['ics'], ats_str, output_pml,'A', './', 0)