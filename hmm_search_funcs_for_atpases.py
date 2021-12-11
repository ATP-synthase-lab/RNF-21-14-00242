from Bio.SearchIO import HmmerIO
from Bio import AlignIO, SeqIO
import re 

import pandas as pd
import subprocess

from glob import glob
import os

from tqdm.notebook import tqdm

import pandas as pd
from addict import Dict
import sys
sys.path.append('../03_acetyl')
from funcs import *

from matplotlib_venn import venn2

def run_hmmsearch(genome_path, hmm_name, outdir, mode = 'tblout'):
    '''mode tblout (default) or domtblout'''
    name = genome_path.split('/')[-1][:-12] # get genome id
    
    outpath   = f'{outdir}/{hmm_name}/{name}'
    os.makedirs(f'{outdir}/{hmm_name}/', exist_ok = True)
    
    if mode == 'domtblout':
        outpath   = f'{outdir}/{hmm_name}_dom/{name}'
        os.makedirs(f'{outdir}/{hmm_name}_dom/', exist_ok = True)
    
    if os.path.isfile(outpath): 
        size = os.path.getsize(outpath)
        if size != 0:  return
        
    subprocess.run(f'hmmsearch --cpu 16 --{mode} {outpath} hmms/{hmm_name}.hmm {genome_path} >> /dev/null',  
                   shell = True)
    return



def run_hmmsearch_1(hmm_name, outdir, mode = 'tblout'):
    '''mode -- default tblout or domtblout '''
    for n, path in enumerate(glob('/media/pc208/es/litvinanna/hmmdb/*')):
        outpath = f'{outdir}/{hmm_name}'
        chunk = path.split('/')[-1].split('.')[0]
        r = subprocess.run(f'hmmsearch --cpu 16 --{mode} {outpath}_{chunk} hmms/{hmm_name}.hmm {path} >> /dev/null',  
                   shell = True, capture_output = True)
        print(r.stdout, end = '')
        print(r.stderr, end = '')
        
    if mode == 'tblout':
        final = pd.concat([parse_hmmresults_1(path) for path in glob(f'{outpath}_*')])  
        final.to_csv(outpath, index = None)
        
    return

def parse_description(line):
    '''parse description line from hmm search results
    return dict
    line example:
    # 3614 # 5830 # -1 # ID=178_5;partial=00;start_type=ATG;rbs_motif=AGGAGG;rbs_spacer=5-10bp;gc_cont=0.524''' 
    items = re.split("[# ;]+",line)
    items = list(filter(None, items))
    desc = {}
    for i, k in enumerate(['start', 'end', 'strand']):
        desc[k] = int(items[i])
    for i in items:
        if '=' in i:
            (k, v) = i.split('=')
            desc[k] = v
    return desc

def parse_hmmresults_dom(path):
    r = HmmerIO.hmmer3_domtab.Hmmer3DomtabHmmhitParser(open(path, 'r') )
    name = path.split('/')[-1] #genome name
    data = {}
    for i in r:
        for j in i.hits:
            seq_description = parse_description(j.description)
            for k in j.hsps:
                hit_info = parse_hit(k)
                hit_info['fasta_id'] = j.id
                hit_info['dom_n'] = k.domain_index
                
                hit_info['genome'] = name
                hit_info.update(seq_description)
                data[(j.id, k.domain_index)] = hit_info 

    data = pd.DataFrame(data).T
    return data


def parse_hit(hit):
    data = {}
    data['bitscore'] = hit.bitscore
    data['evalue'] = hit.evalue
    data['hit_start'] = hit.hit_start
    data['hit_end'] = hit.hit_end
    data['hit_span'] = hit.hit_span
    data['query_start'] = hit.query_start
    data['query_end'] = hit.query_end
    data['query_span'] = hit.query_span
    data['query_span'] = hit.query_span

    return data

def parse_hmmresults_1(path):
    full_data = Dict()
    r = HmmerIO.hmmer3_tab.Hmmer3TabParser(open(path, 'r'))
    for i in r:
        for j in i.hits:
            data = parse_description(j.description)
            data['fasta_id'] = j.id
            data['bitscore'] = j.bitscore
            data['evalue'] = j.evalue
            full_data[j.id] = data       
    df = pd.DataFrame(full_data).T
    return df


''''acc_avg', 'aln', 'aln_all', 'aln_annotation', 'aln_annotation_all', 'aln_span', 'bias', 'bitscore', 'domain_index', 'env_end', 'env_start', 'evalue', 'evalue_cond', 'fragment', 'fragments', 'hit', 'hit_all', 'hit_description', 'hit_end', 'hit_end_all', 'hit_features', 'hit_features_all', 'hit_frame', 'hit_frame_all', 'hit_id', 'hit_inter_ranges', 'hit_inter_spans', 'hit_range', 'hit_range_all', 'hit_span', 'hit_span_all', 'hit_start', 'hit_start_all', 'hit_strand', 'hit_strand_all', 'is_fragmented', 'molecule_type', 'output_index', 'query', 'query_all', 'query_description', 'query_end', 'query_end_all', 'query_features', 'query_features_all', 'query_frame', 'query_frame_all', 'query_id', 'query_inter_ranges', 'query_inter_spans', 'query_range', 'query_range_all', 'query_span', 'query_span_all', 'query_start', 'query_start_all', 'query_strand', 'query_strand_all']'''
    