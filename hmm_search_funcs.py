from Bio.SearchIO import HmmerIO
from Bio import AlignIO, SeqIO
from Bio.AlignIO import StockholmIO
import re 

import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import numpy as np

from glob import glob
from joblib import Parallel, delayed
import os

from tqdm import tqdm

import pandas as pd
from addict import Dict

import pfam
import random

        
def iterate_seq_records(genome_ids, fasta_ids):
    for genome_id, fasta_id in zip(genome_ids, fasta_ids):
        yield get_seq_record(genome_id, fasta_id)

# Load PFAM models and search        
def load_pfam(name):
    if os.path.isfile(f'hmms/{name}.hmm'):
        print(f'file hmms/{name}.hmm already exists')
    else:
        subprocess.run(f'wget https://pfam.xfam.org/family/{name}/hmm -O hmms/{name}.hmm', shell = True)
    if os.path.isfile(f'hmms/{name}_seed.fasta'):
        print(f'file hmms/{name}_seed.fasta already exists')
    else:
        options = 'format=fasta&alnType=seed&order=t&case=l&gaps=dashes&download=1'
        subprocess.run(f'wget "https://pfam.xfam.org/family/{name}/alignment/seed/format?{options}" -O hmms/{name}_seed.fasta', 
                   shell = True)

def run_hmmsearch(genome_path, hmm_name, outdir):
    name = genome_path.split('/')[-1][:-12] # get genome id
    outpath = f'{outdir}/{hmm_name}/{name}'
    os.makedirs(f'{outdir}/{hmm_name}/', exist_ok = True)
    if os.path.isfile(outpath): 
        size = os.path.getsize(outpath)
        if size != 0:
            return
    subprocess.run(f'hmmsearch --cpu 16 --tblout {outpath} hmms/{hmm_name}.hmm {genome_path} >> /dev/null',  
                   shell = True)
    return


# Read search results

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

def parse_hmmresults(path):
    '''parse hmm search results in multiple files in directory 'path' 
    (path contains *)
    returns dataframe'''
    full_data = Dict()
    for path1 in tqdm(glob(path), ncols = 100):
        name = path1.split('/')[-1] #genome name
        r = HmmerIO.hmmer3_tab.Hmmer3TabParser(open(path1, 'r'))
        for i in r:
            for j in i.hits:
                data = parse_description(j.description)
                data['fasta_id'] = j.id
                data['bitscore'] = j.bitscore
                data['evalue'] = j.evalue
                data['genome'] = name
                full_data[(name, j.id)] = data              
    df = pd.DataFrame(full_data).T
    return df


def plot_hmm_results(df,cutoff,name, ax = None):
    if ax is None:
        _, ax = plt.subplots()
        print('new ax')
    if df.empty:
        return
    
    mn = df.bitscore.min()
    mx = df.bitscore.max()
    step = 2
    df[df.bitscore<cutoff] .bitscore.plot.hist(color = 'red',  
                                               bins = np.arange(mn + abs(mn - cutoff)%step - step, cutoff + 0.01, step), 
                                               ax = ax)
    df[df.bitscore>=cutoff].bitscore.plot.hist(color = 'green',
                                               bins = np.arange(cutoff, mx - abs(mx - cutoff)%step+ 0.01, 2), 
                                               ax = ax)
#     ax.set_xlim(0, 150)
    ax.set_title(name + ' ' + pfam_id)
#     ax.set_xlabel('bitscore')
    return ax


# Filter findings
def filter_findings(group, default, keys):
    for pfam in keys:
        cutoff_name = group[pfam].cutoff_name
        cutoff      = group[pfam].cutoff
        print(pfam, cutoff_name, '>', cutoff)
        group[pfam].df_filt    = default[pfam].df[default[pfam].df[cutoff_name] >= cutoff].copy()
        group[pfam].df_filt_no = default[pfam].df[default[pfam].df[cutoff_name] <  cutoff].copy()

def make_final(group, keys):
    # В inct для каждого семейства я записала число находок в каждом геноме
    inct = Dict()
    for key in keys:
        inct[key] = group[key].df_filt.genome.value_counts()
        inct[key].name = key

    # Делаю датасет с числами вхождений всех семейств
    final_incl = pd.DataFrame()
    for key in keys:
        final_incl = final_incl.merge(inct[key], how = 'outer', left_index = True, right_index = True)
    final_incl.fillna(0, inplace = True)
    
    group.inclusion = final_incl
    return

def logical_selection(group, include_keys = [], exclude_keys = []):
    df = group.inclusion
    
    if (len(include_keys) == 0) and (len(exclude_keys) == 0):
        print('no keys'); return
    
    if len(include_keys) == 0: index = df[exclude_keys[0]] == 0
    else:                      index = df[include_keys[0]] > 0
        
    for i in include_keys: 
        index = index & (df[i] > 0)
    for i in exclude_keys: 
        index = index & (df[i] == 0)
        
    group.final_list = df[index]
    
    return 


# Alignments

def get_seq_record(path, fasta_id):
    if type(fasta_id) == str:
        fasta_id = [fasta_id]
    for i in SeqIO.parse(path, 'fasta'):
        if i.name in fasta_id:
            yield i

def write_records(df, fasta_id, path, filename, description = None):
    with open(filename, 'w') as file: 
        for path, chunk in tqdm(df.groupby(path), ncols = 100):
            try:
                records = get_seq_record(path, chunk[fasta_id].to_list())
                for i in records:
                    if description is not None:
                        i.description = chunk[chunk[fasta_id] == i.id][description].iloc[0]

                    SeqIO.write(i, file, 'fasta') 
            except Exception as e:
                print(e)
                print(path)

    subprocess.run(f"sed -i 's/\*//' {filename}", shell = True)# заменить звездочки на ничего
    

def create_sample_hmm(in_file, hmm_file, sample = 50, sample_alignment = None):
    seqs =list(SeqIO.parse(in_file, format = 'fasta'))
    seqs_new = random.sample(seqs, sample)
    if sample_alignment is None:
        sample_alignment = 'temp.fasta'
    with open(sample_alignment, 'w') as file:
        SeqIO.write(seqs_new, file, 'fasta')
    subprocess.run(f'muscle -quiet -in {sample_alignment} -out {sample_alignment}', shell = True)
    subprocess.run(f'hmmbuild {hmm_file} {sample_alignment} >> /dev/null',  shell = True)
    subprocess.run(f'rm temp.fasta', shell = True)
    print(f'created hmm-model {hmm_file}')

    
def align_hmm(hmm_file, in_file, out_file):
    subprocess.run(f"sed -i 's/-//g' {in_file}", shell = True) #Убрать гепы
    subprocess.run(f"sed -i '/^$/d' {in_file}", shell = True) #Убрать пустые линии
    subprocess.run(f'hmmalign {hmm_file} {in_file}  >  temp.stockholm', shell = True) #Выровнять по hmm
    a = list(StockholmIO.StockholmIterator(open(f'temp.stockholm', 'r')))[0]
    SeqIO.write(a, open(out_file, 'w'), 'fasta') #Записать fasta
    os.remove("temp.stockholm") # убрать stockholm
    
    
def sort_alignment(file_name, field = -1):
    a = AlignIO.read(file_name, 'fasta')
    sort_order = dict(zip([i.id for i in a], [float(i.description.split(' ')[field]) for i in a]))
    a.sort(key = lambda i: sort_order[i.id])
    AlignIO.write(a, open(file_name, 'w'), 'fasta')
    


def sample_alignment(in_file, out_file, n = 20):
    random.seed(0)
    a = AlignIO.read(in_file, 'fasta')
    a1 = AlignIO.MultipleSeqAlignment([a[i] for i in random.sample(range(len(a)), n)])
    with open(out_file, 'w') as file:
        AlignIO.write(a1, file, 'fasta')
        
def make_description(df, fields, col = None, assign = True):
    '''Makes description line from concatenating dataframe(df) columns (fields)
    by default(assign) assigns new description to dataframe column(col)'''
    description = df[fields[0]].astype(str)
    if len(fields) == 1:
        return description
    for f in fields[1:]:
        description += ' ' + df[f].astype(str)
    if assign:
        df[col] = description
        return
    return description


def replace_description(in_file, new_desc, out_file):
    a = AlignIO.read(in_file, 'fasta')
    for i in a:
        i.description = new_desc[i.id]
    with open(out_file, 'w') as file:
        AlignIO.write(a, file, 'fasta')
        
        
def rescore_sequences(results_file, in_file = None, hmm_file = None):
    '''Specify either results_file to parse results or all parameters'''
    
    if not in_file is None:
        subprocess.run(f"sed 's/-//g'  {in_file} > temp/hmm_search_file", shell = True) #Убрать гепы
        subprocess.run(f"sed -i '/^$/d' temp/hmm_search_file", shell = True) #Убрать пустые линии
        subprocess.run(f'hmmsearch --cpu 8 --tblout {results_file} {hmm_file} temp/hmm_search_file >> /dev/null', shell = True)
        subprocess.run('rm temp/hmm_search_file', shell = True) 
        
    r = HmmerIO.hmmer3_tab.Hmmer3TabParser(open(results_file, 'r'))
    new_bitscore = {j.id : j.bitscore for i in r for j in i.hits}
    
    return new_bitscore 


def print_families(FAMILIES, keys):
    df = pd.DataFrame(index = keys)
    for k in keys:
        df.at[k, 'name'] = FAMILIES.default[k].name
#         df.at[k, 'cutoff'] = FAMILIES.default[k].cutoff
        df.at[k, 'findings > cutoff'] = FAMILIES.default[k].df_filt.shape[0]
        for group in list(set(FAMILIES.keys()) - set(['default'])):
            df.at[k, f'{group} cutoff']   = FAMILIES[group][k].cutoff
            df.at[k, f'{group} findings'] = FAMILIES[group][k].df_filt.shape[0]
    display(df)
    
def combine_domains(default_group, cutoff_group,  domains):
    '''Assigns genomes df and fasta df'''
    # Собираем в кучу данные о скоре каждого домена
    for i in domains:
        df = default_group[i].df[['fasta_id', 'genome', 'bitscore']].rename({'bitscore': f'{i}_score'}, axis = 1)
        if i == domains[0]: domains_df = df.copy()
        domains_df = domains_df.merge(df, how = 'outer')
    
    # Проверяем проход по cutoff    
    for i in domains:
        domains_df[i] = (domains_df[f'{i}_score'] > cutoff_group[i].cutoff).astype(int)
        
    final_col = 'and'.join(domains)
    domains_df[final_col] = domains_df[domains[0]] & domains_df[domains[1]]# токо на два белка костыль
    
    # записать информацию о последовательностях в df
    cutoff_group[final_col].df = domains_df
        
    return 