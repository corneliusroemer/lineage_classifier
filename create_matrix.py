#%%
from typing import OrderedDict
from sqlalchemy import create_engine,Column, String
from sqlalchemy.orm import sessionmaker, declarative_base, Session
from sqlalchemy.sql.sqltypes import BLOB
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyzstd import ZstdDict, decompress
import pandas as pd
import numpy as np
import re
import tempfile
import ray
import traceback
import psutil
import gc
from collections import OrderedDict

Base = declarative_base()

class Fasta(Base):
    __tablename__ = 'fasta'

    strain = Column(String, primary_key=True)
    sequence = Column(BLOB)

    def __repr__(self):
        return f"Sequence(strain={self.strain!r}, name={self.sequence!r})"

class Zstd_dict_table(Base):
    __tablename__ = 'zstd_dict'

    id = Column(String, primary_key=True)
    dictionary = Column(BLOB)

    def __repr__(self):
        return f"Dictionary(id={self.id!r}, name={ZstdDict(self.dictionary)!r})"

def connect_to_db(path, debug=False):
    return create_engine(f"sqlite+pysqlite:///{path}",echo=debug)

def start_session(engine)->Session:
    Session = sessionmaker(bind=engine)
    return Session()

def close_session(session):
    session.close()

def nextstrain_from_pango_strain_name(nextstrain_name):
    """Convert a nextstrain name to a pango name"""
    # Remove all underscores until the first slash using a Regex
    split_name = nextstrain_name.split('/')
    split_name[0] = split_name[0].replace('_','')
    return '/'.join(split_name)

def get_seq_from_db(db_path, strain_names: list) -> list:
    """Read in a db and strain names and returns list of Seq"""
    engine = connect_to_db(db_path)
    session = start_session(engine)
    provided_dict = session.query(Zstd_dict_table).one_or_none()
    zd = ZstdDict(provided_dict.dictionary)
    nextstrain_strain_names = [nextstrain_from_pango_strain_name(i) for i in strain_names]
    sequences = session.query(Fasta).filter(Fasta.strain.in_(nextstrain_strain_names)).all()
    output_sequences = []
    for sequence in sequences:
        record = SeqRecord(
            Seq(decompress(sequence.sequence,zstd_dict=zd).decode('UTF-8')),
            id=sequence.strain,
            description=''
        )
        output_sequences.append(record)
    close_session(session)
    return output_sequences

alpha_digits_regex = re.compile('[^A-Z0-9/]')

def canonicalize_string(in_string):
    """Convert a string to a canonical form"""
    return alpha_digits_regex.sub('',in_string.upper())

def create_canonical_to_pango_lookup() -> dict:
    """Create lookup table from canonical to pango designation"""
    pango_to_canonical = {strain:canonicalize_string(strain) for strain in pd.read_csv('lineages.csv').taxon}
    canonical_to_pango = {canonical:strain for strain, canonical in pango_to_canonical.items()}
    return canonical_to_pango 

canonical_to_pango = create_canonical_to_pango_lookup()

def pango_from_nextstrain(nextstrain_name):
    """Convert a nextstrain name to a pango name"""
    return canonical_to_pango[canonicalize_string(nextstrain_name)]

def get_seq_not_in_db(db_path='gisaid.fasta.db') -> list:
    """Return strains missing from db"""
    pango_strains = {strain:canonicalize_string(strain) for strain in pd.read_csv('lineages.csv').taxon}
    engine = connect_to_db(db_path)
    session = start_session(engine)
    all_strains_in_db = {canonicalize_string(i[0]):i[0] for i in session.query(Fasta.strain).all()}
    close_session(session)
    not_in_db = []
    counter = 0
    strain_translation = {'not_in_db': []}
    for strain in pango_strains.keys():
        canonical = pango_strains[strain]
        counter += 1
        if counter % 1000 == 0:
            print(f"{counter} processed, {len(strain_translation['not_in_db'])} missing so far")
        if canonical in all_strains_in_db:
            strain_translation[canonical] = all_strains_in_db[canonical]
        else:
            strain_translation['not_in_db'].append(strain)
            print(strain)
    return strain_translation

def designated_strains_from_pango(lineage_name,lineages_path='lineages.csv') -> list:
    """Read in designation csv and return a list of strains matching pango designation"""
    designations = pd.read_csv(lineages_path)
    return designations[designations.lineage == lineage_name].taxon.tolist()

def create_nuc_offset_array() -> np.array:
    nuc_dict = {'A':0, 'C':1, 'G':2, 'T':3, '-':4}
    alphabet_length = len(nuc_dict)
    offset_base = np.zeros(6,dtype=np.uint32)
    offset_base[alphabet_length] = 1
    offset_dict = np.full((255,6),offset_base)
    for key, value in nuc_dict.items():
        offset_dict[ord(key)][value] = 1
        offset_dict[ord(key)][alphabet_length] = 0
    return offset_dict

offset_array = create_nuc_offset_array()
alphabet_length = len(offset_array[0])

def to_array(c:np.uint8) -> np.array:
    return offset_array[c]

def seq_to_array(seq: SeqRecord) -> np.array:
    # Maybe speed up by using sparse notation then converting?
    sequence = np.frombuffer(seq.seq.encode(encoding='ascii'),dtype=np.uint8).copy()

    array = np.apply_along_axis(to_array,0,sequence)

    count = 0
    while (array[count][4]):
        try:
            array[count][4]=0
            array[count][5]=1
        except:
            print(f"{seq.id} @ {count} : {array[0:count+2]} : {str(sequence[0:count+2])}")
            traceback.print_stack()
        count += 1
    count = len(sequence) - 1
    while (array[count][4]):
        try:
            array[count][4]=0
            array[count][5]=1
        except:
            # print(f"{seq.id} @ {count} : {sequence[count-5:]}")
            pass
        count -= 1
    return array

@ray.remote
def aggregate_seqs(seqs: list) -> np.ndarray:
    """Aggregate a list of SeqRecords into a single numpy array"""
    aggregate = np.zeros((29903,6), dtype=np.uint32)
    for count, seq in enumerate(seqs):
        # if count % 10000 == 0:
            # print(f"{seq.id:'%6s'}:{count:'%6u'} of {len(seqs):'%6u'} processed")
        # if count > 1000:
        #     break
        aggregate += seq_to_array(seq)
    # print(f"{seq.id:'%6s'}:{count:'%6u'} of {len(seqs):'%6u'} processed")
    return aggregate

def get_lineages_from_designation(lineages_path='lineages.csv'):
    """Read in a csv and return a list of lineages"""
    df = pd.read_csv(lineages_path)
    return df['lineage']

def auto_garbage_collect(pct=80.0):
    """
    auto_garbage_collection - Call the garbage collection if memory used is greater than 80% of total available memory.
                            This is called to deal with an issue in Ray not freeing up used memory.

        pct - Default value of 80%.  Amount of memory in use that triggers the garbage collection call.
    """
    if psutil.virtual_memory().percent >= pct:
        gc.collect()
    return

def generate_series_from_pango(lineage_name,db_path='gisaid.fasta.db'):
    """Generate a pandas series from a pango designation"""
    strains = designated_strains_from_pango(lineage_name)
    seqs = get_seq_from_db(db_path, strains)
    print(f"{lineage_name:<10}:{len(strains)-len(seqs):6d} missing of {len(strains):6d} designated strains")
    return aggregate_seqs.remote(seqs)

def generate_df_from_pangos(lineage_names,db_path='gisaid.fasta.db') -> dict:
    """Generate a pandas dataframe from a list of pango designations"""
    matrix = {'index': [('index','lineage','counts')]}
    matrix['array'] = np.zeros((len(lineage_names),29903,6), dtype=np.uint32)
    parallel_count = 0
    futures = []
    for i, lineage_name in enumerate(lineage_names):
        futures.append(generate_series_from_pango(lineage_name))
    for i, result in enumerate(ray.get(futures)):
        matrix['array'][i] = result
    for i, lineage_name in enumerate(lineage_names):
        matrix['index'].append((i,lineage_name,np.sum(matrix['array'][i][0])))
    return matrix

def generate_df_from_all(db_path='gisaid.fasta.db'):
    """Generate a pandas dataframe from all pango designations"""
    return generate_df_from_pangos(get_lineages_from_designation().unique().tolist())

def export_np_matrix(df,filename=tempfile.TemporaryFile()):
    """Export a pandas dataframe to a numpy matrix"""
    np.savez_compressed(filename,array=df['array'],index=df['index'])

def import_np_matrix(filename) -> dict:
    """Import a numpy matrix from a file"""
    return np.load(filename)

def sequence_to_likelihood(seq: SeqRecord, m=None) -> OrderedDict:
    """Generate a likelihood for each pango lineage"""
    if m is not None:
        matrix = m['matrix']
        index = m['index']
    else:
        matrix_dict = import_np_matrix('full_run.npz')
        index = matrix_dict['index']
        matrix = matrix_dict['array'].astype(np.float32)
        #Flatten matrix and array for multiplication
        #Normalise for anything but Ns
        #Set all Ns to 0, multiply by 1/(count-N)
        #Add constant to matrix to account for new mutations
        matrix[:,:,5] = 0
        matrix += 0.001
        matrix = matrix / matrix.sum(axis=2)[:,:,np.newaxis]
        # Masking bad sites S:95 and S:142
        #matrix[:,[21845,21986],:] = 1
        matrix = np.log(matrix)
        np.save('matrix.npy',matrix)
        np.save('index.npy',index)
    seq_vec = seq_to_array(seq).astype(bool)
    seq_vec[:100] = 0
    seq_vec[-100:] = 0
    seq_vec[:,5] = 0
    # ref = seq_to_array(get_seq_from_db('gisaid.fasta.db', ['Wuhan/Hu-1/2019'])[0]).astype(bool)
    # Set to 1 if >0.5 remain at reference
    # mask = np.logical_and((matrix > 0.1).astype(bool), np.repeat([ref],len(matrix),axis=0))
    # np.putmask(matrix,mask,1)
    # Mask with reference
    # Then check if >0.5 to refine mask
    # Then set to 1 for mask
    #Normalise by dividing by total count
    #Or simpler: just use seq_array as mask then sum up 
    res = pd.Series(np.add.reduce(matrix,axis=(1,2),where=seq_vec),index=index[1:,1]).sort_values(ascending=False)[:].to_dict(OrderedDict)

    #Maybe do postchecking for parent? Does it have the necessary defining mutations? If not -> it's not this lineage. -> Local sanity checking
    #Find mutations not present in parent: If subset of these present -> it's child

    #Penalize missing mutations more -> not the real solution here

    return res

@ray.remote
def sequences_to_pango(seqs: list) -> list:
    df = pd.DataFrame(columns=['strain','designated_lineage','likelihood'])
    series = [] 
    m = {'matrix': np.load('matrix.npy'), 'index' : np.load('index.npy')}
    for seq in seqs:
        series.append((pango_from_nextstrain(seq.id), sequence_to_likelihood(seq,m)))
        # print(f"{seq.id:<40}  {series[-1].index[0]:<10}\n")  
    return series

def multi_sequences_to_pango(seqs: list,threads=1) -> pd.DataFrame:
    n = len(seqs)
    futures = []
    for i in range(threads):
        futures.append(sequences_to_pango.remote(seqs[i*n//threads:(i+1)*n//threads]))
    results = []
    for result in ray.get(futures):
        results.extend(result)
    ray.shutdown()
    df = pd.DataFrame.from_records(results,index='strain',columns=['strain','assignment'])
    return df

def test_assign() -> pd.DataFrame:
    lin = pd.read_csv('lineages.csv')
    samples = lin.groupby('lineage').agg(pd.DataFrame.sample)[:]
    df = multi_sequences_to_pango(get_seq_from_db('gisaid.fasta.db',samples.taxon.tolist()),threads=6)
    return df

def test_run() -> pd.DataFrame:
    lin = pd.read_csv('lineages.csv')
    l = lin.set_index('taxon')
    joined = test_assign().join(l)
    joined['assigned_lineage'] = joined.assignment.apply(lambda x: list(x)[0])
    joined['second'] = joined.assignment.apply(lambda x: list(x)[1])
    joined['likelihood'] = joined.assignment.apply(lambda x: x[list(x)[0]])
    joined['diff_second'] = joined.assignment.apply(lambda x: x[list(x)[1]]-x[list(x)[0]])
    try:
        joined['diff_true'] = joined.apply(lambda x: x[0][x[1]]-x[0][x[2]] if x[1] != x[2] else np.NAN,axis=1)
    except:
        pass
    print(joined[joined.diff_true.notna()])
    return joined
    # for sample in samples:
    #     df.loc[('Australia/QLD1080/2020',slice(None)),'designation'] = sample.lineage
    # for row in df.loc[(:
    # print(df.iloc[:,0])
#%%
if __name__ == '__main__':
    print(test_run())

# lin = pd.read_csv('lineages.csv')
# samples = lin.groupby('lineage').agg(pd.DataFrame.sample)[0:50]

# df = pd.DataFrame(multi_sequences_to_pango(get_seq_from_db('gisaid.fasta.db',samples.taxon.tolist()),threads=6),columns=['likelihood'])


# #%%
# lin = pd.read_csv('lineages.csv',index='taxon')
# samples = lin.groupby('lineage').agg(pd.DataFrame.sample)[0:20]
# df = multi_sequences_to_pango(get_seq_from_db('gisaid.fasta.db',samples.taxon.tolist()),threads=6)
# #%%
# df.join(lin,on='strain',how='left')
# # %%

# %%

# Generate metrics for certainty based on diff_second, and maybe even diff_third
# Predict if true or not based on likelihood,diff_second, diff_third etc -> nice separator
# joined[joined.diff_true.isna() & (joined.diff_second > -50)]['diff_second'].hist(bins=25)
# joined[joined.diff_true.notna()]['diff_second'].hist()