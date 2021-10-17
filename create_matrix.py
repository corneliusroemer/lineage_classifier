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

def get_seq_from_db(db_path, strain_names: list[str]) -> list[SeqRecord]:
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

def get_seq_not_in_db(db_path='gisaid.fasta.db') -> list[str]:
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

def designated_strains_from_pango(lineage_name,lineages_path='lineages.csv') -> list[str]:
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
    sequence = np.frombuffer(bytes(seq.seq),dtype=np.uint8).copy()

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

def aggregate_seqs(seqs: list[SeqRecord]) -> np.ndarray:
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

@ray.remote
def generate_series_from_pango(lineage_name,db_path='gisaid.fasta.db') -> pd.Series:
    """Generate a pandas series from a pango designation"""
    strains = designated_strains_from_pango(lineage_name)
    seqs = get_seq_from_db(db_path, strains)
    print(f"{lineage_name:<10}:{len(strains)-len(seqs):6d} missing of {len(strains):6d} designated strains")
    aggregate = aggregate_seqs(seqs)
    return aggregate

def generate_df_from_pangos(lineage_names,db_path='gisaid.fasta.db') -> dict:
    """Generate a pandas dataframe from a list of pango designations"""
    matrix = {'index': [('index','lineage','counts')]}
    matrix['array'] = np.zeros((len(lineage_names),29903,6), dtype=np.uint32)
    parallel_count = 0
    futures = []
    for i, lineage_name in enumerate(lineage_names):
        futures.append(generate_series_from_pango.remote(lineage_name,db_path))
        auto_garbage_collect()
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