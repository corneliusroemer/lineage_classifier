## Plan
# Filter sequences down to designated pangos using tsv-filter
# Read line by line of aligned sequence
# Build 3D numpy array with counts at each position [lineage,position,sum of(A,C,G,T,-,N)]
#%%
import zstandard
import io
import pandas as pd
from collections import defaultdict
#%%
lineages = pd.read_csv('lineages.csv')
lineages.set_index('taxon',inplace=True)

#%%
def add_to_index(sequence,lineage,index):
    """Go through sequence base by base and add to index"""
    i = 0
    try:
        for step in 1,-1:
            while sequence[i]=='-':
                sequence[i] = 'N'
                i += step 
            i = len(sequence)-1
    except IndexError:
        print(i,sequence)

    for i, base in enumerate(sequence):
        if base in ['A','C','G','T','-']:
            index[lineage][i][base] += 1

    return index

#%%
index = defaultdict(lambda : defaultdict(lambda: defaultdict(int)))
path='filtered.fasta.zst'
with open(path, 'rb') as fh:
    dctx = zstandard.ZstdDecompressor()
    stream_reader = dctx.stream_reader(fh)
    text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
    counter = 0
    sequence_counter = 0
    lineage = None
    for line in text_stream:
        if line.startswith('>'):
            if lineage is not None:
                index = add_to_index(list(sequence),lineage,index)
            sequence = ""
            name = line.strip('>\n')
            # print(name)
            try: 
                lineage = lineages.loc[name]['lineage']
                counter += 1
            except KeyError:
                lineage = None
            sequence_counter += 1
        else:
            sequence += line.strip('\n')
        if sequence_counter % 10000 == 0:
            print(f"sequences: {sequence_counter}, pangos:{counter}")

        if counter > 100000:
            break



# %%
index['B.1.617.2'][2000]
# %%
index.keys()

# %%
result = []
storage = 0
for key in index.keys():
    storage += sys.getsizeof(index[key])
    sum = 0
    for i in index[key][2000].keys():
        sum += index[key][2000][i]
    result.append((key,sum))
sorted(result,key=lambda x: x[1])
print(storage)
# %%
