#%%
from create_matrix import *

#%%
matrix_dict = import_np_matrix('full_run.npz')
index = matrix_dict['index']
matrix = matrix_dict['array'].astype(np.float32)
matrix += 0.0001
matrix[:,:,5] = 0
matrix = matrix / matrix.sum(axis=2)[:,:,np.newaxis]

# %%
df = pd.DataFrame(index[1:,[0,2]], columns=index[0,[0,2]], index=index[1:,1].astype(str))
# %%
ref = get_seq_from_db('gisaid.fasta.db', ['Wuhan/Hu-1/2019'])[0]
# %%
ref_vec = seq_to_array(ref)
# %%
