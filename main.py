#%%
from create_matrix import *

#%%
matrix_dict = import_np_matrix('full_run.npz')
index = matrix_dict['index']
matrix = matrix_dict['array'].astype(np.float32)
matrix += 0.0001
matrix[:,:,5] = 0

matrix = matrix / matrix.sum(axis=2)[:,:,np.newaxis]
matrix[:,:100,:] = 0
matrix[:,-100:,:] = 0

# %%
df = pd.DataFrame(index[1:,[0,2]].astype(int), columns=['pos','count'], index=index[1:,1].astype(str))
# %%
ref = get_seq_from_db('gisaid.fasta.db', ['Wuhan/Hu-1/2019'])[0]



# %%
ref_vec = seq_to_array(ref).astype(np.int8)
ref_vec[:100]=0
ref_vec[-100:]=0
# %%
test = seq_to_array(get_seq_from_db('gisaid.fasta.db', ['Denmark/DCGC-156462/2021'])[0]).astype(np.int8)

# %%
ay33 = matrix[df.loc['AY.33'].pos]

# %%
(ay33 - ref_vec).sum()
# %%
np.set_printoptions(precision=3,suppress=True)

# %%
d33 = pd.DataFrame(matrix[df.loc['AY.33'].pos])
d6 = pd.DataFrame(matrix[df.loc['AY.6'].pos])
dref = pd.DataFrame(ref_vec)
dd = pd.DataFrame(matrix[df.loc['AY.4'].pos])
## Extract mutation added in hierarchy: what's different?
## What is diversity within lineage?
## What is diversity in parent?
diff = d33 - dd
for row in range(len(dref)):
    diff_val = diff.iloc[row].max()
    if diff_val > 0.3:
        d_max = diff.iloc[row].argmax()
        a_max = d33.iloc[row].argmax()
        b_max = dd.iloc[row].argmax()
        c_max = dref.iloc[row].argmax()
        if b_max == c_max:
            if a_max == c_max:
                if d33.iloc[row][a_max] < dd.iloc[row][b_max]:
                    mut_type = 'descendent partial'
                else:
                    mut_type = 'sibling'

            else:
                if diff_val > 0.9:
                    mut_type = "defining"
                else:
                    if abs((dd.iloc[row]-dref.iloc[row])[d_max]) > 0.1: 
                        mut_type = "shared parental"
                    else:
                        mut_type = "descendent"
        else:
            if a_max == c_max:
                mut_type = 'back'
            else:
                mut_type = 'shared parental'
        print(f"{row:>5}: {d_max}/{diff.iloc[row].max():.2f}, child: {a_max}/{d33.iloc[row][a_max]:.2f}, parent: {b_max}/{dd.iloc[row][b_max]:.2f}, ref: {dref.iloc[row].argmax()} - {mut_type}")
# %%
# Where does it go in one rather than another direction?
# Losses for one, and losses for another