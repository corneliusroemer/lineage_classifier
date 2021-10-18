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
s1 = get_seq_from_db('gisaid.fasta.db', ['England/MILK-1526901/2021'])[0] 


# %%
ref_vec = seq_to_array(ref).astype(np.int8)
s1_vec = seq_to_array(s1).astype(bool)
ref_vec[:100]=0
ref_vec[-100:]=0
s1_vec[:100] = False
s1_vec[-100:] = False
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
# Ignore backmutations -> uninformative
# If it doesn't have the 99% ones, kick it back

# If it has all key mutations, bingo
# If it lacks one, can forgive if it has another of the somewhat common ones

# Always check against the parent

dp = matrix[df.loc['B.1.617.2'].pos].astype(np.float32)
dc = matrix[df.loc['AY.4'].pos]
diff = dc - dp
new_muts = (dc > 0.8).astype(bool) & (dp < 0.1).astype(bool)
muts_where = np.argwhere(new_muts)
for mut in muts_where:
    print(f"{mut}\n{dp[mut][0]}\n{dc[mut][0]}\n")

# %%

# Alternative
# Find common mutations in each lineage
# Score presence of mutation by each lineage
# Calculate number of hits vs. misses
# Make matrix that distinguishes nicely
# Present in query, present in lineage -> hit (somewhat good)
# Key: Not in query, sometimes in lineage -> miss (very bad) -> minimize
# Key: Not in query, 100% in lineage -> miss (very bad) -> minimize
# Present in query, never in lineage -> extra or other [shouldn't be too many]

# What could the metric be? c_miss * log(1-p-1/n) + c_extra * log(p+1/n) 
print((np.log(dp+0.001) * s1_vec).sum())
print((np.log(dc+0.001) * s1_vec).sum())

#%%
##
# Tabulate mutations in lineage
s = pd.Series([1,2,3],index=['a','b','c'])
# %%
df = pd.concat([df,s],keys=['a','b','c'])

# %%
df.index

# %%
