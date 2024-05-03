import os

import pandas as pd

'''
genetic correlation between t2dm and brain health (Brain disorder, cognitive)
using LDSC
# ldsc
# --rg /path/to/pheno1.sumstats.gz,/path/to/pheno2.sumstats.gz 
# --ref-ld-chr /mnt/d/OneDrive/work/Data/eur_w_ld_chr
# --w-ld-chr /mnt/d/OneDrive/work/Data/eur_w_ld_chr
# --out /path/to/result/pheno1_pheno2_h2
'''

# path to reference data for windows subsystem for linux (WSL), D: drive to /mnt/d
ref_data_path = '/mnt/d/OneDrive/work/Data/eur_w_ld_chr/'

# pheno
df_t2dm = pd.read_csv('data/t2dm_gwas_files.csv')
df_brain = pd.read_csv('data/disease_gwas_ids.csv')
df_cognitive = pd.read_csv('data/cognitive_gwas_files.csv')
df_idp = pd.read_csv('data/brain_idps_587used.csv')
df_cmr = pd.read_csv('data/Bai82_names_ukb_v2.csv')

t2dm = df_t2dm['pheno'].tolist()
brain = df_brain['pheno'].tolist()
cognitive = df_cognitive['pheno'].tolist()
idp = df_idp['Pheno'].tolist()
idp = [f'{i:04d}' for i in idp]
cmrs = df_cmr['ID'].tolist()

# path to sumstats
sumstats_t2dm = 'data/sumstats/t2dm'
sumstats_brain = 'data/sumstats/disease'
sumstats_cognitive = 'data/sumstats/cognitive'
sumstats_idp = 'data/sumstats/idp'
sumstats_cmr = 'data/sumstats/cmr'

# result path
out_path = 'results/ldsc'
out_path_bd = f'{out_path}/t2dm_disease'
out_path_cog = f'{out_path}/t2dm_cognitive'
out_path_idp = f'{out_path}/t2dm_idp'
out_path_cmr = f'{out_path}/t2dm_cmr'

os.makedirs(out_path_bd, exist_ok=True)
os.makedirs(out_path_cog, exist_ok=True)
os.makedirs(out_path_idp, exist_ok=True)
os.makedirs(out_path_cmr, exist_ok=True)

# create script
with open('ldsc.sh', 'w') as f:
    f.write('# genetic correlation between t2dm and Diseases\n')
    for pheno1 in t2dm:
        for pheno2 in brain:
            f.write(f'ldsc --rg {sumstats_t2dm}/{pheno1}.sumstats.gz,{sumstats_brain}/{pheno2}.sumstats.gz'
                    f' --ref-ld-chr {ref_data_path}'
                    f' --w-ld-chr {ref_data_path}'
                    f' --out {out_path}/{pheno1}_{pheno2}_h2\n')
            # f.write('\n')
    f.write('# genetic correlation between t2dm and cognitive\n')
    for pheno1 in t2dm:
        for pheno2 in cognitive:
            f.write(f'ldsc --rg {sumstats_t2dm}/{pheno1}.sumstats.gz,{sumstats_cognitive}/{pheno2}.sumstats.gz'
                    f' --ref-ld-chr {ref_data_path}'
                    f' --w-ld-chr {ref_data_path}'
                    f' --out {out_path}/{pheno1}_{pheno2}_h2\n')
            # f.write('\n')
    f.write('echo "ALL done!\n"')

f.close()

# create script for idp
with open('ldsc_idp.sh', 'w') as f:
    f.write('# genetic correlation between t2dm and IDP\n')
    for pheno1 in t2dm:
        for pheno2 in idp:
            f.write(f'ldsc --rg {sumstats_t2dm}/{pheno1}.sumstats.gz,{sumstats_idp}/{pheno2}.sumstats.gz'
                    f' --ref-ld-chr {ref_data_path}'
                    f' --w-ld-chr {ref_data_path}'
                    f' --out {out_path}/t2dm_idp/{pheno1}_{pheno2}_h2\n')
    f.write('echo "ALL done!\n"')

with open('ldsc_cmr.sh', 'w') as f:
    f.write('# genetic correlation between t2dm and CMR\n')
    for pheno1 in t2dm:
        for pheno2 in cmrs:
            f.write(f'ldsc --rg {sumstats_t2dm}/{pheno1}.sumstats.gz,{sumstats_cmr}/{pheno2}.sumstats.gz'
                    f' --ref-ld-chr {ref_data_path}'
                    f' --w-ld-chr {ref_data_path}'
                    f' --out {out_path}/t2dm_cmr/{pheno1}_{pheno2}_h2\n')
    f.write('echo "ALL done!\n"')