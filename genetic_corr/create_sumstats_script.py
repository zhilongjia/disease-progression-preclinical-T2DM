import os

import pandas as pd
'''
dos (Windows) file convert to unix: 
    sed -i 's/\r$//' filename
    or in vim: set ff=unix
    or use dos2unix: dos2unix filename
    or use tr: tr -d '\r' < inputfile > outputfile
    or use awk: awk '{ sub("\r$", ""); print }' inputfile > outputfile
'''

# path to reference data for windows subsystem for linux (WSL), D: drive to /mnt/d
ref_data_path = '/mnt/d/OneDrive/work/Data/eur_w_ld_chr'

# sumstats for brain disorder
df = pd.read_csv('data/disease_gwas_ids.csv')
sample_sizes = df['Sample_size'].tolist()
phenos = df['pheno'].tolist()
gwas_files = df['fullsummary'].tolist()

out_path = 'data/sumstats/disease'
os.makedirs(out_path, exist_ok=True)

with open('data/sumstats/sumstats.sh', 'w') as f:
    f.write('# munge sumstats for disease\n')
    for i in range(len(gwas_files)):
        gwas_file = gwas_files[i]
        pheno = phenos[i]
        # use ldsc to create sumstats
        f.write(f'munge_sumstats --sumstats ../gwas_summary/disease/{gwas_file} --N {sample_sizes[i]}'
                f' --out disease/{pheno}'
                f' --merge-alleles {ref_data_path}/w_hm3.snplist\n')
    f.write('\n')

# sumstats for t2dm
df_diabetes = pd.read_csv('data/t2dm_gwas_files.csv')
sample_sizes = df_diabetes['Sample_size'].tolist()
gwas_files = df_diabetes['gwas_id'].tolist()
phenos = df_diabetes['pheno'].tolist()

out_path = 'data/sumstats/t2dm'
os.makedirs(out_path, exist_ok=True)

# append to the same file
with open('data/sumstats/sumstats.sh', 'a') as f:
    f.write('# munge sumstats for t2dm\n')
    for i in range(len(gwas_files)):
        gwas_file = gwas_files[i]
        pheno = phenos[i]
        # use ldsc to create sumstats
        f.write(f'munge_sumstats --sumstats ../gwas_summary/t2dm/{gwas_file} --N {sample_sizes[i]}'
                f' --out t2dm/{pheno}'
                f' --merge-alleles {ref_data_path}/w_hm3.snplist\n')
    f.write('\n')


# sumstats for cognitive
df_cognitive = pd.read_csv('data/cognitive_gwas_files.csv')
sample_sizes = df_cognitive['Sample_size'].tolist()
gwas_files = df_cognitive['gwas_id'].tolist()
phenos = df_cognitive['pheno'].tolist()

out_path = 'data/sumstats/cognitive'
os.makedirs(out_path, exist_ok=True)

with open('data/sumstats/sumstats.sh', 'a') as f:
    f.write('# munge sumstats for cognitive\n')
    for i in range(len(gwas_files)):
        gwas_file = gwas_files[i]
        pheno = phenos[i]
        # use ldsc to create sumstats
        f.write(f'munge_sumstats --sumstats ../gwas_summary/cognitive/{gwas_file} --N {sample_sizes[i]}'
                f' --out cognitive/{pheno}'
                f' --merge-alleles {ref_data_path}/w_hm3.snplist\n')
    f.write('echo "All done!"\n')
    f.close()


# IDP
df_idp = pd.read_csv('data/brain_idps_587used.csv')
idp = df_idp['Pheno'].tolist()
sample_sizes = df_idp['N(all)'].tolist()
# idx to 4 digits
idp = [f'{i:04d}' for i in idp]

# path to sumstats
gwas_path_idp = '/mnt/e/data/idp'
os.makedirs(out_path, exist_ok=True)

# create script
with open('sumstats_idp.sh', 'w') as f:
    f.write('# munge sumstats for IDP\n')
    for i in range(len(idp)):
        gwas_file = f'{idp[i]}_idp.txt'
        pheno = idp[i]
        # use ldsc to create sumstats
        f.write(f'munge_sumstats --sumstats {gwas_path_idp}/{gwas_file} --N {sample_sizes[i]}'
                f' --out data/sumstats/idp/{pheno}'
                f' --merge-alleles {ref_data_path}/w_hm3.snplist\n')
    f.write('echo "All done!"\n')
    f.close()

# sumstats for CMR
df = pd.read_csv('data/Bai82_names_ukb_v2.csv')
phenos = df['ID'].tolist()

gwas_path_cmr = '/mnt/e/download/CMR/unzipped'

with open('sumstats_cmr.sh', 'w') as f:
    f.write('# munge sumstats for IDP\n')
    for i in range(len(phenos)):
        gwas_file = f'{idp[i]}_idp.txt'
        pheno = phenos[i]
        # use ldsc to create sumstats
        f.write(f'munge_sumstats --sumstats {gwas_path_cmr}/ukb_phase1to3_heart_may_2022_pheno{pheno}.fastGWA'
                f' --N-col N --frq AF1'
                f' --out data/sumstats/cmr/{pheno}'
                f' --merge-alleles {ref_data_path}/w_hm3.snplist\n')
    f.write('echo "All done!"\n')
    f.close()