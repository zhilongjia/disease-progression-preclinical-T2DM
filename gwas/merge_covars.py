import pandas as pd


def main(args):
    subtype = args.subtype

    # pca results
    pca_results = pd.read_csv('pca.eigenvec', sep=' ', header=None)
    pca_results.columns = ['FID', 'IID'] + ['PC' + str(i) for i in range(1, 11)]

    # phenotype
    df_pheno = pd.read_csv(f'{subtype}.pheno', sep=' ')

    # merge
    df = pd.merge(df_pheno, pca_results, on=['FID', 'IID'])

    # save
    df.to_csv(f'{subtype}.phenos', sep=' ', index=False, na_rep='NA')

    print(df['t2dm'].value_counts())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--subtype', type=str, default='subtype1')
    args = parser.parse_args()
    main(args)
