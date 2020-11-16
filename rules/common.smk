def is_single_end(sample, rep):
    return pd.isnull(units.loc[(sample, rep), 'fq2'])
