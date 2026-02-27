import pandas as pd

uniprot_df = pd.read_csv('/home/angelis/thesis/reference_proteomes/active_sites.csv')
mafft_df = pd.read_csv('/home/angelis/thesis/reference_proteomes/active_sites_mafft.csv')

merged_df = pd.merge(uniprot_df, mafft_df, on='Accession', how='inner')

merged_df.to_csv('/home/angelis/thesis/reference_proteomes/final_active_sites_comparison.csv', index=False)

print('Merged data saved to final_active_sites_comparison.csv')