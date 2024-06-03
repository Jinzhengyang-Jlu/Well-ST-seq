import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#sc.logging.print_header()
#sc.set_figure_params(facecolor='white', figsize=(3,3), dpi_save=600, dpi=300)
#sc.settings.verbosity = 3

adata = sc.read_10x_mtx(
    '/filtered_feature_bc_matrix/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
adata

from PIL import Image
image = Image.open('')

coor = pd.read_csv('barcode_x_y.csv', index_col=0, header=0)
coor.index =coor['barcode']
coor.columns = ['in_tissue','array_row', 'array_col']
coor.in_tissue = 1
coor

adata = adata[adata.obs.index.isin(coor.index)]
coor = coor.loc[adata.obs.index,:]
adata.obs = coor
adata.obsm['spatial'] = coor.loc[:,['array_row','array_col']].values
adata

spatial_key = "spatial"
library_id = "tissue"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 103, "spot_diameter_fullres": 20}

adata.var_names_make_unique()
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    bins=60,
    ax=axs[3],
)

sc.pl.spatial(adata, img_key='hires', color=['n_genes_by_counts','total_counts','log1p_n_genes_by_counts','log1p_total_counts'], 
              cmap='Spectral_r', ncols=2, basis='spatial', img=None, size=0.2, spot_size=5, bw=False, 
              alpha_img=1)

sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.tl.leiden(adata, key_added="clusters",resolution=1)

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4,palette="tab20")

sc.tl.rank_genes_groups(adata, "clusters", method="wilcoxon")

result = adata.uns["rank_genes_groups"]
groups = result["names"].dtype.names
markers = pd.DataFrame({
        group + "_" + key: result[key][group]
        for group in groups
        for key in ["names", "pvals",'logfoldchanges']}
)
markers.head(10)
markers.to_csv('markers.csv')

marker_genes = ['']
sc.pl.dotplot(adata, marker_genes, groupby="clusters");

adata.uns['spatial'] = []
adata.write('Sample_name.h5ad')
adata = sc.read_h5ad('Sample_name.h5ad')
meta = pd.read_csv('Sample_name_meta.csv', index_col=0)
adata.obs['Anno'] = meta['anno']
adata

from PIL import Image
image = Image.open('Sample_name.png')
spatial_key = "spatial"
library_id = "tissue_Sample_name"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 106.5, "spot_diameter_fullres": 20}

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.spatial(adata, img_key="hires", color="Anno", size=0.04, alpha_img=1, save='_Sample_name_Anno.pdf'
             )

sc.pl.spatial(adata, img_key='hires',color=[''],
              cmap='Spectral_r', ncols=3,
              basis='spatial', img=None, size=0.8, spot_size=1, bw=False, alpha_img=0.8,
             show=False,save='_Sample_name_genesExpr.pdf'
             )
genes = ['']
for gene in genes:
    sc.pl.spatial(adata, img_key='hires',color=gene,
                  cmap='Spectral_r', ncols=1,
              basis='spatial', img=None, size=0.8, bw=False, alpha_img=0.8,
              spot_size=1,show=False,save=f'_Sample_name_{gene}.pdf'
             )