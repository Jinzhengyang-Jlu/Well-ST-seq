import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
import scvi
import scvelo as scv
from PIL import Image

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"

import torch 
from torch import nn 

# 1，查看gpu信息
if_cuda = torch.cuda.is_available()
print("if_cuda=",if_cuda)

gpu_count = torch.cuda.device_count()
print("gpu_count=",gpu_count)

results_folder = ''

# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

adata = sc.read_10x_mtx(
    '',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
count = scv.DataFrame(adata.X)
count.index = adata.obs.index
count.columns = adata.var.index
coor = pd.read_csv('barcode_x_y.csv', index_col=0, header=None)
#coor = coor.drop(labels=0)
coor.index =coor[2]
coor[1] = 1
coor = coor.loc[:,[1,3,4]]
coor.columns = [1,2,3]
coor[4] = coor[2]
coor[5] = coor[3]
s = count.index.intersection(coor.index)
count= count.loc[s]
coor = coor.loc[count.index]
coord1 = coor.loc[:, 0:3]
coord1.columns = ['in_tissue','array_row', 'array_col']
coord2 = coor.loc[:, [4,5]]
coord2 = coord2.values
adata=adata[s,:]
adata.obs = coord1
adata.obsm['spatial'] = coord2
adata.var_names_make_unique()

image = Image.open('')
spatial_key = "spatial"
library_id = "tissue_Sample_name"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 19, "spot_diameter_fullres": 0.5}
adata

adata.var['SYMBOL'] = adata.var_names
#adata.var.set_index('gene_ids', drop=True, inplace=True)
# find mitochondria-encoded (MT) genes
adata.var['MT_gene'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata.obsm['MT'] = adata[:, adata.var['MT_gene'].values].X.toarray()
adata = adata[:, ~adata.var['MT_gene'].values]
adata

sc.pl.spatial(adata, img_key='hires', color=['in_tissue'], 
              cmap='Spectral_r',spot_size=0.1, bw=False, alpha_img=1)

adata_ref = sc.read('scRNA_merge/mouse_data/mouseunlog2.h5ad')
adata_ref.var['SYMBOL'] = adata_ref.var.index
# rename 'GeneID-2' as necessary for your data
#adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)

# delete unnecessary raw slot (to be removed in a future version of the tutorial)
#del adata_ref.raw
adata_ref

scv.DataFrame(adata_ref.X)

sc.settings.set_figure_params(figsize=(5, 5))
sc.pl.umap(adata_ref,  color=['orig.ident','celltype'], 
             ncols=1,gene_symbols='SYMBOL',size=1.0)

adata_ref.var['SYMBOL'] = adata_ref.var.index
# rename 'GeneID-2' as necessary for your data
#adata_ref.var.set_index('GeneID-2', drop=True, inplace=True)

# delete unnecessary raw slot (to be removed in a future version of the tutorial)
del adata_ref.raw

from cell2location.utils.filtering import filter_genes

sc.set_figure_params(facecolor='white', figsize=(7,7),dpi_save=300, dpi=100)
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# filter the object
adata_ref = adata_ref[:, selected].copy()

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key='orig.ident',
                        # cell type, covariate used for constructing signatures
                        labels_key='celltype',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=None
                       )
# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# view anndata_setup as a sanity check
mod.view_anndata_setup()

mod.train(max_epochs=250, use_gpu=True)

mod.plot_history(20)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

sc.set_figure_params(facecolor='white', figsize=(6,6),dpi_save=300, dpi=100)
mod.plot_QC()

adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata.var_names, inf_aver.index)
adata = adata[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
adata.obs['sample'] = ''
# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata, batch_key='sample')

# create and train the model
mod = cell2location.models.Cell2location(
    adata, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=30,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()

mod.train(max_epochs=10000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True)

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata = mod.export_posterior(
    adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
#adata.write(adata_file)
adata_file

adata_file = f"{run_name}/sp.h5ad"
#adata = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata)

adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
sc.pl.spatial(adata, cmap='magma',
                  # show first 8 cell types
                  color=[],
                  ncols=5, size=2,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.2'
                 )
adata_save = adata
del adata_save.uns['spatial']
adata_file = f"{run_name}/sp.h5ad"
adata_save.write(adata_file)
adata_file
adata_file = f"{run_name}/sp.h5ad"
adata_save = sc.read_h5ad(adata_file)
adata_save

image = Image.open('')
spatial_key = "spatial"
library_id = "tissue_Sample_name"
adata_save.uns[spatial_key] = {library_id: {}}
adata_save.uns[spatial_key][library_id]["images"] = {}
adata_save.uns[spatial_key][library_id]["images"] = {"hires": image}
adata_save.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 19, "spot_diameter_fullres": 0.5}
adata_save

# Now we use cell2location plotter that allows showing multiple cell types in one panel
from cell2location.plt import plot_spatial

# select up to 6 clusters
clust_labels = []
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

#slide = select_slide(adata_save, 'jilin_st_10um')

with mpl.rc_context({'figure.figsize': (11, 11)}):
    fig = plot_spatial(
        adata=adata_save,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.992,
        # size of locations (adjust depending on figure size)
        circle_diameter=8,
        colorbar_position='right'
    )