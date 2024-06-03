import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import scvi

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams["font.sans-serif"] = "Arial"

adata_file = f"/cell2location_map/sp.h5ad"
adata_save = sc.read_h5ad(adata_file)
from PIL import Image
image = Image.open('')
spatial_key = "spatial"
library_id = "tissue_Sample_name"
adata_save.uns[spatial_key] = {library_id: {}}
adata_save.uns[spatial_key][library_id]["images"] = {}
adata_save.uns[spatial_key][library_id]["images"] = {"hires": image}
adata_save.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 19, "spot_diameter_fullres": 0.5}
adata_save

adata = adata_save
adata_save.obs

dfsave = pd.DataFrame(adata_save.obs)
dfsave.to_csv('/scRNA_merge/pieplot/metadata.csv')

with mpl.rc_context({'axes.facecolor':  'white',#black white
                     'figure.figsize': [4.5, 5]}):

    sc.pl.spatial(adata, cmap='magma',
                  color=[],
                  ncols=2, size=1.6,
                  img_key='hires',
                  vmin=0, vmax='p99.2',
                         save='celltype_plot_spatial_split.pdf'
                 )
#plt.savefig('/celltype_plot_spatial_split.pdf')

from cell2location.plt import plot_spatial
# select up to clusters
clust_labels = ['']
clust_col = ['' + str(i) for i in clust_labels] # in case column names differ from labels

slide = adata#select_slide(adata_vis, 'V1_Human_Lymph_Node')

#flatui = ['#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF']

with mpl.rc_context({'figure.figsize': (15, 15)}):
    fig = plot_spatial(
        adata=slide,
        # labels to show on a plot
        color=clust_col, labels=clust_labels,
        show_img=True,
        reorder_cmap = [0,1,2,3,4,6],
        #image_cmap = 'PRGn',
        # 'fast' (white background) or 'dark_background'
        style='fast',
        # limit color scale at 99.2% quantile of cell abundance
        max_color_quantile=0.9,
        # size of locations (adjust depending on figure size)
        circle_diameter=8,
        colorbar_position='right',
        white_spacing = 10,
    )
fig.savefig('/celltype_plot_spatial.pdf')

adata_ref = sc.read('/scRNA_merge/mouse_data/mouseunlog.h5ad')#mouse.h5ad
adata_ref

sc.settings.set_figure_params(figsize=(5, 5))
sc.pl.umap(adata_ref, color='celltype',save='sc_celltype_umap.pdf')

adata_ref

sc.settings.set_figure_params(figsize=(5, 5))
sc.pl.umap(adata_ref, color='seurat_clusters')#,save='sc_celltype_umap.pdf'