import SPACEL
from SPACEL import Scube
import scanpy as sc
import pandas as pd
import numpy as np
import warnings
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

adata = sc.read_h5ad('Sample1.h5ad')
meta = pd.read_csv('Sample1_meta.csv', index_col=0)
adata.obs['Anno'] = meta['anno']
from PIL import Image
image = Image.open('Sample1.png')
spatial_key = "spatial"
library_id = "tissue_Sample1"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 106.5, "spot_diameter_fullres": 20}
ad1 = adata
ad1

adata = sc.read_h5ad('Sample2.h5ad')
meta = pd.read_csv('Sample2_meta.csv', index_col=0)
adata.obs['Anno'] = meta['anno']
from PIL import Image
image = Image.open('Sample2.png')
spatial_key = "spatial"
library_id = "tissue_Sample2"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 103, "spot_diameter_fullres": 20}
ad2 = adata
ad2

adata = sc.read_h5ad('Sample3.h5ad')
meta = pd.read_csv('Sample3_meta.csv', index_col=0)
adata.obs['Anno'] = meta['anno']
from PIL import Image
image = Image.open('Sample3.png')
spatial_key = "spatial"
library_id = "tissue_Sample3"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 107, "spot_diameter_fullres": 20}
ad3 = adata
ad3

adata = sc.read_h5ad('Sample4.h5ad')
meta = pd.read_csv('Sample4_meta.csv', index_col=0)
adata.obs['Anno'] = meta['anno']
from PIL import Image
image = Image.open('Sample4.png')
spatial_key = "spatial"
library_id = "tissue_Sample4"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 99, "spot_diameter_fullres": 20}
ad4 = adata
ad4

adata_list=[ad1,ad2,ad3,ad4]

Scube.align(adata_list,
      cluster_key='Anno',
      n_neighbors = 15,
      n_threads=10,
      p=1,
      write_loc_path='Scube_outputs/aligned_coordinates.csv'
    )

adata = sc.concat(adata_list)
adata.write(f'all_slices.h5ad')

adata.obsm['spatial_aligned'] = np.array(adata.obsm['spatial_aligned'])
sc.pl.embedding(adata,basis='spatial_aligned',color='Anno',show=False)

ad1.obsm['spatial_aligned'] = np.array(ad1.obsm['spatial_aligned'])
sc.pl.embedding(ad1,basis='spatial_aligned',color='Anno',show=False)
ad2.obsm['spatial_aligned'] = np.array(ad2.obsm['spatial_aligned'])
sc.pl.embedding(ad2,basis='spatial_aligned',color='Anno',show=False)
ad3.obsm['spatial_aligned'] = np.array(ad3.obsm['spatial_aligned'])
sc.pl.embedding(ad3,basis='spatial_aligned',color='Anno',show=False)
ad4.obsm['spatial_aligned'] = np.array(ad4.obsm['spatial_aligned'])
sc.pl.embedding(ad4,basis='spatial_aligned',color='Anno',show=False)

from SPACEL.setting import set_environ_seed
set_environ_seed()
from SPACEL import Scube
from SPACEL.Scube.utils_3d import create_mesh, smooth_mesh, sample_in_mesh, get_surface_colors, save_view_parameters, load_view_parameters
import scanpy as sc
import numpy as np
import pandas as pd
np.random.seed(42)

adata = sc.read_h5ad('./all_slices.h5ad')
adata
adata.obs
adata.obsm['spatial_aligned'] = np.array(adata.obsm['spatial_aligned'])
sc.pl.embedding(adata,basis='spatial_aligned',color='Anno',show=False)

st_ad = sc.read_h5ad('./all_slices_update.h5ad')
sc.pp.normalize_total(st_ad,target_sum=1e4)
sc.pp.log1p(st_ad)

# Normalized expression data (rows as spots/cells, columns as genes)
norm_expr = st_ad.to_df()
# 3D location (rows as spots/cells, columns as x,y,z axis coordinate)
loc = st_ad.obsm['spatial_aligned'].copy()

st_ad.obs.loc[st_ad.obs['sample']==7,'sample'] = 7
st_ad.obs.loc[st_ad.obs['sample']==8,'sample'] = 5
st_ad.obs.loc[st_ad.obs['sample']==9,'sample'] = 3
st_ad.obs.loc[st_ad.obs['sample']==10,'sample'] = 1
st_ad.obs
loc['Z'] = st_ad.obs['sample']
loc

color_mapping = {
    "Cavity": '#1F77B4FF',
    'Cb' : '#86BCB6',
    "DM": '#FF7F0EFF',
    "IE": '#D62728FF',
    'Mb': '#9467BDFF',
    'Mb/Hb VZ': '#8C564B',
    'PM' : '#E377C2FF',
    "PPH": '#7F7F7FFF',
    "PMH/MH": '#2CA02CFF'
}

st_ad.obs['anno_color'] = st_ad.obs['Anno'].map(color_mapping)

Scube.plot_3d(loc=loc.values,val=st_ad.obs['anno_color'],s=3,show=True,elev=35,azim=75,zlim=(0,7),save_dpi=300, save_path='3d.pdf')

color_mapping = {
    "Cavity": '#CDCDCD',
    'Cb' : '#CDCDCD',
    "DM": '#CDCDCD',
    "IE": '#CDCDCD',
    'Mb': '#CDCDCD',
    'Mb/Hb VZ': '#CDCDCD',
    'PM' : '#CDCDCD',
    "PPH": '#CDCDCD',
    "PMH/MH": '#2CA02CFF'
}

st_ad.obs['anno_color'] = st_ad.obs['Anno'].map(color_mapping)
Scube.plot_3d(loc=loc.values,val=st_ad.obs['anno_color'],s=3,show=True,elev=35,azim=75,zlim=(0,7),save_dpi=300, save_path='3d_PMH_MH.pdf')

