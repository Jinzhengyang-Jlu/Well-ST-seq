import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import interpolate

test = sc.read_h5ad("Sample_name.h5ad")
test.uns['log1p']["base"] = None
from PIL import Image
image = Image.open('')
spatial_key = "spatial"
library_id = "tissue_Sample_name"
test.uns[spatial_key] = {library_id: {}}
test.uns[spatial_key][library_id]["images"] = {}
test.uns[spatial_key][library_id]["images"] = {"hires": image}
test.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 19, "spot_diameter_fullres": 0.5}

adata = test
sc.pl.spatial(adata, img_key='hires',color='annotation',spot_size=0.9,show=False)

plot_color=['#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF','#E377C2FF','#7F7F7FFF','#28757D','#17BECFFF']
domains="annotation"
plt.figure(figsize=(18, 6))
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()

sc.pl.umap(adata, color=['annotation'], wspace=0.5,show = False)
sc.pl.umap(adata, color=['annotation'], wspace=0.5,show = False,save = 'Sample_name_umap.pdf')

sc.pl.spatial(adata, img_key='hires',color='annotation',spot_size=0.9,show=False)
sc.pl.spatial(adata, img_key='hires',color='annotation',spot_size=0.9,show=False,save = 'Sample_name_spatial.pdf')

import matplotlib as mpl
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

clist=['#B9C7DF','#E5B3AA','#F24B38']
# N 表示插值后的颜色数目
newcmp = LinearSegmentedColormap.from_list('hanyincolor',clist,N=256)
#aa = plt.imshow(np.random.randint(1,10,(5,5)),newcmp)
#plt.colorbar(aa)
newcmp
plt.register_cmap(cmap=newcmp)

sc.pl.spatial(adata, img_key='hires', color='total_counts', cmap='hanyincolor',#RdBu_r
             basis='spatial', img=None, spot_size=0.9, bw=False, alpha_img=0.8,save = 'Sample_name_total_counts.pdf')#
sc.pl.spatial(adata, img_key='hires', color='log1p_total_counts',cmap='hanyincolor',#RdBu_r
             basis='spatial', img=None, spot_size=0.9, bw=False, alpha_img=0.8,save = 'Sample_name_log_counts.pdf')#
sc.pl.spatial(adata, img_key='hires', color='n_genes_by_counts', cmap='hanyincolor',
             basis='spatial', img=None, spot_size=0.9, bw=False, alpha_img=0.8,save = 'Sample_name_n_genes_by_counts.pdf')#
sc.pl.spatial(adata, img_key='hires', color='log1p_n_genes_by_counts', cmap='hanyincolor',
             basis='spatial', img=None, spot_size=0.9, bw=False, alpha_img=0.8,save = 'Sample_name_n_genes_by_counts.pdf')#

#featureplot
genename = ['']
adata.X = np.array(pd.DataFrame(adata.X.todense()).replace(0,NaN))# 让0表达的点为透明

sc.set_figure_params(facecolor='white', figsize=(3,3),dpi_save=300, dpi=100)
sc.pl.spatial(adata, img_key='hires', color=genename, cmap='Spectral_r', ncols=5,
             basis='spatial', img=None, spot_size=0.9, bw=False, alpha_img=0.8,save = 'picked_gene_featureplot.pdf')#

#vlnplot
import matplotlib as mpl
from pylab import *
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,LinearSegmentedColormap

df1 = pd.DataFrame(log10(adata.obs['total_counts']))
df2 = pd.DataFrame(log10(adata.obs['n_genes_by_counts']))

df = df1
df['n_genes_by_counts'] = df2['n_genes_by_counts']
df.rename(columns={'total_counts':'UMI'},inplace=True)
df.rename(columns={'n_genes_by_counts':'Genes'},inplace=True)
print(df)

import seaborn as sns
import matplotlib.pyplot as plt
plt.figure(figsize=(8,5))  
sns.violinplot(data=df,inner = None)
plt.ylabel('log10 UMI/ log10 Genes')
#plt.show()
savefig("Sample_name_log10_vln.pdf",dpi=300)

