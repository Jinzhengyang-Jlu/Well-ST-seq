import os,csv,re
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
#In order to read in image data, we need to install some package. Here we recommend package "opencv"
#inatll opencv in python
#!pip3 install opencv-python
import cv2
import scvelo as scv

test = sc.read_h5ad("Sample_name.h5ad")
from PIL import Image
image = Image.open('')
spatial_key = "spatial"
library_id = "tissue_Sample_name"
test.uns[spatial_key] = {library_id: {}}
test.uns[spatial_key][library_id]["images"] = {}
test.uns[spatial_key][library_id]["images"] = {"hires": image}
test.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 19, "spot_diameter_fullres": 0.5}
adata = test

fig, axes = plt.subplots(nrows=1, ncols=3, figsize = (20,4))
myfig = plt.gcf()
sc.pl.umap(adata, color=['refined_pred'], wspace=0.5, ax=axes[0],show = False)
sc.pl.spatial(adata, img_key='hires',color='refined_pred',spot_size=0.9,ax=axes[1],show=False)
sc.pl.spatial(adata, img_key='hires', color='in_tissue',ax=axes[2],
              cmap='Spectral_r',spot_size=0, bw=False, alpha_img=1)
plt.subplots_adjust(wspace=0.5)

plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
plot_color=['#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF','#E377C2FF','#7F7F7FFF',
            '#BCBD22FF','#17BECFFF','#AEC7E8FF','#FFBB78FF','#98DF8AFF','#FF9896FF','#C5B0D5FF','#C49C94FF',
            '#F7B6D2FF','#C7C7C7FF', '#DBDB8DFF','#9EDAE5FF',"#F56867","#FEB915","#C798EE","#59BE86","#6D1A9C",
            "#15821E","#3A84E6","#997273","#787878","#DB4C6C","#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1"]
domains="refined_pred"
plt.figure(figsize=(18, 6))
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()

#rename
old_to_new = {}
adata.obs['annotation'] = (
adata.obs['refined_pred']
.map(old_to_new)
.astype('category')
)

fig, axes = plt.subplots(nrows=1, ncols=3, figsize = (20,4))
myfig = plt.gcf()
sc.pl.umap(adata, color=['annotation'], wspace=0.5, ax=axes[0],show = False)
sc.pl.spatial(adata, img_key='hires',color='annotation',spot_size=0.9,ax=axes[1],show=False)
sc.pl.spatial(adata, img_key='hires', color='in_tissue',ax=axes[2],
              cmap='Spectral_r',spot_size=0, bw=False, alpha_img=1)
plt.subplots_adjust(wspace=0.5)

#heatmap
marker_genes_dict = {'cluster': ['gene_name'],
                    }
sc.pl.dotplot(adata, marker_genes_dict, groupby='annotation')
sc.pl.matrixplot(adata, marker_genes_dict, groupby='annotation', dendrogram=False, standard_scale='var')

del adata.uns[spatial_key][library_id]["images"]
try:
    del adata.uns[spatial_key][library_id]["images"]
except KeyError:
    pass
adata.write("/Sample_name.h5ad")

#top 10 gene
test = sc.read_h5ad("/Sample_name_anno.h5ad")
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

fig, axes = plt.subplots(nrows=1, ncols=3, figsize = (20,4))
myfig = plt.gcf()
sc.pl.umap(adata, color=['annotation'], wspace=0.5, ax=axes[0],show = False)
sc.pl.spatial(adata, img_key='hires',color='annotation',spot_size=0.9,ax=axes[1],show=False)
sc.pl.spatial(adata, img_key='hires', color='in_tissue',ax=axes[2],
              cmap='Spectral_r',spot_size=0, bw=False, alpha_img=1)
plt.subplots_adjust(wspace=0.5)

sc.pl.umap(adata, color=['annotation'], wspace=0.5,show = False)
sc.pl.umap(adata, color=['annotation'], wspace=0.5,show = False,save = 'Sample_name_annotation_umap.pdf')
sc.pl.spatial(adata, img_key='hires',color='annotation',spot_size=0.9,show=False)
sc.pl.spatial(adata, img_key='hires',color='annotation',spot_size=0.9,show=False,save = 'Sample_name_annotation_spatial.pdf')
sc.pl.spatial(adata, img_key='hires', color='in_tissue',
              cmap='Spectral_r',spot_size=0, bw=False, alpha_img=1)
sc.pl.spatial(adata, img_key='hires', color='in_tissue',
              cmap='Spectral_r',spot_size=0, bw=False, alpha_img=1,save = 'Sample_name_tissue.pdf')

sc.tl.rank_genes_groups(adata, 'annotation')
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)

top20 = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(20)
top20.to_csv("/Sample_name_Anno_DEGs_top20.csv",index=False,sep=',')
allmk = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(100)
allmk.to_csv("/Sample_name_Anno_DEGs_all.csv",index=False,sep=',')

#top heatmap
top5 = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
marker_genes_dict = {'cluster': list(top5['gene_name']),
                    }
sc.pl.dotplot(adata, marker_genes_dict, groupby='annotation')
sc.pl.matrixplot(adata, marker_genes_dict, groupby='annotation', dendrogram=False, standard_scale='var')

#picked heatmap
top5 =  pd.read_csv('/Sample_name_top20_final.csv', index_col=None, header=0)
marker_genes_dict = {'cluster': list(top5['cluster'].str.capitalize())[0:],
                    }
sc.pl.matrixplot(adata, marker_genes_dict, groupby='annotation', dendrogram=False, standard_scale='var',save = 'Sample_name_final_picked_gene_heatmap.pdf')
allgene = pd.DataFrame(adata.var_names)
allgene.to_csv("/Sample_name/all_gene.csv",index=False,sep=',')

#featureplot
for k in marker_genes_dict:
    print(marker_genes_dict[k])
    
genename = ['']
sc.set_figure_params(facecolor='white', figsize=(3,3),dpi_save=300, dpi=100)
sc.pl.spatial(adata, img_key='hires', color=genename, cmap='Spectral_r', ncols=5,
             basis='spatial', img=None, size=2.3, spot_size=None, bw=False, alpha_img=0.5,save = 'Sample_name_picked_gene_featureplot.pdf')

#UMI feature
sc.pl.spatial(adata, img_key='hires', color='total_counts', cmap='Spectral_r', 
             basis='spatial', img=None, size=2.3, spot_size=None, bw=False, alpha_img=0.5,save = 'Sample_name_total_counts.pdf')
sc.pl.spatial(adata, img_key='hires', color='log1p_total_counts', cmap='Spectral_r', 
             basis='spatial', img=None, size=2.3, spot_size=None, bw=False, alpha_img=0.5,save = 'Sample_name_log1p_total_counts.pdf')
#UMI gene vln
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
import scvelo as scv
from brokenaxes import brokenaxes

df1 = pd.DataFrame(adata.obs['total_counts'])
df2 = pd.DataFrame(adata.obs['n_genes_by_counts'])

df = df1
df['n_genes_by_counts'] = df2['n_genes_by_counts']
df.rename(columns={'total_counts':'UMI'},inplace=True)
df.rename(columns={'n_genes_by_counts':'Genes'},inplace=True)
print(df)

import seaborn as sns
fig, ax = plt.subplots(figsize=(6,6))
sns.violinplot(data=df.iloc[:,0:2])
bax = brokenaxes(#xlims=((0, 10), (11, 20)), #设置x轴裂口范围
                 ylims=((0, 12000), (0.4, 2)), #设置y轴裂口范围
                 hspace=0.25,#y轴裂口宽度
                 wspace=0.2,#x轴裂口宽度                 
                 despine=False,#是否y轴只显示一个裂口
                 diag_color='r',#裂口斜线颜色                
                
                )
bax.plot(pts)
plt.show()
df.to_csv("/Sample_name_UMI_gene.csv",index=True,sep=',')

#save data
adata.obs
dfsave = pd.DataFrame(adata.obs)
dfsave.to_csv('/Sample_name_metadata.csv')