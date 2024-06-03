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

sc.logging.print_header()

adata = sc.read_10x_mtx(
    '',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
adata

count = scv.DataFrame(adata.X)
count.index = adata.obs.index
count.columns = adata.var.index
count

coor = pd.read_csv('/barcode_x_y.csv', index_col=0, header=None)
#coor = coor.drop(labels=0)
coor.index =coor[2]
coor[1] = 1
coor = coor.loc[:,[1,3,4]]
coor.columns = [1,2,3]
coor[4] = coor[2]
coor[5] = coor[3]
coor

s = count.index.intersection(coor.index)
count= count.loc[s]
coor = coor.loc[count.index]
coor

coord1 = coor.loc[:, 0:3]
coord1.columns = ['in_tissue','array_row', 'array_col']
coord2 = coor.loc[:, [4,5]]
coord2 = coord2.values

adata=adata[s,:]
adata.obs = coord1
adata.obsm['spatial'] = coord2
adata.var_names_make_unique()
adata

from PIL import Image
image = Image.open('')
image

spatial_key = "spatial"
library_id = "Sample_name"
adata.uns[spatial_key] = {library_id: {}}
adata.uns[spatial_key][library_id]["images"] = {}
adata.uns[spatial_key][library_id]["images"] = {"hires": image}
adata.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 19, "spot_diameter_fullres": 0.5}

sc.pl.spatial(adata, img_key='hires', color=['in_tissue'], 
              cmap='Spectral_r',spot_size=0.1, bw=False, alpha_img=1)

adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
sc.pl.highest_expr_genes(adata, n_top=20)

sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters")
sc.tl.tsne(adata)

plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)

fig, axes = plt.subplots(nrows=1, ncols=2, figsize = (8,4))
sc.pl.umap(adata, color=['clusters'], wspace=0.5, ax=axes[0],show = False)
sc.pl.spatial(adata, img_key='hires',color='clusters',spot_size=0.9,ax=axes[1],show=False)
plt.subplots_adjust(wspace=0.5)

adata.obs["x_array"]=coor[2]
adata.obs["y_array"]=coor[3]
adata.obs["x_pixel"]=coor[5]
adata.obs["y_pixel"]=coor[4]

#Set coordinates
x_array=adata.obs["x_array"].tolist()
y_array=adata.obs["y_array"].tolist()
x_pixel=adata.obs["x_pixel"].tolist()
y_pixel=adata.obs["y_pixel"].tolist()

#Read in hitology image
img=cv2.imread('')
#Test coordinates on the image
img_new=img.copy()
for i in range(len(x_pixel)):
    x=x_pixel[i]
    y=y_pixel[i]
    img_new[int(x-20):int(x+20), int(y-20):int(y+20),:]=0

cv2.imwrite('', img_new)

#Calculate adjacent matrix
s=1
b=49
adj=spg.calculate_adj_matrix(x=x_pixel,y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=img, beta=b, alpha=s, histology=True)
#If histlogy image is not available, SpaGCN can calculate the adjacent matrix using the fnction below
#adj=calculate_adj_matrix(x=x_pixel,y=y_pixel, histology=False)
#np.savetxt('', adj, delimiter=',')
#adj=np.loadtxt('', delimiter=',')
adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)
#Normalize and take log for UMI
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

#Find the l value given p
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

#For this toy data, we set the number of clusters=7 since this tissue has 7 layers
n_clusters=15
#Set seed
r_seed=t_seed=n_seed=100
#Search for suitable resolution
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

clf=spg.SpaGCN()
clf.set_l(l)
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
#Do cluster refinement(optional)
#shape="hexagon" for Visium data, "square" for ST data.
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
adata.obs["refined_pred"]=refined_pred
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')

#adata=sc.read("./sample_results/results.h5ad")
#adata.obs should contain two columns for x_pixel and y_pixel
#Set colors used
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
#Plot spatial domains
domains="pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig("", dpi=600)

#Plot refined spatial domains
domains="refined_pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig("", dpi=600)

plot_color=['#1F77B4FF','#FF7F0EFF','#2CA02CFF','#D62728FF','#9467BDFF','#8C564BFF','#E377C2FF','#7F7F7FFF',
            '#BCBD22FF','#17BECFFF','#AEC7E8FF','#FFBB78FF','#98DF8AFF','#FF9896FF','#C5B0D5FF','#C49C94FF',
            '#F7B6D2FF','#C7C7C7FF', '#DBDB8DFF','#9EDAE5FF',"#F56867","#FEB915","#C798EE","#59BE86","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C"]
r_seed=t_seed=n_seed=100

for l in [0.3, 0.5, 0.7, 0.9]:# 越多cluster越少
    for res in [0.5, 1.0, 1.5, 2.0]:#越大 cluster越多 默认1
        clf=spg.SpaGCN()
        clf.set_l(l)
        #Set seed
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)
        #Run
        clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob=clf.predict()
        adata.obs["pred"]= y_pred
        adata.obs["pred"]=adata.obs["pred"].astype('category')
        #Do cluster refinement(optional)
        #shape="hexagon" for Visium data, "square" for ST data.
        adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
        refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
        adata.obs["refined_pred"]=refined_pred
        adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
        domains="pred"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        plt.savefig("" % (l,res), dpi=600)
        domains="refined_pred"
        num_celltype=len(adata.obs[domains].unique())
        adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
        ax=sc.pl.scatter(adata,alpha=1,x="y_pixel",y="x_pixel",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
        ax.set_aspect('equal', 'box')
        ax.axes.invert_yaxis()
        plt.savefig("" % (l,res), dpi=600)  
#保存UMAP坐标
cord=pd.DataFrame(data=adata.obs[domains])
cord.to_csv("" % (l,res)) 

del adata.uns[spatial_key][library_id]["images"]
try:
    del adata.uns[spatial_key][library_id]["images"]
except KeyError:
    pass
adata.write("/Sample_name_spaGCN.h5ad")

