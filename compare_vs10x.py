import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(5, 5))
sc.settings.verbosity = 3

ad1 = ad_st = sc.read_visium('../../dataset/',
                       count_file= 'V1_Adult_Mouse_Brain_Coronal_Section_2_filtered_feature_bc_matrix.h5',
                       source_image_path='../../dataset/spatial/detected_tissue_image.jpg')

ad1.var_names_make_unique()
sc.pp.calculate_qc_metrics(
    ad1, percent_top=None, log1p=False, inplace=True
)
ad1

sc.pl.spatial(ad1, img_key='hires', color=['in_tissue'], cmap='Spectral_r',spot_size=0, bw=False, alpha_img=1)

add = ad1[((ad1.obs.array_row >15) &(ad1.obs.array_row <35) & (ad1.obs.array_col>60) &(ad1.obs.array_col<90)),:]
sc.pl.spatial(add, img_key='hires', color=['in_tissue'], cmap='Spectral_r',
              spot_size=0, bw=False, alpha_img=1, save='visium_tissue.pdf')
add

import scvelo as scv
count1 = scv.DataFrame(add.X)
count1.index = add.obs.index
count1.columns = add.var.index
count1

ad2 = sc.read_10x_mtx(
    'matrix',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)
ad2.var_names_make_unique()
ad2

import scvelo as scv
count2 = scv.DataFrame(ad2.X)
count2.index = ad2.obs.index
count2.columns = ad2.var.index
count2

coor = pd.read_csv('barcode_x_y.csv', index_col=0, header=None)
#coor = coor.drop(labels=0)
coor.index =coor[2]
coor[1] = 1
coor = coor.loc[:,[1,3,4]]
coor.columns = [1,2,3]
coor[4] = coor[2]
coor[5] = coor[3]
coor

s = count2.index.intersection(coor.index)
count2= count2.loc[s]
coor = coor.loc[count2.index]
coord1 = coor.loc[:, 0:3]
coord1.columns = ['in_tissue','array_row', 'array_col']
coord1

# new
coord2 = coor.loc[:, [4,5]]
coord2 = coord2.values
coord2

ad2=ad2[s,:]
ad2.obs = coord1
ad2.obsm['spatial'] = coord2

sc.pp.calculate_qc_metrics(
    ad2, percent_top=None, log1p=False, inplace=True
)
ad2

ad2.write_h5ad('10um.h5ad')

import scvelo as scv
count2 = scv.DataFrame(ad2.X)
count2.index = ad2.obs.index
count2.columns = ad2.var.index
count2

dff1 = count1.sum(axis=1)
dff2 = count2.sum(axis=1)

dff1 = pd.DataFrame({'group':['visium']*len(dff1),
                  'umi':dff1.values/round((17.5**(2) * 3.14),0)})
dff2 = pd.DataFrame({'group':['ours']*len(dff2),
                  'umi':dff2.values/100})
dff = pd.concat([dff1, dff2], axis=0)
dff

dff.umi = np.log(dff.umi)

from scipy.stats import mannwhitneyu
plt.figure(figsize=(8,9))
#sns.set(style=None)
sns.set_palette('pastel')
sns.violinplot(x='group', y='umi', data=dff, linewidth=2, palette='muted', saturation=0.75,inner = None)
#statistic, pvalue = mannwhitneyu(dff1.values, dff2.values)
#plt.txt(0.5, max(max(dff1),max(dff2)), f'p-value:{pvalue:0.4f}', ha='center',va='bottom',color='red')
plt.title('Violin Plot', fontsize=18)
plt.xlabel('Datasets', fontsize=16)
plt.ylabel(f'log1p(nUMI) per \u03BCm\u00B2', fontsize=16)
plt.grid(False)
plt.savefig('violin_compar.pdf', format='pdf', dpi=300, transparent=True)
plt.savefig('violin_compar.png', format='png', dpi=300, transparent=True)
plt.show()

ad2.obs

dff1 = pd.DataFrame({'group':['visium']*len(add.obs),
                  'umi':add.obs.n_genes_by_counts/round((17.5**(2) * 3.14),0)})
dff2 = pd.DataFrame({'group':['ours']*len(ad2.obs),
                  'umi':ad2.obs.n_genes_by_counts/100})

dff = pd.concat([dff1, dff2], axis=0)
dff

#from scipy.columnsts import mannwhitneyu
plt.figure(figsize=(8,9))
#sns.set(style=None)
sns.set_palette('pastel')
sns.violinplot(x='group', y='umi', data=dff, linewidth=2, palette='muted', saturation=0.75,inner = None)
#statistic, pvalue = mannwhitneyu(dff1.values, dff2.values)
#plt.txt(0.5, max(max(dff1),max(dff2)), f'p-value:{pvalue:0.4f}', ha='center',va='bottom',color='red')
plt.title('Violin Plot', fontsize=18)
plt.xlabel('Datasets', fontsize=16)
plt.ylabel(f'nums of genes per \u03BCm\u00B2', fontsize=16)
plt.grid(False)
plt.savefig('violin_compar_ngenes.pdf', format='pdf', dpi=300, transparent=True)
plt.savefig('violin_compar_ngenes.png', format='png', dpi=300, transparent=True)
plt.show()

dff1 = np.log1p(count1.sum(axis=0))
dff2 = np.log1p(count2.sum(axis=0))
dff1
(max(dff2)/max(dff1))*23

plt.figure(figsize=(10,6))
plt.hist(dff1, bins=23, alpha=0.7, label='visium',color='skyblue',edgecolor='skyblue',linewidth=1.2)
plt.hist(dff2, bins=25, alpha=0.7, label='ours',color='salmon',edgecolor='salmon',linewidth=1.2)
plt.legend()
plt.title('Histogram Plot', fontsize=18)
plt.xlabel('log1p(UMIs) per Feature', fontsize=16)
plt.ylabel('Number of Features', fontsize=16)
#plt.xlim([0,2000])
plt.grid(False)
#plt.grid(axis='y',linestyle='--',alpha=0.7)
plt.savefig('hist_compar.pdf', format='pdf', dpi=300, transparent=True)
plt.savefig('hist_compar.png', format='png', dpi=300, transparent=True)
#plt.grid(axis='y',linestyles='--',alpha=0.7)
plt.show()

ad2 = sc.read_h5ad('10um.h5ad')
ad2

from PIL import Image
image = Image.open('/tissue.png')
image

spatial_key = "spatial"
library_id = "tissue_10um"
ad2.uns[spatial_key] = {library_id: {}}
ad2.uns[spatial_key][library_id]["images"] = {}
ad2.uns[spatial_key][library_id]["images"] = {"hires": image}
ad2.uns[spatial_key][library_id]["scalefactors"] = {"tissue_hires_scalef": 19, "spot_diameter_fullres": 0.5}

ad2.var_names_make_unique()
sc.pl.spatial(ad2, img_key='hires', color=['in_tissue'], 
              cmap='Spectral_r',spot_size=0.1, bw=False, alpha_img=1, save='ours_tissue.pdf')

sc.pp.normalize_total(ad2, inplace=True)

sc.pl.spatial(ad2, img_key='hires', color=['gene_name'], cmap='Spectral_r', ncols=1,vmax=20,
             basis='spatial', img=None, size=0.8, spot_size=1, bw=False, alpha_img=0.8, save='ours_gene_name.pdf')

ad1 = ad_st = sc.read_visium('../../dataset/',
                       count_file= 'V1_Adult_Mouse_Brain_Coronal_Section_2_filtered_feature_bc_matrix.h5',
                       source_image_path='../../dataset/spatial/detected_tissue_image.jpg')

ad1.var_names_make_unique()
sc.pp.calculate_qc_metrics(
    ad1, percent_top=None, log1p=False, inplace=True
)
sc.pp.normalize_total(ad1, inplace=True)
add = ad1[((ad1.obs.array_row >15) &(ad1.obs.array_row <35) & (ad1.obs.array_col>60) &(ad1.obs.array_col<90)),:]
add

sc.pl.spatial(add, img_key='hires', color=['gene_name'], cmap='Spectral_r', ncols=2,
             basis='spatial', img=None, size=1.0, spot_size=None, bw=False, alpha_img=0.3, save='visium.pdf')

add

ad2
sc.pl.spatial(ad2, img_key='hires', color=['total_counts'], cmap='Spectral_r', ncols=2,vmax=30000,
             basis='spatial', img=None, size=0.8, spot_size=1, bw=False, alpha_img=1, save='ours_umi.pdf')

sc.pl.spatial(add, img_key='hires', color=['total_counts'], cmap='Spectral_r', ncols=2,vmax=30000,
             basis='spatial', img=None, size=1.0, spot_size=None, bw=False, alpha_img=0.3, save='visium_umi.pdf')