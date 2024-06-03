library(dplyr)
library(stringr)
library(readxl)


get_post <- function(num_samples):

barcode = as.data.frame(read_excel('', sheet=1,col_names = paste('V',rep(1：num_samples),sep = '')))
dim(barcode)
barcode

post = as.data.frame(matrix(nrow=num_samples * num_samples,ncol=3))
l = 1
for(i in 1:num_samples){
    for(j in 1:num_samples){
        post[l,1] = barcode[i,j]
        post[l,2] = j
        post[l,3] = num_samples-i+1
        l = l+1
    }
}
colnames(post) = c('barcode', 'x','y')
	return(post)

data_post <- get_post()

write.csv(post,'')
