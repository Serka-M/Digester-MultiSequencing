Merge-multi
================

### Description: R code for combining outputs of multiple metagenomics sequencing approaches

### Load dependencies

``` r
library(stringr)
library(dplyr)
```

### Load dataframes from multiple sequencing approaches.

#### cov/rel - merged bin coverage/relative abundace profiles for different read datasets by Cover-M

#### bin\_\* and gen\_\* are “results\_bins.tsv”,“results\_general.tsv” files made by “merge-single” script. “Cdb.csv” is the bin clustering file from dRep

``` r
cov_pb = read.table("abund/cov_pbccs.tsv", sep="\t", header=T)
rel_pb = read.table("abund/rel_pbccs.tsv", sep="\t", header=T)
bins_pb = read.table("bins-df/pbccs_bins.tsv", sep="\t", header=T)
gen_pb = read.table("bins-df/pbccs_general.tsv", sep="\t", header=T)

cov_il = read.table("abund/cov_ilmn.tsv", sep="\t", header=T)
rel_il = read.table("abund/rel_ilmn.tsv", sep="\t", header=T)
bins_il = read.table("bins-df/ilmn_bins.tsv", sep="\t", header=T)
gen_il = read.table("bins-df/ilmn_general.tsv", sep="\t", header=T)

cov_r9 = read.table("abund/cov_r9.tsv", sep="\t", header=T)
rel_r9 = read.table("abund/rel_r9.tsv", sep="\t", header=T)
bins_r9 = read.table("bins-df/r9_bins.tsv", sep="\t", header=T)
gen_r9 = read.table("bins-df/r9_general.tsv", sep="\t", header=T)

cov_r103 = read.table("abund/cov_r10.tsv", sep="\t", header=T)
rel_r103 = read.table("abund/rel_r10.tsv", sep="\t", header=T)
bins_r103 = read.table("bins-df/r10_bins.tsv", sep="\t", header=T)
gen_r103 = read.table("bins-df/r10_general.tsv", sep="\t", header=T)

cov_r104 = read.table("abund/cov_r104.tsv", sep="\t", header=T)
rel_r104 = read.table("abund/rel_r104.tsv", sep="\t", header=T)
bins_r104 = read.table("bins-df/r104_bins.tsv", sep="\t", header=T)
gen_r104 = read.table("bins-df/r104_general.tsv", sep="\t", header=T)

cov_r9_il = read.table("abund/cov_r9-il.tsv", sep="\t", header=T)
rel_r9_il = read.table("abund/rel_r9-il.tsv", sep="\t", header=T)
bins_r9_il = read.table("bins-df/r9-il_bins.tsv", sep="\t", header=T)
gen_r9_il = read.table("bins-df/r9-il_general.tsv", sep="\t", header=T)

cov_r103_il = read.table("abund/cov_r10-il.tsv", sep="\t", header=T)
rel_r103_il = read.table("abund/rel_r10-il.tsv", sep="\t", header=T)
bins_r103_il = read.table("bins-df/r10-il_bins.tsv", sep="\t", header=T)
gen_r103_il = read.table("bins-df/r10-il_general.tsv", sep="\t", header=T)

cov_r104_il = read.table("abund/cov_r104-il.tsv", sep="\t", header=T)
rel_r104_il = read.table("abund/rel_r104-il.tsv", sep="\t", header=T)
bins_r104_il = read.table("bins-df/r104-il_bins.tsv", sep="\t", header=T)
gen_r104_il = read.table("bins-df/r104-il_general.tsv", sep="\t", header=T)

drep = read.table("Cdb.csv", sep=",", header=T)
gen <- rbind(gen_il,gen_r9,gen_r9_il,gen_r103,gen_r103_il,gen_pb)
```

### Index bin datasets and remove irrelevant columns

``` r
bins_id <- function(bins,name) {
bins$mags_workflow_mode = name
bins$cov_ilm <- NULL
bins$cov_long <- NULL
bins$Rabund_ilm <- NULL
bins$Rabund_long <- NULL
return(bins)}

bins_pb = bins_id(bins_pb,"PacBio CCS")
bins_il = bins_id(bins_il,"Illumina")

bins_r9 = bins_id(bins_r9,"Nanopore R9.4.1")
bins_r103 = bins_id(bins_r103,"Nanopore R10.3")
bins_r104 = bins_id(bins_r104,"Nanopore R10.4")

bins_r9_il = bins_id(bins_r9_il,"Nanopore R9.4.1 + Illumina")
bins_r103_il = bins_id(bins_r103_il,"Nanopore R10.3 + Illumina")
bins_r104_il = bins_id(bins_r104_il,"Nanopore R10.4 + Illumina")
```

### Process bin abundace data and get a brief summary

#### For the abundace profiles, the column positions are fixed to a specific read dataset.

``` r
abund <- function(bins,cov,rel) {
colnames(cov) <- c("bin","cov_il","cov_pb","cov_r9","cov_r103","cov_r104")
colnames(rel) <- c("bin","r_abund_il","r_abund_pb","r_abund_r9","r_abund_r103","r_abund_r104")
bins = merge(bins,cov,by="bin")
bins = merge(bins,rel,by="bin")
return(bins)}

bins_pb  <- abund(bins_pb,cov_pb,rel_pb)
bins_pb$bin = paste("pbccs",bins_pb$bin,sep="_")
message("\nPB all: ",sum(bins_pb$r_abund_pb), " %, PB HQ: ",
        sum(bins_pb[(bins_pb$MAG_status == "HQ"), ]$r_abund_pb)," %")
```

    ## 
    ## PB all: 82.871278094 %, PB HQ: 47.866132025 %

``` r
bins_il  <- abund(bins_il,cov_il,rel_il )
bins_il$bin = paste("ilmn",bins_il$bin,sep="_")
message("IL all: ",sum(bins_il$r_abund_il), " %, IL HQ: ",
        sum(bins_il[(bins_il$MAG_status == "HQ"), ]$r_abund_il)," %")
```

    ## IL all: 76.07319384 %, IL HQ: 16.14419332 %

``` r
bins_r9  <- abund(bins_r9,cov_r9,rel_r9 )
bins_r9$bin = paste("r9",bins_r9$bin,sep="_")
message("\nR9 all: ",sum(bins_r9$r_abund_r9), " %, R9 HQ: ",
        sum(bins_r9[(bins_r9$MAG_status == "HQ"), ]$r_abund_r9)," %")
```

    ## 
    ## R9 all: 86.03866902 %, R9 HQ: 46.289110831 %

``` r
bins_r9_il  <- abund(bins_r9_il,cov_r9_il,rel_r9_il )
bins_r9_il$bin = paste("r9-il",bins_r9_il$bin,sep="_")
message("R9-IL all: ",sum(bins_r9_il$r_abund_r9), " %, R9-IL HQ: ",
        sum(bins_r9_il[(bins_r9_il$MAG_status == "HQ"), ]$r_abund_r9)," %")
```

    ## R9-IL all: 85.212034859 %, R9-IL HQ: 48.790658012 %

``` r
bins_r103  <- abund(bins_r103,cov_r103,rel_r103 )
bins_r103$bin = paste("r103",bins_r103$bin,sep="_")
message("\nR103 all: ",sum(bins_r103$r_abund_r103), " %, R103 HQ: ",
        sum(bins_r103[(bins_r103$MAG_status == "HQ"), ]$r_abund_r103)," %")
```

    ## 
    ## R103 all: 85.855040325 %, R103 HQ: 38.80907224 %

``` r
bins_r103_il  <- abund(bins_r103_il,cov_r103_il,rel_r103_il )
bins_r103_il$bin = paste("r103-il",bins_r103_il$bin,sep="_")
message("R103-IL all: ",sum(bins_r103_il$r_abund_r103), " %, R103-IL HQ: ",
        sum(bins_r103_il[(bins_r103_il$MAG_status == "HQ"), ]$r_abund_r103)," %")
```

    ## R103-IL all: 86.254879576 %, R103-IL HQ: 40.7046396 %

``` r
bins_r104  <- abund(bins_r104,cov_r104,rel_r104)
bins_r104$bin = paste("r104",bins_r104$bin,sep="_")
message("\nR104 all: ",sum(bins_r104$r_abund_r104), " %, R104 HQ: ",
        sum(bins_r104[(bins_r104$MAG_status == "HQ"), ]$r_abund_r104)," %")
```

    ## 
    ## R104 all: 82.854020354 %, R104 HQ: 38.58349084 %

``` r
bins_r104_il  <- abund(bins_r104_il,cov_r104_il,rel_r104_il)
bins_r104_il$bin = paste("r104-il",bins_r104_il$bin,sep="_")
message("R104-IL all: ",sum(bins_r104_il$r_abund_r104), " %, R104-IL HQ: ",
        sum(bins_r104_il[(bins_r104_il$MAG_status == "HQ"), ]$r_abund_r104)," %")
```

    ## R104-IL all: 83.65310614 %, R104-IL HQ: 39.62454212 %

### Make QUAST command file based on bin clustering data

``` r
drep1 <- drep
drep1$extra <- drep1$genome
drep1=cbind(drep1, data.frame(do.call('rbind',strsplit(as.character(drep1$extra), "_",fixed=TRUE))))

n_occur = as.data.frame(table(drep1$secondary_cluster,drep1$X1))

cluster = unique(drep1$secondary_cluster)

write(paste("",sep=""),file="quast_cmds.txt",append=FALSE)
for(i in cluster){
  drep1_tmp = drep1[(drep1$secondary_cluster == i), ]
  if (nrow(drep1_tmp) > 1 && "pbccs" %in% drep1_tmp$X1) {
    ref=drep1_tmp[(drep1_tmp$X1 == "pbccs"), ]
    ref=as.character(ref$genome)
    query=drep1_tmp[(drep1_tmp$X1 != "pbccs"), ]
    query=paste0(query$genome, collapse=" ")
    output=paste("--output-dir $output/",i, sep="")
    write(paste("quast.py",query,"--no-icarus --no-html --no-plots -t 5 -R",ref,output,sep=" "),file="quast_cmds.txt",append=TRUE)
  }}
```

### After running QUAST commands and aggregating the results to “quast.tsv” file, wrangle the data and fix a possible dataframe shift bug for some rows (always inspect manually just to be sure)

``` r
quast=read.delim("quast.tsv", sep="\t", header=T)

quast1=quast[grep("*part",quast$unalignedcontigs), ]
quast2=subset(quast, !(quast$Assembly %in% quast1$Assembly))

quast1=quast1[,c(1,34:39)]
quast2=quast2[,c(1,32:37)]

colnames(quast1) <- c("bin","unalligned_bp","genome_frac","dup_ratio","Ns_per_100kb","MMs_per_100kb","Indels_per_100kb")
colnames(quast2) <- c("bin","unalligned_bp","genome_frac","dup_ratio","Ns_per_100kb","MMs_per_100kb","Indels_per_100kb")
quast2 <- quast2[!is.na(quast2$Indels_per_100kb), ]

quast <- rbind(quast1, quast2)
```

### Merge the multiple dataframes and save file

``` r
merge <- rbind(bins_il, bins_pb, bins_r9, bins_r103, bins_r104, bins_r9_il, bins_r103_il, bins_r104_il)
drep <- drep[,c(1,2)]
colnames(drep) <- c("bin","cluster")
drep$bin <- gsub(".fa", "", drep$bin)

merge <- merge(merge,drep, by="bin", all=TRUE)
merge <- merge(merge,quast, by="bin", all=TRUE)

write.table(merge,"bins.tsv", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
```

### Based on clustering results and coverage filtering, make lists of bins, which are shared (with some differences) between different sequencing approaches

``` r
bins_list <- merge[(merge$Genome_size_bp > 0), ]
bins_list$bin <- paste(bins_list$bin,".fa",sep="")

bins_list <- bins_list[((bins_list$mags_workflow_mode=="PacBio CCS" & bins_list$cov_pb >= 10) |
                     (bins_list$mags_workflow_mode=="Nanopore R9.4.1" & bins_list$cov_r9 >= 10) |
                     (bins_list$mags_workflow_mode=="Nanopore R10.3" & bins_list$cov_r103 >= 10) |
                       (bins_list$mags_workflow_mode=="Nanopore R10.4" & bins_list$cov_r104 >= 10) | 
            (bins_list$mags_workflow_mode=="Nanopore R9.4.1 + Illumina" & bins_list$cov_r9 >= 10 & bins_list$cov_il >= 5) |
            (bins_list$mags_workflow_mode=="Nanopore R10.3 + Illumina" & bins_list$cov_r103 >= 10 & bins_list$cov_il >= 5) |
              (bins_list$mags_workflow_mode=="Nanopore R10.4 + Illumina" & bins_list$cov_r104 >= 10 & bins_list$cov_il >= 5) |
              (bins_list$mags_workflow_mode=="Illumina" & bins_list$cov_il >= 10)), ]

bins_list <- bins_list[!is.na(bins_list$cluster), ]
bins_list <- bins_list %>% group_by(bins_list$cluster) %>% filter(n() == 8)

list_ilm <- bins_list[(bins_list$mags_workflow_mode == "Illumina"), ]
list_ilm <- list_ilm[,1]

list_pbccs <- bins_list[(bins_list$mags_workflow_mode == "PacBio CCS"), ]
list_pbccs <- list_pbccs[,1]

list_r9 <- bins_list[(bins_list$mags_workflow_mode == "Nanopore R9.4.1"), ]
list_r9 <- list_r9[,1]

list_r103 <- bins_list[(bins_list$mags_workflow_mode == "Nanopore R10.3"), ]
list_r103 <- list_r103[,1]

list_r104 <- bins_list[(bins_list$mags_workflow_mode == "Nanopore R10.4"), ]
list_r104 <- list_r104[,1]

list_r9_ilm <- bins_list[(bins_list$mags_workflow_mode == "Nanopore R9.4.1 + Illumina"), ]
list_r9_ilm <- list_r9_ilm[,1]

list_r103_ilm <- bins_list[(bins_list$mags_workflow_mode == "Nanopore R10.3 + Illumina"), ]
list_r103_ilm <- list_r103_ilm[,1]

list_r104_ilm <- bins_list[(bins_list$mags_workflow_mode == "Nanopore R10.4 + Illumina"), ]
list_r104_ilm <- list_r104_ilm[,1]
```

### Save lists

``` r
write.table(list_ilm,"bin-clust/list_ilm.csv", quote=F,row.names=FALSE,col.names=FALSE)

write.table(list_r9,"bin-clust/list_r9.csv", quote=F,row.names=FALSE,col.names=FALSE)

write.table(list_r103,"bin-clust/list_r103.csv", quote=F,row.names=FALSE,col.names=FALSE)

write.table(list_r104,"bin-clust/list_r104.csv", quote=F,row.names=FALSE,col.names=FALSE)

write.table(list_r9_ilm,"bin-clust/list_r9_ilm.csv", quote=F,row.names=FALSE,col.names=FALSE)

write.table(list_r103_ilm,"bin-clust/list_r103_ilm.csv", quote=F,row.names=FALSE,col.names=FALSE)

write.table(list_r104_ilm,"bin-clust/list_r104_ilm.csv", quote=F,row.names=FALSE,col.names=FALSE)

write.table(list_pbccs,"bin-clust/list_pbccs.csv", quote=F,row.names=FALSE,col.names=FALSE)
```
