Merge-single-il
================

### Description: alternative “merge-single” script, meant to be used with Illumina-only metagenomics data.

#### Since the main purpose is dataframe compatibility with downstream processes, some features will be empty/erroneous and will have to be re-calculated, if relevant for the analysis

``` r
name=read.table(file="name.txt", sep=",", header=F)

bin_contig=read.table("contig_bin.tsv", sep="\t", header=F)
colnames(bin_contig) <- c("contig","bin")

coverage=read.table("cov.tsv", sep="\t", header=T)
coverage[,3]<-NULL
colnames(coverage) <- c("contig","contig_len_bp","cov_long","cov_long_var","cov_ilm","cov_ilm_var")

assembly_general=read.table("assembly.csv", sep=",", header=T)
ilm_general=read.table("nanoq_reads_ilm.ssf", sep=" ", header=T)
np_general=read.table("nanoq_reads_long.ssf", sep=" ", header=T)

checkm=read.table("checkm.tsv", sep="\t", header=T)
colnames(checkm)[1] <- "bin"
colnames(checkm)[2] <- "checkm_marker_lineage"
checkm=checkm[,c(1,2,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,23)]
colnames(checkm) <- gsub("\\.", "_", colnames(checkm))
colnames(checkm) <- gsub("__", "_", colnames(checkm))
colnames(checkm) <- gsub("bp_", "bp", colnames(checkm))

gtdb=read.table("gtdbtk.tsv", sep="\t", header=T)
gtdb=gtdb[,c(1,2,6,7,8,11,12,13,14,17,19,20)]
colnames(gtdb)[1] <- "bin"
colnames(gtdb)[2] <- "gtdb_classification"
colnames(gtdb)[12] <- "gtdb_warning"

cmseq_output=read.table("cmseq_output.csv", sep=",", header=T)
cmseq_output=cmseq_output[,c(1,2,3,4,5,6)]
colnames(cmseq_output)[1] <- "contig"
colnames(cmseq_output)[2] <- "cmseq_total_covered_bases"


silva=read.table("silva_16s.tsv", sep="\t", header=F)
silva=cbind(silva, data.frame(do.call('rbind', strsplit(as.character(silva$V1),':',fixed=TRUE))))
silva=silva[c("X3","V2","V3")]
colnames(silva) <- c("contig","silva_16s","identity")
silva=merge(silva,bin_contig, by="contig")
silva$bin <- as.factor(silva$bin)
silva <- do.call(rbind, lapply(split(silva,silva$bin), function(x) {return(x[which.max(x$identity),])}))
silva <- silva[c("bin","silva_16s")]
colnames(silva) <- c("bin","silva_taxonomy")
silva <- silva[!duplicated(silva$bin),]


rrna=read.table("rrna_stats.csv", sep=",", header=T)

trna=read.table("trna_stats_summary.csv", sep=",", header=T)

prokka=read.table("prokka.csv", sep=",", header=T)
prokka$bin <- gsub("_", "\\.", prokka$bin)


kaiju=read.delim("kaiju.out.names", sep="\t", header=F)
kaiju=kaiju[(kaiju$V1 == "C"), ]
kaiju=cbind(kaiju, data.frame(do.call('rbind', strsplit(as.character(kaiju$V4),';',fixed=TRUE))))
kaiju=kaiju[c("V2","X1","X2","X3","X4","X5","X6","X7")]
colnames(kaiju) <- c("contig","domain","phylum","class","order","family","genus","species")
kaiju=as.data.frame(kaiju, stringsAsFactors = FALSE)
names <- c("domain","phylum","class","order","family","genus","species")
kaiju[names] <- lapply(kaiju[names], gsub, pattern = "NA", replacement = "Unclassified")


AMR=read.delim("abricate_ncbi.tsv", sep="\t", header=T)
VIR=read.delim("abricate_vfdb.tsv", sep="\t", header=T)
colnames(AMR) <- c("contig","amr_gene","amr_identity","amr_name","amr_resistance")
colnames(VIR) <- c("contig","vir_gene","vir_identity","vir_name","vir_resistance")

abricate_general=data.frame(nrow(AMR),nrow(VIR))
colnames(abricate_general) <- c("AMR_count","VIR_count")

if (nrow(VIR) == 0) {VIR[1,] = list("mock_vir","mock_vir","mock_vir","mock_vir","mock_vir")}
if (nrow(AMR) == 0) {AMR[1,] = list("mock_amr","mock_amr","mock_amr","mock_amr","mock_amr")}

AMR_count <- as.data.frame(table(AMR$contig))
VIR_count <- as.data.frame(table(VIR$contig))

colnames(AMR_count) <- c("contig","AMR_count")
colnames(VIR_count) <- c("contig","VIR_count")


general <- cbind(assembly_general,ilm_general,np_general,abricate_general)

contigs <- coverage
contigs$status_circular <- "N"
contigs$status_repeat <- "N"
contigs$graph_path <- "N"

contigs <- merge(contigs,cmseq_output,by="contig")
contigs <- merge(contigs,kaiju,by="contig", all=TRUE)
contigs <- merge(contigs,AMR_count,by="contig", all=TRUE)
contigs <- merge(contigs,VIR_count,by="contig", all=TRUE)

contigs$domain[is.na(contigs$domain)] <- "Unclassified"
contigs$phylum[is.na(contigs$phylum)] <- "Unclassified"
contigs$class[is.na(contigs$class)] <- "Unclassified"
contigs$order[is.na(contigs$order)] <- "Unclassified"
contigs$family[is.na(contigs$family)] <- "Unclassified"
contigs$genus[is.na(contigs$genus)] <- "Unclassified"
contigs$species[is.na(contigs$species)] <- "Unclassified"
contigs$AMR_count[is.na(contigs$AMR_count)] <- "0"
contigs$VIR_count[is.na(contigs$VIR_count)] <- "0"

contigs <- contigs[contigs$contig != "mock_amr", ]
contigs <- contigs[contigs$contig != "mock_vir", ] 

bins <- merge(prokka,rrna, by="bin")
bins <- merge(bins,trna, by="bin")
bins <- merge(bins,checkm, by="bin")
bins <- merge(bins,gtdb, by="bin", all=TRUE)
bins <- merge(bins,silva, by="bin", all=TRUE)

general$mags_workflow_name <- name$V1
general$mags_workflow_date <- name$V2
general$mags_workflow_mode <- name$V3

contigs$mags_workflow_name <- name$V1
contigs$mags_workflow_date <- name$V2
contigs$mags_workflow_mode <- name$V3

bins$mags_workflow_name <- name$V1
bins$mags_workflow_date <- name$V2
bins$mags_workflow_mode <- name$V3


bins$MAG_status <- ifelse((bins$Completeness >= 90 & bins$Contamination <= 5 & bins$unique_tRNA >= 18 & 
              bins$X5S >= 1 & bins$X16S >= 1 & bins$X23S >= 1),"HQ",
              ifelse(bins$Completeness >= 50 & bins$Contamination <= 10, "MQ",
              ifelse(bins$Completeness <= 50 & bins$Contamination <= 10, "LQ", "Contaminated")))


cov_long=read.table("bin_cov_long.tsv", sep="\t", header=T)
cov_short=read.table("bin_abund_long.tsv", sep="\t", header=T)
abund_long=read.table("bin_cov_short.tsv", sep="\t", header=T)
abund_short=read.table("bin_abund_short.tsv", sep="\t", header=T)

cov_abund <- merge(cov_long,cov_short, by="Genome")
cov_abund <- merge(cov_abund,abund_long, by="Genome")
cov_abund <- merge(cov_abund,abund_short, by="Genome")
colnames(cov_abund)<- c("bin","cov_long","cov_ilm","Rabund_long","Rabund_ilm")
bins = merge(bins,cov_abund, by="bin")

bin_calc=contigs
bin_calc <- merge(bin_calc,bin_contig,by="contig")
bin_aggregate <- aggregate(list(bin_calc$cmseq_total_covered_bases,bin_calc$total_polymorphic_bases), by=list(Category=bin_calc$bin), FUN=sum)
colnames(bin_aggregate)<- c("bin","cmseq_covered_bases","cmseq_polymorphic_bases")
bins <- merge(bins,bin_aggregate,by="bin")
bins$percent_snp_rate <- bins$cmseq_polymorphic_bases/bins$cmseq_covered_bases*100

bins$cMAG_status <- "N"


kaiju_domain <- aggregate(bin_calc$contig_len_bp, by=list(bin_calc$bin,bin_calc$domain), FUN=sum)
colnames(kaiju_domain)<- c("bin","kaiju_domain","contig_len_bp")

kaiju_domain_sort <- aggregate(kaiju_domain$contig_len_bp, by = list(kaiju_domain$bin), max)
colnames(kaiju_domain_sort)<- c("bin","contig_len_bp")

kaiju_domain_sort <- merge(kaiju_domain,kaiju_domain_sort, by=c("contig_len_bp","bin"))
kaiju_domain_sort$contig_len_bp <- NULL
bins <- merge(bins,kaiju_domain_sort, by="bin", all=TRUE)


kaiju_phylum <- aggregate(bin_calc$contig_len_bp, by=list(bin_calc$bin,bin_calc$phylum), FUN=sum)
colnames(kaiju_phylum)<- c("bin","kaiju_phylum","contig_len_bp")

kaiju_phylum_sort <- aggregate(kaiju_phylum$contig_len_bp, by = list(kaiju_phylum$bin), max)
colnames(kaiju_phylum_sort)<- c("bin","contig_len_bp")

kaiju_phylum_sort <- merge(kaiju_phylum,kaiju_phylum_sort, by=c("contig_len_bp","bin"))
kaiju_phylum_sort$contig_len_bp <- NULL
bins <- merge(bins,kaiju_phylum_sort, by="bin", all=TRUE)


contigs <- merge(contigs,bin_contig,by="contig", all=TRUE)
contigs <- contigs[(contigs$contig_len_bp > 0), ]
contigs <- contigs[!is.na(contigs$contig),]
contigs$AMR_count <- as.numeric(contigs$AMR_count)
contigs$VIR_count <- as.numeric(contigs$VIR_count)
contigs <- contigs[!duplicated(contigs$contig),]

AMR_count_bin <- aggregate(contigs$AMR_count, by=list(contigs$bin), FUN=sum)
VIR_count_bin <- aggregate(contigs$VIR_count, by=list(contigs$bin), FUN=sum)
colnames(AMR_count_bin)<- c("bin","AMR_count")
colnames(VIR_count_bin)<- c("bin","VIR_count")
bins <- merge(bins,AMR_count_bin, by="bin", all=TRUE)
bins <- merge(bins,VIR_count_bin, by="bin", all=TRUE)
bins <- bins[!duplicated(bins$bin),]


general$contigs_circ <- nrow(contigs[(contigs$status_circular == "Y"), ])
general$contigs_circ_above_HalfMb <- nrow(contigs[(contigs$status_circular == "Y" & contigs$contig_len_bp >= 1000000), ])

general$contigs_abund_bac <- sum(contigs[(contigs$domain == "Bacteria"), ]$contig_len_bp * contigs[(contigs$domain == "Bacteria"), ]$cov_long) / sum(contigs$contig_len_bp * contigs$cov_long) * 100

general$contigs_abund_arc <- sum(contigs[(contigs$domain == "Archaea"), ]$contig_len_bp * contigs[(contigs$domain == "Archaea"), ]$cov_long) / sum(contigs$contig_len_bp * contigs$cov_long) * 100

general$all_mags <- nrow(bins)
general$circ_mags <- nrow(bins[(bins$cMAG_status == "Y"), ])
general$hq_mags <- nrow(bins[(bins$MAG_status == "HQ"), ])
general$mq_mags <- nrow(bins[(bins$MAG_status == "MQ"), ])
general$lq_mags <- nrow(bins[(bins$MAG_status == "LQ"), ])
general$contaminated_mags <- nrow(bins[(bins$MAG_status == "Contaminated"), ])

general$asm_binned <- sum(bins$Genome_size_bp)/general$assembly_size_bp*100
general$Rabund_binned_long <- sum(bins$Rabund_long)
general$Rabund_binned_ilm <- sum(bins$Rabund_ilm)

general$mag_cov_long_mean <- mean(bins$cov_long)
general$mag_cov_ilm_mean <- mean(bins$cov_ilm)
general$contigs_per_mag_mean <- mean(bins$contigs)


write.table(general,"results_general.tsv", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
write.table(contigs,"results_contigs.tsv", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
write.table(bins,"results_bins.tsv", quote=F,row.names=FALSE, col.names=TRUE, sep = "\t")
```
