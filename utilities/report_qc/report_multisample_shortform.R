#!/usr/bin/env Rscript

#####################################
##### SQANTI3 report generation ######
#####################################



### Author: Lorena de la Fuente, Elizabeth Tseng & Francisco J Pardo-Palacios
### Last Modified: 09/23/2020 by etseng@pacb.com

#********************** Taking arguments from python script

# args <- commandArgs(trailingOnly = TRUE)
# class.file <- args[1]
# junc.file <- args[2]
# utilities.path <- args[4]
# saturation.curves <- args[5]
# report.format <- args[6]

get_script_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    
    # First, check if there's a "--file" argument
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg[1])
    } else {
        # If no "--file", then remove arguments that start with a dash and R itself
        args <- args[!(grepl("^-", args) | grepl("Rscript", args))]
        
        # If there's anything left, the first one should be the script
        if (length(args) > 0) {
            script_path <- args[1]
        } else {
            stop("Couldn't determine the script path.")
        }
    }
    
    # Normalize the path
    script_path <- normalizePath(script_path, mustWork = FALSE)
    
    # If the path is not absolute, prepend the current working directory
    if (!startsWith(script_path, "/")) {
        script_path <- file.path(getwd(), script_path)
    }
    
    return(script_path)
}

args <- commandArgs(trailingOnly = TRUE)
sample.names = strsplit(args[1], split = ",")[[1]]
class.files = strsplit(args[2], split = ",")[[1]]
junc.files = strsplit(args[3], split = ",")[[1]]
utilities.path = dirname(get_script_path())
report.prefix = args[4]
saturation.curves = args[5]
report.format = "pdf"

if (length(args) < 4) {
    stop("Incorrect number of arguments! Script usage is: [classification files] [junction files] [output prefix] [True/False for saturation curves]. Abort!")
}

if (!(saturation.curves %in% c('True', 'False'))) {
    stop("Saturation curve argument needs to be 'True' or 'False'. Abort!")
}

if (!(report.format %in% c('pdf', 'html', 'both'))) {
    stop("Report format needs to be: pdf, html, or both. Abort!")
}

# source(paste(utilities.path, "/report_qc/generatePDFreport.R", sep = "/"))
# source(paste(utilities.path, "/report_qc/generatePDFmultireport.R", sep = "/"))

if (saturation.curves=='True'){
    source(paste(utilities.path, "/saturation/plot_saturation.R", sep = "/"))
    source(paste(utilities.path, "/saturation/LR_saturation.R", sep = "/"))
    source(paste(utilities.path, "/saturation/data_prep_saturation.R", sep = "/"))
}

# report.prefix <- strsplit(class.file, "_classification.txt")[[1]][1];
# output_directory <- dirname(class.file)
# output_name <- basename(report.prefix)
pdf.report.file <- paste0(report.prefix, "_SQANTI3_report.pdf");
# html.report.file <- paste0(output_name, "_SQANTI3_report.html");
# class.file2 <- paste(report.prefix, "_classification_TPM.txt", sep='');


#********************** Packages 

library(ggplot2)
library(ggplotify)
library(scales)
library(reshape)
library(grid)
library(gridExtra)
library(NOISeq)
library(rmarkdown)
library(htmltools)
library(DT)
library(plyr)
library(plotly)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(purrr)
# library(R.utils) # needs to be installed for gz reading with data.table::fread


####### SRTM and SNTM functions

STM_function <- function(x){
    five=FALSE
    three=FALSE
    sj=FALSE
    ev_sj <- !is.na(x["min_cov"]) & as.numeric(x["exons"])>1
    ref_TSS <- FALSE
    ref_TTS <- FALSE
    
    if (!is.na(x["diff_to_gene_TSS"])){
        if (abs(as.numeric(x["diff_to_gene_TSS"]))<=50){
            ref_TSS <- TRUE
        }
    }
    
    if (!is.na(x["diff_to_gene_TTS"])){
        if (abs(as.numeric(x["diff_to_gene_TTS"]))<=50){
            ref_TTS <- TRUE
        }
    }
    
    w_cage <- !is.na(x["within_CAGE_peak"]) & x["within_CAGE_peak"]=="True"
    if ( ref_TSS | w_cage  ){
        five=TRUE
    }
    w_polya <- !is.na(x["within_polyA_site"]) & x["within_polyA_site"]=="True"
    if (ref_TTS | w_polya | !is.na(x["polyA_motif"])){
        three=TRUE
    }
    if (x["structural_category"]=="FSM" | x["structural_category"]=="ISM"){
        sj=TRUE
    }else if ( ev_sj ){
        if (as.numeric(x["min_cov"])>0){
            sj=TRUE
        }
    }
    
    if (five & three & sj){
        return("Fully supported")
    }else{
        return("Not fully supported")
    }
}



# ***********************
# ***********************
# PLOTS INDEX


# PLOT length of isoforms
# p.length.all: length of all isoforms, regardless of category
# p.length.cat: length of isoforms, by category
# p.length.exon: length of isoforms, mono- vs mult-exon
# (optional) p.length.all.sample
# (optional) p.length.exon.sample


# p0: Splicing complexity (X) Isoforms per Gene (Y) Number of Genes
# p1: Distribution of Isoform Classification
# p2: Distribution of Ref Lengths (FSM ISM only)
# p3: Distribution of Ref Exons   (FSM ISM only)
# p4: Distribution of Isoform Lengths, by Classification
# p5: Distribution of Exon Counts, by Classification
# p6: Distribution of Mono- vs Multi-Exons, Novel vs Annotated Genes
# p7: Splicing complexity, Novel vs Annotated Genes
# p8: Transcript Expression (SR) by Structural Category { log2(TPM+1) }
# p9: Transcript Expression (FL) by Structural Category { log2(TPM+1) }
# p10: Gene Expression (SR), Annotated vs Novel genes { log2(Gene_TPM+1) }
# p11: Gene Expression (SR), Annotated vs Novel genes { log2(Gene_TPM+1) }
# p13: Gene Expression level in NNC/FSM containing genes
# p13.c: Gene Expression level in NNC/FSM containing genes

# p.classByLen.a: Structural categories with increasing transcript length, absolute
# p.classByLen.b: Structural categories with increasing transcript length, normalized
#
# p21.a: Distance to polyA site for FSM, absolute
# p21.b: Distance to polyA site for FSM, percentage
# p21.dist3.ISM.a: Distance to polyA site for ISM, absolute
# p21.dist3.ISM.b: Distance to polyA site for ISM, percentage
# p22.a: Distance to start site for FSM, absolute
# p22.b: Distance to start site for FSM, percentage

# p23.a: Splice Junctions by Classification (known canonical, known non-canonical, novel canonical, novel non-canonical)
# p23.b: Splice Junctions by Classification (canonical vs non-canonical)
#
# p28.a: Good Quality Control Attributes Across Structural Categories (annot, SJ, coverage)
# p28.aa: Good Quality Control Attributes Across Structural Categories (polyA, Cage)
# p28.a.SJ: Percentage of  All Canonical Junctions
# p28.a.Cov: Percentage of Splice Junctions With Short Reads Coverage
# p28.a.Cage: Percentage of Cage Support
# p28.a.polyA : Percentage of PolyA Support
# p28.a.annot : Percentage of Annotation Support
# p29.a: Splice Junction, % of RT switching, all junctions
# p29.b: Splice Junction, % of RT switching, unique junctions
#
# p30: intra-priming, by Classification
# p31: intra-priming, Mono- vs Multi-Exons
# p32: intra-priming, Coding vs Non-Coding


# ***********************


########## Generating plots

#*** Global plot parameters

xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
subc.levels=c("alternative_3end",'alternative_3end5end', "alternative_5end","reference_match", "3prime_fragment","internal_fragment", "5prime_fragment","combination_of_known_junctions", "combination_of_known_splicesites", "intron_retention","no_combination_of_known_junctions", "mono-exon_by_intron_retention", "at_least_one_novel_splicesite", "mono-exon", "multi-exon")
subc.labels=c("Alternative 3'end", "Alternative 3'5'end", "Alternative 5'end", "Reference match", "3' fragment", "Internal fragment", "5' fragment", "Comb. of annot. junctions", "Comb. of annot. splice sites", "Intron retention", "Not comb. of annot. junctions", "Mono-exon by intron ret.", "At least 1 annot. don./accept.", "Mono-exon", "Multi-exon")
coding.levels=c("coding", "non_coding")
coding.labels=c("Coding", "Non coding")

myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
subcat.palette = c("Alternative 3'end"='#02314d',
                   "Alternative 3'5'end"='#0e5a87',
                   "Alternative 5'end"='#7ccdfc',
                   'Reference match'='#c4e1f2',
                   "3' fragment"='#c4531d',
                   "Internal fragment"='#e37744',  
                   "5' fragment"='#e0936e', 
                   "Comb. of annot. junctions"='#014d02',
                   "Comb. of annot. splice sites"='#379637',  
                   "Intron retention"='#81eb82', 
                   "Not comb. of annot. junctions"='#6ec091',
                   "Mono-exon by intron ret."='#4aaa72',
                   "At least 1 annot. don./accept."='#32734d',
                   "Mono-exon"='#cec2d2',
                   "Multi-exon"='#876a91')



cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")
sample.palette = brewer.pal(length(sample.names), "Set3")
names(sample.palette) = sample.names


mytheme <- theme_classic(base_family = "Helvetica") +
    theme(axis.line.x = element_line(color="black", size = 0.4),
          axis.line.y = element_line(color="black", size = 0.4)) +
    theme(axis.title.x = element_text(size=13),
          axis.text.x  = element_text(size=12),
          axis.title.y = element_text(size=13),
          axis.text.y  = element_text(vjust=0.5, size=12) ) +
    theme(legend.text = element_text(size = 11), legend.title = element_text(size=12), legend.key.size = unit(0.5, "cm")) +
    theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
    theme(plot.margin = unit(c(2.5,1,1,1), "cm"))


# Function to calculate the grid layout
calculate_layout <- function(N) {
    assign("title_height_global", 1.5, envir = .GlobalEnv)
    if(N <= 3) {
        assign("nrow_global", 1, envir = .GlobalEnv)
        assign("ncol_global", N, envir = .GlobalEnv)
    } else if(N <= 8) {
        assign("nrow_global", 2, envir = .GlobalEnv)
        assign("ncol_global", ceiling(N / 2), envir = .GlobalEnv)
    } else {
        assign("nrow_global", 3, envir = .GlobalEnv)
        assign("ncol_global", ceiling(N / 3), envir = .GlobalEnv)
    }
}

# organize_in_grid <- function(objects_list) {
#   # Use the global nrow and ncol
#   do.call(grid.arrange, c(objects_list, list(nrow = nrow_global, ncol = ncol_global)))
# }

# organize_in_grid_with_title <- function(objects_list) {
#     # heights = c(1, rep(6.5, nrow_global))
#     layout_matrix <- matrix(c(rep(1, ncol_global), 2:(length(objects_list))), ncol = ncol_global, byrow = TRUE)
#     grid.arrange(grobs = objects_list, layout_matrix = layout_matrix, heights = grid::unit(heights, "inches"))
# }

organize_in_grid_with_title <- function(title, objects_list) {
    heights = c(1, rep(6.5, nrow_global))
    layout_matrix <- matrix(c(rep(1, ncol_global), 2:(length(objects_list)+1)), ncol = ncol_global, byrow = TRUE)
    grid.arrange(grobs = c(list(title), objects_list), layout_matrix = layout_matrix, heights = grid::unit(heights, "inches"))
}


# Calculate the layout one time based on the number of inputs
calculate_layout(length(class.files))

#********************** Reading Information

########## Classification information

data.class.list = vector("list", length(class.files))

for (i in seq_along(class.files)) {
    print(paste0(Sys.time(), " starting to read input file ", class.files[i]))
    #class.file = class.files[i]
    # if (tools::file_ext(class.file) == "gz") {
    #     data.class <- read.table(connection <- gzfile(class.file, 'rt'), header=T, as.is=T, sep="\t")
    #     close(connection)
    # } else {
    #     data.class = read.table(class.file, header=T, as.is=T, sep="\t")
    # }
    # data.class = fread(class.files[i], drop = c("min_sample_cov", "sd_cov", "n_indels", "n_indels_junc", "bite", "ratio_exp", "ORF_length", "CDS_length", "CDS_start", "CDS_end", "CDS_genomic_start", "CDS_genomic_end", "ORF_seq", "dist_to_polyA_site", "within_polyA_site"))
    data.class = fread(class.files[i], drop = c("associated_transcript", "ref_exons", "diff_to_TSS", "diff_to_TTS", "RTS_stage", "all_canonical", "min_sample_cov", "sd_cov", "FL", "n_indels", "n_indels_junc", "bite", "iso_exp", "gene_exp", "ratio_exp", "FSM_class", "coding", "ORF_length", "CDS_length", "CDS_start", "CDS_end", "CDS_genomic_start", "CDS_genomic_end", "predicted_NMD", "perc_A_downstream_TTS", "seq_A_downstream_TTS", "dist_to_CAGE_peak", "within_CAGE_peak", "dist_to_polyA_site", "within_polyA_site", "polyA_motif", "polyA_dist", "polyA_motif_found", "ORF_seq", "ratio_TSS"))
    # can drop "coding" and related columns since it's not being output correctly in the table to begin with because --skipORF (old buggy dependency)

    rownames(data.class) <- data.class$isoform
    
    # data.class$sample_name = factor(sample.names[i])
    
    data.class$structural_category = factor(data.class$structural_category,
                                            labels = xaxislabelsF1,
                                            levels = xaxislevelsF1,
                                            ordered=TRUE)
    data.class$subcategory = factor(data.class$subcategory,
                                    labels = subc.labels,
                                    levels = subc.levels,
                                    ordered=TRUE)
    #data.class$coding = factor(data.class$coding,
    #                           labels = coding.labels,
    #                           levels = coding.levels,
    #                           ordered=TRUE)
    #legendLabelF1 <- levels(as.factor(data.class$coding));
    
    
    ## data.class$STM <- apply(data.class,1, STM_function)
    
    # Create a new attribute called "novelGene"
    
    data.class$novelGene <- "Annotated Genes"
    data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
    data.class$novelGene = factor(data.class$novelGene,
                                  levels = c("Novel Genes","Annotated Genes"),
                                  ordered=TRUE)
    
    # Create a new attribute called "exonCat"
    
    data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
    data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
    data.class$exonCat = factor(data.class$exonCat,
                                levels = c("Multi-Exon","Mono-Exon"),
                                ordered=TRUE)
    
    # canonical.labels=c("Canonical", "Non-canonical")
    # data.class$all_canonical = factor(data.class$all_canonical,
    #                                   labels=canonical.labels,
    #                                   levels = c("canonical","non_canonical"),
    #                                   ordered=TRUE)
    
    data.class$lenCat <- as.factor(as.integer(data.class$length %/% 1000))
 
    ################################
    data.class.list[[i]] = data.class
    print(paste0(Sys.time(), " done reading ", class.files[i]))
}

rm(data.class)
gc(reset=TRUE, full=TRUE)



pdf.report.file <- paste0(report.prefix, "_SQANTI3_report.pdf");
pdf(file=pdf.report.file, width = 7.5*ncol_global, height = 6.5*nrow_global + title_height_global)


#### ADD stacked version of this data from plot above in the format of plot below where instead of structural categories we use samples


compute_bins <- function(data_tbl, sample_name, bin_width = 100) {
    breaks <- seq(0 - bin_width/2, 
                  to = ceiling((max(data_tbl$length) - bin_width/2) / bin_width) * bin_width + bin_width/2, 
                  by = bin_width)
    
    bins <- cut(data_tbl$length, 
                breaks=breaks, 
                include.lowest=TRUE, 
                labels=FALSE) * bin_width
    
    binned_data <- data.frame(bin = bins, length = data_tbl$length) %>%
        group_by(bin) %>%
        summarize(count=n()) %>%
        ungroup()
    
    all_bins <- data.frame(bin = breaks[-length(breaks)] + bin_width/2)  # Take mid points of the breaks for bins
    binned_data <- all_bins %>%
        left_join(binned_data, by = "bin") %>%
        mutate(count = ifelse(is.na(count), 0, count))
    
    total_count <- nrow(data_tbl)
    binned_data$normalized_count <- binned_data$count / total_count
    binned_data$sample_name <- sample_name
    
    return(binned_data)
}



print(paste0(Sys.time(), " starting to generate transcript length distribution plots"))
binned_data_list <- map2(data.class.list, sample.names, compute_bins)
combined_binned_data = bind_rows(binned_data_list)
rm(binned_data_list)
p.tmp <- ggplot(combined_binned_data, aes(x=bin, y=count, color=sample_name)) +
    geom_line(size=1) +
    scale_color_manual(values = sample.palette, name="Sample") +
    labs(x="Transcript length", y="Read Count", title="Transcript Lengths Distribution by Sample (w/ supplementary and secondary alignments counted for each instance)") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme +
    theme(legend.position="bottom", legend.title=element_blank())

p.tmp2 <- ggplot(combined_binned_data, aes(x=bin, y=normalized_count, color=sample_name)) +
    geom_line(size=1) +
    scale_color_manual(values = sample.palette, name="Sample") +
    labs(x="Transcript length", y="Normalized Read Count", title="Normalized Transcript Lengths Distribution by Sample (w/ supplementary and secondary alignments counted for each instance)") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme +
    theme(legend.position="bottom", legend.title=element_blank())

grid.arrange(p.tmp, p.tmp2, ncol=2)

x_limit = max(sapply(data.class.list, function(dt) quantile(dt$length, 0.99)))

p.tmp <- p.tmp + 
    coord_cartesian(xlim = c(0, x_limit)) +
    labs(title = paste("Transcript Lengths Distribution by Sample (w/ supplementary and secondary alignments counted for each instance)\n\n(Top 1% at >", round(x_limit), "excluded)"))
p.tmp2 <- p.tmp2 + 
    coord_cartesian(xlim = c(0, x_limit)) +
    labs(title = paste("Normalized Transcript Lengths Distribution by Sample  (w/ supplementary and secondary alignments counted for each instance)\n\n(Top 1% at >", round(x_limit), "excluded)"))

grid.arrange(p.tmp, p.tmp2, ncol=2)
rm(p.tmp)
rm(p.tmp2)
rm(combined_binned_data)
gc(reset=TRUE, full=TRUE)
print(paste0(Sys.time(), " done generating transcript length distribution plots"))





### testing
print(paste0(Sys.time(), " starting to calculate x_limit"))
full_data_tbl = NULL
x_limit2 = 0
for (data_tbl in data.class.list) {
    data_tbl_tmp <- data_tbl %>%
        filter(!grepl("_sec\\d+", isoform))
    
    
    # Remove "_supXX" suffix from "isoform"
    data_tbl_tmp$isoform_processed <- gsub("_sup\\d+", "", data_tbl_tmp$isoform)
    
    # Sum lengths by processed isoform
    data_tbl_tmp <- data_tbl_tmp %>%
        group_by(isoform_processed) %>%
        summarize(length = sum(length)) %>%
        ungroup()
    
    x_limit2 = max(x_limit2, quantile(data_tbl_tmp$length, 0.99))
    # full_data_tbl = rbind(full_data_tbl, data_tbl_tmp)
}
rm(data_tbl_tmp)
gc(reset=TRUE, full=TRUE)
# x_limit2 = quantile(full_data_tbl$length, 0.99)
x_limit = x_limit2
print(paste0(Sys.time(), " done calculating x_limit"))


compute_bins <- function(data_tbl, sample_name, bin_width = 100) {
    # Filter out entries with "_supXX" suffix in "isoform"
    data_tbl <- data_tbl %>%
        filter(!grepl("_sup\\d+|_sec\\d+", isoform))
    
    breaks <- seq(0 - bin_width/2, 
                  to = ceiling((max(data_tbl$length) - bin_width/2) / bin_width) * bin_width + bin_width/2, 
                  by = bin_width)
    
    bins <- cut(data_tbl$length, 
                breaks=breaks, 
                include.lowest=TRUE, 
                labels=FALSE) * bin_width
    
    binned_data <- data.frame(bin = bins, length = data_tbl$length) %>%
        group_by(bin) %>%
        summarize(count=n()) %>%
        ungroup()
    
    all_bins <- data.frame(bin = breaks[-length(breaks)] + bin_width/2)  # Take mid points of the breaks for bins
    binned_data <- all_bins %>%
        left_join(binned_data, by = "bin") %>%
        mutate(count = ifelse(is.na(count), 0, count))
    
    total_count <- nrow(data_tbl)
    binned_data$normalized_count <- binned_data$count / total_count
    binned_data$sample_name <- sample_name
    
    return(binned_data)
}
print(paste0(Sys.time(), " starting to generate transcript length distribution plots part2"))
binned_data_list <- map2(data.class.list, sample.names, compute_bins)
combined_binned_data = bind_rows(binned_data_list)
rm(binned_data_list)
p.tmp <- ggplot(combined_binned_data, aes(x=bin, y=count, color=sample_name)) +
    geom_line(size=1) +
    scale_color_manual(values = sample.palette, name="Sample") +
    labs(x="Transcript length", y="Read Count", title="Transcript Lengths Distribution by Sample (w/o supplementary alignments)") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme +
    theme(legend.position="bottom", legend.title=element_blank())

p.tmp2 <- ggplot(combined_binned_data, aes(x=bin, y=normalized_count, color=sample_name)) +
    geom_line(size=1) +
    scale_color_manual(values = sample.palette, name="Sample") +
    labs(x="Transcript length", y="Normalized Read Count", title="Normalized Transcript Lengths Distribution by Sample (w/o supplementary alignments)") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme +
    theme(legend.position="bottom", legend.title=element_blank())

grid.arrange(p.tmp, p.tmp2, ncol=2)

# x_limit = max(sapply(data.class.list, function(dt) quantile(dt$length, 0.99)))

p.tmp <- p.tmp + 
    coord_cartesian(xlim = c(0, x_limit)) +
    labs(title = paste("Transcript Lengths Distribution by Sample (w/o supplementary alignments)\n\n(Top 1% at >", round(x_limit), "excluded)"))
p.tmp2 <- p.tmp2 + 
    coord_cartesian(xlim = c(0, x_limit)) +
    labs(title = paste("Normalized Transcript Lengths Distribution by Sample (w/o supplementary alignments)\n\n(Top 1% at >", round(x_limit), "excluded)"))

grid.arrange(p.tmp, p.tmp2, ncol=2)
rm(p.tmp)
rm(p.tmp2)
rm(combined_binned_data)
gc(reset=TRUE, full=TRUE)
print(paste0(Sys.time(), " done generating transcript length distribution plots part2"))


compute_bins <- function(data_tbl, sample_name, bin_width = 100) {
    # Filter out entries with "_secXX" suffix in "isoform"
    data_tbl <- data_tbl %>%
        filter(!grepl("_sec\\d+", isoform))
    
    # Remove "_supXX" suffix from "isoform"
    data_tbl$isoform_processed <- gsub("_sup\\d+", "", data_tbl$isoform)
    
    # Sum lengths by processed isoform
    data_tbl <- data_tbl %>%
        group_by(isoform_processed) %>%
        summarize(length = sum(length)) %>%
        ungroup()

    breaks <- seq(0 - bin_width/2, 
                  to = ceiling((max(data_tbl$length) - bin_width/2) / bin_width) * bin_width + bin_width/2, 
                  by = bin_width)
    
    bins <- cut(data_tbl$length, 
                breaks=breaks, 
                include.lowest=TRUE, 
                labels=FALSE) * bin_width
    
    binned_data <- data.frame(bin = bins, length = data_tbl$length) %>%
        group_by(bin) %>%
        summarize(count=n()) %>%
        ungroup()
    
    all_bins <- data.frame(bin = breaks[-length(breaks)] + bin_width/2)  # Take mid points of the breaks for bins
    binned_data <- all_bins %>%
        left_join(binned_data, by = "bin") %>%
        mutate(count = ifelse(is.na(count), 0, count))
    
    total_count <- nrow(data_tbl)
    binned_data$normalized_count <- binned_data$count / total_count
    binned_data$sample_name <- sample_name
    
    return(binned_data)
}

print(paste0(Sys.time(), " starting to generate transcript length distribution plots part3"))

binned_data_list <- map2(data.class.list, sample.names, compute_bins)
combined_binned_data = bind_rows(binned_data_list)
rm(binned_data_list)
p.tmp <- ggplot(combined_binned_data, aes(x=bin, y=count, color=sample_name)) +
    geom_line(size=1) +
    scale_color_manual(values = sample.palette, name="Sample") +
    labs(x="Transcript length", y="Read Count", title="Transcript Lengths Distribution by Sample (w/ supplementary alignments added to primary)") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme +
    theme(legend.position="bottom", legend.title=element_blank())

p.tmp2 <- ggplot(combined_binned_data, aes(x=bin, y=normalized_count, color=sample_name)) +
    geom_line(size=1) +
    scale_color_manual(values = sample.palette, name="Sample") +
    labs(x="Transcript length", y="Normalized Read Count", title="Normalized Transcript Lengths Distribution by Sample (w/ supplementary alignments added to primary)") +
    scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
    mytheme +
    theme(legend.position="bottom", legend.title=element_blank())

grid.arrange(p.tmp, p.tmp2, ncol=2)

# x_limit = max(sapply(data.class.list, function(dt) quantile(dt$length, 0.99)))

p.tmp <- p.tmp + 
    coord_cartesian(xlim = c(0, x_limit)) +
    labs(title = paste("Transcript Lengths Distribution by Sample (w/ supplementary alignments added to primary)\n\n(Top 1% at >", round(x_limit), "excluded)"))
p.tmp2 <- p.tmp2 + 
    coord_cartesian(xlim = c(0, x_limit)) +
    labs(title = paste("Normalized Transcript Lengths Distribution by Sample (w/ supplementary alignments added to primary)\n\n(Top 1% at >", round(x_limit), "excluded)"))

grid.arrange(p.tmp, p.tmp2, ncol=2)
rm(p.tmp)
rm(p.tmp2)
rm(combined_binned_data)
gc(reset=TRUE, full=TRUE)
print(paste0(Sys.time(), " done generating transcript length distribution plots part3"))

for (i in seq_along(data.class.list)) {
    data.class.list[[i]] <- data.class.list[[i]] %>%
        filter(!grepl("_sup\\d+|_sec\\d+", isoform))
}
### end of testing


#####################

compute_bins <- function(data_tbl, sample_name, bin_width = 100) {
    breaks <- seq(0 - bin_width/2, 
                  to = ceiling((max(data_tbl$length) - bin_width/2) / bin_width) * bin_width + bin_width/2, 
                  by = bin_width)
    
    bins <- cut(data_tbl$length, 
                breaks=breaks, 
                include.lowest=TRUE, 
                labels=FALSE) * bin_width
    
    binned_data <- data.frame(bins = bins, length = data_tbl$length, structural_category = data_tbl$structural_category) %>%
        group_by(bins, structural_category) %>%
        summarize(count=n(), .groups = 'drop') %>%
        ungroup()
    
    all_bins <- data.frame(bins = breaks[-length(breaks)] + bin_width/2)  # Take mid points of the breaks for bins
    
    # Combine all possible combinations of bins and structural categories
    complete_data <- expand.grid(bins = all_bins$bins, structural_category = unique(data_tbl$structural_category))
    
    binned_data <- complete_data %>%
        left_join(binned_data, by = c("bins", "structural_category")) %>%
        mutate(count = ifelse(is.na(count), 0, count))
    
    total_count <- nrow(data_tbl)
    binned_data$normalized_count <- binned_data$count / total_count
    binned_data$sample_name <- sample_name
    
    return(binned_data)
}
print(paste0(Sys.time(), " starting to generate normalized transcript length distribution plots"))
results <- bind_rows(map2(data.class.list, sample.names, compute_bins))


y_limit = max(results$normalized_count)
for (subcat in unique(results$structural_category)) {
    p.tmp = ggplot(results[results$structural_category == subcat, ], aes(x=bins, y=normalized_count, color=sample_name)) +
        geom_line(size=1) +
        scale_color_manual(values = sample.palette, name="Sample") +
        labs(x="Transcript length", y="Normalized Count", title=paste0("Normalized Transcript Lengths Distribution by Sample for ", subcat)) +
        scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
        mytheme +
        theme(legend.position="bottom", legend.title=element_blank())
    p.tmp2 = p.tmp + coord_cartesian(ylim = c(0, y_limit))
    p.tmp3 = p.tmp + coord_cartesian(xlim = c(0, x_limit))
    p.tmp4 = p.tmp + coord_cartesian(xlim = c(0, x_limit), ylim = c(0, y_limit))
    grid.arrange(p.tmp, p.tmp2, p.tmp3, p.tmp4, ncol=2)
}
rm(results)
rm(p.tmp)
rm(p.tmp2)
rm(p.tmp3)
rm(p.tmp4)
gc(reset=TRUE, full=TRUE)
print(paste0(Sys.time(), " done generating normalized transcript length distribution plots"))

######################

# df.tmp = bind_rows(mapply(function(dt, sample_name) {
#     summary_dt = dt %>% 
#         group_by(structural_category) %>% 
#         summarize(sample = sample_name, 
#                   n = n(),
#                   coding_prop = sum(coding == "Coding")) %>% 
#         mutate(coding_prop = coding_prop / sum(n) * 100,
#                n = n / sum(n) * 100
#         )
#     return(summary_dt)
# }, data.class.list, sample.names, SIMPLIFY = FALSE))
print(paste0(Sys.time(), " starting to generate structural_category plots"))
df.tmp = bind_rows(mapply(function(dt, sample_name) {
    summary_dt = dt %>% 
        group_by(structural_category) %>% 
        summarize(sample = sample_name, 
                  n = n()) %>%
        mutate(total = sum(n)) %>%
        mutate(n = (n / total) * 100) %>%
        select(-total)
    return(summary_dt)
}, data.class.list, sample.names, SIMPLIFY = FALSE))


p.tmp = ggplot(data = df.tmp, aes(x = structural_category, group = sample)) +
    # Full height as lower alpha
    geom_bar(aes(y = n, fill = structural_category, alpha = 0.3), stat = 'identity', position = "dodge", width = 0.7, color = "black", size = 0.3) +
    # Overlap the coding proportion with full alpha
    # geom_bar(aes(y = coding_prop, fill = structural_category, alpha = 1), stat = 'identity', position = "dodge", width = 0.7, color = "black", size = 0.3) +
    geom_text(aes(y = n, label = sample), position = position_dodge(width = 0.7), angle = 90, hjust = -0.1) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = cat.palette) +
    scale_alpha_identity() +
    ylab("Transcripts, %") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45)) +
    ggtitle(paste0("Isoform Distribution Across Structural Categories and Samples", "\n\n")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x  = element_text(margin = ggplot2::margin(17, 0, 0, 0), size = 12)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1))

p.tmp2 = ggplot(data = df.tmp, aes(x = structural_category, group = sample)) +
    # Full height as lower alpha
    geom_bar(aes(y = n, fill = sample, alpha = 0.3), stat = 'identity', position = "dodge", width = 0.7, color = "black", size = 0.3) +
    # Overlap the coding proportion with full alpha
    # geom_bar(aes(y = coding_prop, fill = sample, alpha = 1), stat = 'identity', position = "dodge", width = 0.7, color = "black", size = 0.3) +
    geom_text(aes(y = n, label = sample), position = position_dodge(width = 0.7), angle = 90, hjust = -0.1) +
    scale_x_discrete(drop = FALSE) +
    scale_alpha_identity() +
    ylab("Transcripts, %") +
    mytheme +
    theme(axis.text.x = element_text(angle = 45)) +
    ggtitle(paste0("Isoform Distribution Across Structural Categories and Samples", "\n\n")) +
    theme(axis.title.x = element_blank()) +
    theme(axis.text.x  = element_text(margin = ggplot2::margin(17, 0, 0, 0), size = 12)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(legend.justification = c(1, 1), legend.position = c(1, 1))

grid.arrange(p.tmp, p.tmp2, ncol=2)
rm(p.tmp)
rm(p.tmp2)
rm(df.tmp)
gc(reset=TRUE, full=TRUE)
print(paste0(Sys.time(), " done generating structural_category plots"))

#########################

print(paste0(Sys.time(), " starting to generate subcategory per category plots"))
for (current_category in c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic", "Antisense", "Fusion", "Intergenic", "Genic\nIntron")) {
    df.tmp = bind_rows( mapply(function(dt, sample_name) {
        subset_dt <- dt[(structural_category == current_category & exons > 1), .(subcategory)]
        summary_dt = subset_dt %>% 
            group_by(subcategory) %>% 
            summarize(sample = sample_name,
                      n = n()) %>%
            mutate(total = sum(n)) %>%
            mutate(n = (n / total) * 100) %>%
            select(-total)
        return(summary_dt)
    }, data.class.list, sample.names, SIMPLIFY = FALSE))

#        df.tmp = bind_rows( mapply(function(dt, sample_name) {
#        subset_dt <- dt[(structural_category == current_category & exons > 1), .(subcategory, coding)]
#        summary_dt = subset_dt %>% 
#            group_by(subcategory) %>% 
#            summarize(sample = sample_name,
#                      n = n(),
#                      coding_prop = sum(coding == "Coding")) %>% 
#            mutate(coding_prop = coding_prop / sum(n) * 100,
#                   n = n / sum(n) * 100)
#        return(summary_dt)
#    }, data.class.list, sample.names, SIMPLIFY = FALSE))
    
    if (dim(df.tmp)[1] > 1) {
        p.tmp <- ggplot(data=df.tmp, aes(x = subcategory, group = sample)) +
            # Full height as lower alpha
            geom_bar(aes(y = n, fill = subcategory, alpha = 0.3), stat = 'identity', position = "dodge", width = 0.7, color = "black", size = 0.3) +
            # Overlap the coding proportion with full alpha
            # geom_bar(aes(y = coding_prop, fill = subcategory, alpha = 1), stat = 'identity', position = "dodge", width = 0.7, color = "black", size = 0.3) +
            geom_text(aes(y=n, label=sample, angle=90), hjust=-0.1, position=position_dodge(width=0.7)) +
            scale_x_discrete(drop=TRUE) +
            scale_alpha_identity() +
            # scale_alpha_manual(values=c(1,0.3), name = "Coding prediction", labels = legendLabelF1)+
            ylab("Transcripts, %") +
            mytheme +
            # geom_blank(aes(y=n), stat = "identity") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_fill_manual(values = subcat.palette) +
            ggtitle(paste0("Isoform Distribution Across ",  gsub("\n", " ", current_category), " and Samples", "\n\n"))+
            scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
            theme(axis.title.x=element_blank())
        
        # maybe need to adjust where the legend is shown
        p.tmp2 <- ggplot(data=df.tmp, aes(x = subcategory, groupe = sample)) +
            geom_bar(aes(y = n, fill = sample, alpha = 0.3), stat="identity", color="black", size=0.3, width=0.7, position="dodge") +
            # geom_bar(aes(y = coding_prop, fill = sample, alpha = 1), stat = 'identity', position = "dodge", width = 0.7, color = "black", size = 0.3) +
            scale_x_discrete(drop=TRUE) +
            scale_alpha_identity() +
            # scale_alpha_manual(values=c(1,0.3), name = "Coding prediction", labels = legendLabelF1)+
            ylab("Transcripts, %") +
            mytheme +
            # geom_blank(aes(y=n), stat = "identity") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            # scale_fill_manual(values = subcat.palette), guide='none') +
            ggtitle(paste0("Isoform Distribution Across ",  gsub("\n", " ", current_category), " and Samples", "\n\n"))+
            scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
            theme(axis.title.x=element_blank())
        grid.arrange(p.tmp, p.tmp2, ncol=2)
        rm(p.tmp)
        rm(p.tmp2)
    }
    rm(df.tmp)
    gc(reset=TRUE, full=TRUE)
}
print(paste0(Sys.time(), " done generating subcategory per category plots"))


dev.off()

