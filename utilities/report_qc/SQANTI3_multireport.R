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
source(paste(utilities.path, "/report_qc/generatePDFmultireport.R", sep = "/"))

if (saturation.curves=='True'){
    source(paste(utilities.path, "/report_qc/saturation/plot_saturation.R", sep = "/"))
    source(paste(utilities.path, "/report_qc/saturation/LR_saturation.R", sep = "/"))
    source(paste(utilities.path, "/report_qc/saturation/data_prep_saturation.R", sep = "/"))
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
sample.palette = brewer.pal(length(sample.names), "Set1")
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
isoPerGene.list = vector("list", length(class.files))
data.FSMISM.list = vector("list", length(class.files))
data.NICNNC.list = vector("list", length(class.files))
data.other.list = vector("list", length(class.files))
data.FSM.list = vector("list", length(class.files))
data.ISM.list = vector("list", length(class.files))
data.NNC.list = vector("list", length(class.files))
data.NIC.list = vector("list", length(class.files))
data.GenicGenomic.list = vector("list", length(class.files))
data.Antisense.list = vector("list", length(class.files))
data.Fusion.list = vector("list", length(class.files))
data.Intergenic.list = vector("list", length(class.files))
data.GenicIntron.list = vector("list", length(class.files))
# data.alt3end.list = vector("list", length(class.files))
# data.alt35end.list = vector("list", length(class.files))
# data.alt5end.list = vector("list", length(class.files))
data.refmatch.list = vector("list", length(class.files))
# data.3fragment.list = vector("list", length(class.files))
# data.int_fragment.list = vector("list", length(class.files))
# data.5fragment.list = vector("list", length(class.files))
# data.intron_ret_ISM.list = vector("list", length(class.files))
# data.comb_annot_js_NIC.list = vector("list", length(class.files))
# data.comb_annot_ss_NIC.list = vector("list", length(class.files))
# data.intron_ret_NIC.list = vector("list", length(class.files))
# data.mono_ex_intron_ret_NIC.list = vector("list", length(class.files))
# data.comb_annot_js_NNC.list = vector("list", length(class.files))
# data.comb_annot_ss_NNC.list = vector("list", length(class.files))
# data.intron_ret_NNC.list = vector("list", length(class.files))
# data.mono_ex_intron_ret_NNC.list = vector("list", length(class.files))
# data.one_don_acc.list = vector("list", length(class.files))
subcategories.FSM.list = vector("list", length(class.files))
subcategories.ISM.list = vector("list", length(class.files))
subcategories.NIC.list = vector("list", length(class.files))
subcategories.NNC.list = vector("list", length(class.files))
data.junction.list = vector("list", length(class.files))
uniqJunc.list = vector("list", length(class.files))
uniqJuncRTS.list = vector("list", length(class.files))
FL_multisample_indices.list = vector("list", length(class.files))
FL_multisample_names.list = vector("list", length(class.files))

for (i in seq_along(class.files)) {

    class.file = class.files[i]
    if (tools::file_ext(class.file) == "gz") {
        data.class <- read.table(connection <- gzfile(class.file, 'rt'), header=T, as.is=T, sep="\t")
        close(connection)
    } else {
        data.class = read.table(class.file, header=T, as.is=T, sep="\t")
    }

    rownames(data.class) <- data.class$isoform

    data.class$sample_name = factor(sample.names[i])

    data.class$structural_category = factor(data.class$structural_category,
                                            labels = xaxislabelsF1,
                                            levels = xaxislevelsF1,
                                            ordered=TRUE)
    data.class$subcategory = factor(data.class$subcategory,
                                    labels = subc.labels,
                                    levels = subc.levels,
                                    ordered=TRUE)
    data.class$coding = factor(data.class$coding,
                               labels = coding.labels,
                               levels = coding.levels,
                               ordered=TRUE)
    legendLabelF1 <- levels(as.factor(data.class$coding));


    data.class$STM <- apply(data.class,1, STM_function)

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

    canonical.labels=c("Canonical", "Non-canonical")
    data.class$all_canonical = factor(data.class$all_canonical,
                                      labels=canonical.labels,
                                      levels = c("canonical","non_canonical"),
                                      ordered=TRUE)

    data.class$lenCat <- as.factor(as.integer(data.class$length %/% 1000))

    # ----------------------------------------------------------
    # Make "isoPerGene" which is aggregated information by gene
    #  $associatedGene - either the ref gene name or novelGene_<index>
    #  $novelGene      - either "Novel Genes" or "Annotated Genes"
    #  $FSM_class      - "A", "B", or "C"
    #  $geneExp        - gene expression info
    #  $nIso           - number of isoforms associated with this gene
    #  $nIsoCat        - splicing complexity based on number of isoforms
    # ----------------------------------------------------------
    if (!all(is.na(data.class$gene_exp))){
        isoPerGene = aggregate(data.class$isoform,
                               by = list("associatedGene" = data.class$associated_gene,
                                         "novelGene" = data.class$novelGene,
                                         "FSM_class" = data.class$FSM_class,
                                         "geneExp"=data.class$gene_exp),
                               length)
    } else {
        isoPerGene = aggregate(data.class$isoform,
                               by = list("associatedGene" = data.class$associated_gene,
                                         "novelGene" = data.class$novelGene,
                                         "FSM_class" = data.class$FSM_class),
                               length)
    }
    # assign the last column with the colname "nIso" (number of isoforms)
    colnames(isoPerGene)[ncol(isoPerGene)] <- "nIso"


    isoPerGene$FSM_class2 = factor(isoPerGene$FSM_class,
                                   levels = c("A", "B", "C"),
                                   labels = c("MonoIsoform Gene", "MultiIsoform Genes\nwithout expression\nof a FSM", "MultiIsoform Genes\nexpressing at least\none FSM"),
                                   ordered=TRUE)

    isoPerGene$novelGene = factor(isoPerGene$novelGene,
                                  levels = c("Annotated Genes", "Novel Genes"),
                                  ordered=TRUE)

    max_iso_per_gene <- max(isoPerGene$nIso)
    if (max_iso_per_gene >= 6) {
            isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5,max_iso_per_gene+1), labels = c("1", "2-3", "4-5", ">=6"));
    } else if (max_iso_per_gene >= 5) {
            isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3,5), labels = c("1", "2-3", "4-5"));
    } else if (max_iso_per_gene >= 3) {
            isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1,3), labels = c("1", "2-3"));
    } else {
            isoPerGene$nIsoCat <- cut(isoPerGene$nIso, breaks = c(0,1), labels = c("1"));
    }

    ################################

    data.FSMISM <- subset(data.class, structural_category %in% c('FSM', 'ISM'))
    data.NICNNC <- subset(data.class, structural_category %in% c("NIC", "NNC"))
    data.other <- subset(data.class, structural_category %in% c("Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron"))
    data.FSM <- subset(data.class, (structural_category=="FSM" & exons>1))
    data.ISM <- subset(data.class, (structural_category=="ISM" & exons>1))
    data.NNC <- subset(data.class, (structural_category=="NNC" & exons>1))
    data.NIC <- subset(data.class, (structural_category=="NIC" & exons>1))
    data.GenicGenomic <- subset(data.class, (structural_category=="Genic\nGenomic" & exons>1 ))
    data.Antisense <- subset(data.class, (structural_category=="Antisense" & exons>1))
    data.Fusion <- subset(data.class, (structural_category=="Fusion" & exons>1))
    data.Intergenic <- subset(data.class, (structural_category=="Intergenic" & exons>1))
    data.GenicIntron <- subset(data.class, (structural_category=="Genic\nIntron" & exons>1))

    # subcategories data sets
    #"FSM"
    data.alt3end <- subset(data.FSM, (subcategory=="Alternative 3'end"))
    data.alt35end <- subset(data.FSM, (subcategory=="Alternative 3'5'end"))
    data.alt5end <- subset(data.FSM, (subcategory=="Alternative 5'end"))
    data.refmatch <- subset(data.FSM, (subcategory=="Reference match"))
    #"ISM"
    data.3fragment <- subset(data.ISM, (subcategory=="3' fragment"))
    data.int_fragment <- subset(data.ISM, (subcategory=="Internal fragment"))
    data.5fragment <- subset(data.ISM, (subcategory=="5' fragment"))
    data.intron_ret_ISM <- subset(data.ISM, (subcategory=="Intron retention"))
    #"NIC"
    data.comb_annot_js_NIC <- subset(data.NIC, (subcategory=="Comb. of annot. junctions"))
    data.comb_annot_ss_NIC <- subset(data.NIC, (subcategory=="Comb. of annot. splice sites"))
    data.intron_ret_NIC <- subset(data.NIC, (subcategory=="Intron retention"))
    data.mono_ex_intron_ret_NIC <- subset(data.NIC, (subcategory=="Mono-exon by intron ret."))
    #"NNC"
    data.comb_annot_js_NNC <- subset(data.NNC, (subcategory=="Comb. of annot. junctions"))
    data.comb_annot_ss_NNC <- subset(data.NNC, (subcategory=="Comb. of annot. splice sites"))
    data.intron_ret_NNC <- subset(data.NNC, (subcategory=="Intron retention"))
    data.mono_ex_intron_ret_NNC <- subset(data.NNC, (subcategory=="Mono-exon by intron ret."))
    data.one_don_acc <- subset(data.NNC, (subcategory=="At least 1 annot. don./accept."))

    ########### Junction information
    junc.file = junc.files[i]
    if (tools::file_ext(class.file) == "gz") {
        data.junction <- read.table(connection <- gzfile(junc.file, 'rt'), header=T, as.is=T, sep="\t")
        close(connection)
    } else {
        data.junction = read.table(junc.file, header=T, as.is=T, sep="\t")
    }

    # make a unique identifier using chrom_strand_start_end
    data.junction$junctionLabel = with(data.junction, paste(chrom, strand, genomic_start_coord, genomic_end_coord, sep="_"))

    data.junction$SJ_type <- with(data.junction, paste(junction_category,canonical,"SJ", sep="_"))
    data.junction$SJ_type <- factor(data.junction$SJ_type, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                                           labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T)

    data.junction$structural_category = data.class[data.junction$isoform, "structural_category"]


    # see if there are multple FL samples
    FL_multisample_indices <- which(substring(colnames(data.class), 1,3)=="FL.")

    if (any(grep('^FL\\.', names(data.class)))){
        data.class$FL= rowSums(data.class[,grep('^FL\\.', names(data.class))])
    }

    if (!all(is.na(data.class$FL))){
        total_fl <- sum(data.class$FL, na.rm=T)
        #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
        data.class$FL_TPM <- data.class$FL*(10**6)/total_fl
    }

    if (!all(is.na(data.FSMISM$FL))){
        total_fl <- sum(data.FSMISM$FL, na.rm=T)
        #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
        data.FSMISM$FL_TPM <- data.FSMISM$FL*(10**6)/total_fl
    }

    if (!all(is.na(data.NICNNC$FL))){
        total_fl <- sum(data.NICNNC$FL, na.rm=T)
        #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
        data.NICNNC$FL_TPM <- data.NICNNC$FL*(10**6)/total_fl
    }

    if (!all(is.na(data.other$FL))){
        total_fl <- sum(data.other$FL, na.rm=T)
        #data.class$FL_TPM <- round(data.class$FL*(10**6)/total_fl)
        data.other$FL_TPM <- data.other$FL*(10**6)/total_fl
    }

    if (length(FL_multisample_indices)>0){
        FL_multisample_names <- substring(colnames(data.class)[FL_multisample_indices],4)
        # FL_TPM_multisample_names <- c();

        for (j in 1:length(FL_multisample_indices)) {
            k <- FL_multisample_indices[j]
            name <- paste("FL_TPM", FL_multisample_names[j], sep='.')
            name2 <- paste(name, "_log10", sep='')
            # FL_TPM_multisample_names <- c(FL_TPM_multisample_names, name)
            total_fl <- sum(data.class[k])
            data.class[,name] <- data.class[k]*(10**6)/total_fl
            data.class[,name2] <- log10(data.class[k]*(10**6)/total_fl + 1)
        }
        # data.class$structural_category = factor(data.class$structural_category,
        #                                         labels = xaxislevelsF1,
        #                                         levels = xaxislabelsF1,
        #                                         ordered=TRUE)
        # # write.table(data.class, class.file2, quote=F, sep = "\t", row.names = FALSE);
        # data.class$structural_category = factor(data.class$structural_category,
        #                                         labels = xaxislabelsF1,
        #                                         levels = xaxislevelsF1,
        #                                         ordered=TRUE)
        FL_multisample_names.list[[i]] = FL_multisample_names
    }

    ################################
    data.class.list[[i]] = data.class
    data.FSMISM.list[[i]] = data.FSMISM
    data.NICNNC.list[[i]] = data.NICNNC
    data.other.list[[i]] = data.other
    data.FSM.list[[i]] = data.FSM
    data.ISM.list[[i]] = data.ISM
    data.NNC.list[[i]] = data.NNC
    data.NIC.list[[i]] = data.NIC
    data.GenicGenomic.list[[i]] = data.GenicGenomic
    data.Antisense.list[[i]] = data.Antisense
    data.Fusion.list[[i]] = data.Fusion
    data.Intergenic.list[[i]] = data.Intergenic
    data.GenicIntron.list[[i]] = data.GenicIntron

    # subcategories data sets
    #"FSM"
    # data.alt3end.list[[i]] = data.alt3end
    # data.alt35end.list[[i]] = data.alt35end
    # data.alt5end.list[[i]] = data.alt5end
    data.refmatch.list[[i]] = data.refmatch
    #"ISM"
    # data.3fragment.list[[i]] = data.3fragment
    # data.int_fragment.list[[i]] = data.int_fragment
    # data.5fragment.list[[i]] = data.5fragment
    # data.intron_ret_ISM.list[[i]] = data.intron_ret_ISM
    #"NIC"
    # data.comb_annot_js_NIC.list[[i]] = data.comb_annot_js_NIC
    # data.comb_annot_ss_NIC.list[[i]] = data.comb_annot_ss_NIC
    # data.intron_ret_NIC.list[[i]] = data.intron_ret_NIC
    # data.mono_ex_intron_ret_NIC.list[[i]] = data.mono_ex_intron_ret_NIC
    #"NNC"
    # data.comb_annot_js_NNC.list[[i]] = data.comb_annot_js_NNC
    # data.comb_annot_ss_NNC.list[[i]] = data.comb_annot_ss_NNC
    # data.intron_ret_NNC.list[[i]] = data.intron_ret_NNC
    # data.mono_ex_intron_ret_NNC.list[[i]] = data.mono_ex_intron_ret_NNC
    # data.one_don_acc.list[[i]] = data.one_don_acc

    subcategories.FSM.list[[i]] <- list(data.alt3end, data.alt35end, data.alt5end, data.refmatch)
    subcategories.ISM.list[[i]] <- list(data.3fragment, data.int_fragment, data.5fragment, data.intron_ret_ISM)
    subcategories.NIC.list[[i]] <- list(data.comb_annot_js_NIC, data.comb_annot_ss_NIC, data.intron_ret_NIC, data.mono_ex_intron_ret_NIC)
    subcategories.NNC.list[[i]] <- list(data.comb_annot_js_NNC, data.comb_annot_ss_NNC, data.intron_ret_NNC, data.mono_ex_intron_ret_NNC, data.one_don_acc)

    isoPerGene.list[[i]] = isoPerGene

    data.junction.list[[i]] = data.junction

    uniqJunc.list[[i]] <- unique(data.junction[,c("junctionLabel", "SJ_type", "total_coverage_unique")]);
    uniqJuncRTS.list[[i]] <- unique(data.junction[,c("junctionLabel","SJ_type", "RTS_junction")]);

    FL_multisample_indices.list[[i]] = FL_multisample_indices

}

rm(data.class)
rm(data.FSMISM)
rm(data.NICNNC)
rm(data.other)
rm(data.FSM)
rm(data.ISM)
rm(data.NNC)
rm(data.NIC)
rm(data.GenicGenomic)
rm(data.Antisense)
rm(data.Fusion)
rm(data.Intergenic)
rm(data.GenicIntron)
rm(data.refmatch)
rm(data.alt3end)
rm(data.alt35end)
rm(data.alt5end)
rm(data.3fragment)
rm(data.int_fragment)
rm(data.5fragment)
rm(data.intron_ret_ISM)
rm(data.comb_annot_js_NIC)
rm(data.comb_annot_ss_NIC)
rm(data.intron_ret_NIC)
rm(data.mono_ex_intron_ret_NIC)
rm(data.comb_annot_js_NNC)
rm(data.comb_annot_ss_NNC)
rm(data.intron_ret_NNC)
rm(data.mono_ex_intron_ret_NNC)
rm(data.one_don_acc)
rm(isoPerGene)
rm(data.junction)
gc()




# start the PDF plotting
# pdf(file=pdf.report.file, width = 7.5*ncol_global, height = 6.5*nrow_global)


pdf(file=pdf.report.file, width = 7.5*ncol_global, height = 6.5*nrow_global + title_height_global)


# cover
grid.newpage()
cover <- textGrob("SQANTI3 report", gp=gpar(fontface="italic", fontsize=40, col="orangered"))
grid.draw(cover)


# gt1.list = vector("list", length(class.files)+1)
# gt2.list = vector("list", length(class.files)+1)
# gt3.list = vector("list", length(class.files)+1)
# gt4.list = vector("list", length(class.files)+1)
# gt1.list[[1]] = textGrob("Transcript Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
# gt2.list[[1]] = textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
# gt3.list[[1]] = textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
# gt4.list[[1]] = textGrob("Unique Counts Summary", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
# for (i in seq_along(class.files)) {
#     # TABLE 1: Number of isoforms in each structural category
#     
#     freqCat <- as.data.frame(table(data.class.list[[i]]$structural_category))
#     #freqCat$ranking = order(freqCat$Freq,decreasing = T)
#     table1 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Isoforms, count"))
#     # title1 <- textGrob("Transcript Classification\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
#     title1 <- textGrob(paste0(sample.names[1], "\n"), gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
#     gt1.list[[i+1]] <- gTree(children=gList(table1, title1))
# 
#     # TABLE 2: Number of Novel vs Known Genes
#     freqCat = as.data.frame(table(isoPerGene.list[[i]]$novelGene))
#     table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Genes, count"))
#     #title2 <- textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
#     title2 <- textGrob(sample.names[i], gp=gpar(fontface="italic", fontsize=17), vjust = -4)#
#     gt2.list[[i+1]] <- gTree(children=gList(table2, title2))
#     
#     # TABLE 3: Junction Classification
#     
#     uniq_sj_count <- nrow(uniqJunc.list[[i]])
#     
#     freqCat <- as.data.frame(table(uniqJunc.list[[i]]$SJ_type))
#     freqCat$Var1 <- gsub(" ", "", freqCat$Var1)
#     freqCat$Var1 <- gsub("\n", " ", freqCat$Var1)
#     freqCat$Frac <- round(freqCat$Freq*100 / uniq_sj_count, 2)
#     table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","SJs, count","Percent"))
#     # title2 <- textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
#     title2 <- textGrob(sample.names[i], gp=gpar(fontface="italic", fontsize=17), vjust = -5)
#     gt3.list[[i+1]] <- gTree(children=gList(table2, title2))
#     
#     # TABLE 4: Summary number of Unique Isoforms and Unique Genes
#     nGenes = nrow(isoPerGene.list[[i]])
#     nIso = nrow(data.class.list[[i]])
#     sn = paste(sample.names[i], "\n", "Unique Genes: ", nGenes, "\n", "Unique Isoforms: ", nIso)
#     gt4.list[[i+1]] <- textGrob(sn, gp=gpar(fontface="italic", fontsize=17), vjust = 0)
#     
#     # Plot Table 1 and Table 2
# #    grid.arrange(gt4,gt2,gt1, layout_matrix = cbind(c(1,2),c(1,4)))
# #    grid.arrange(gt3)
# }
# 
# # title1 = textGrob("Transcript Classification\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
# organize_in_grid_with_title(gt1.list)
# # title1 = textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
# organize_in_grid_with_title(gt2.list)
# # title1 = textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
# organize_in_grid_with_title(gt3.list)
# # title1 = textGrob("Unique Counts Summary", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
# organize_in_grid_with_title(gt4.list)
# rm(gt1.list)
# rm(gt2.list)
# rm(gt3.list)
# rm(gt4.list)



gt1.list = vector("list", length(class.files))
gt2.list = vector("list", length(class.files))
gt3.list = vector("list", length(class.files))
gt4.list = vector("list", length(class.files))
# gt1.list[[1]] = textGrob("Transcript Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
# gt2.list[[1]] = textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
# gt3.list[[1]] = textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
# gt4.list[[1]] = textGrob("Unique Counts Summary", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
for (i in seq_along(class.files)) {
    # TABLE 1: Number of isoforms in each structural category
    
    freqCat <- as.data.frame(table(data.class.list[[i]]$structural_category))
    #freqCat$ranking = order(freqCat$Freq,decreasing = T)
    table1 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Isoforms, count"))
    # title1 <- textGrob("Transcript Classification\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
    title1 <- textGrob(paste0(sample.names[1], "\n"), gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
    gt1.list[[i]] <- gTree(children=gList(table1, title1))

    # TABLE 2: Number of Novel vs Known Genes
    freqCat = as.data.frame(table(isoPerGene.list[[i]]$novelGene))
    table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","Genes, count"))
    #title2 <- textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -4)
    title2 <- textGrob(sample.names[i], gp=gpar(fontface="italic", fontsize=17), vjust = -4)#
    gt2.list[[i]] <- gTree(children=gList(table2, title2))
    
    # TABLE 3: Junction Classification
    
    uniq_sj_count <- nrow(uniqJunc.list[[i]])
    
    freqCat <- as.data.frame(table(uniqJunc.list[[i]]$SJ_type))
    freqCat$Var1 <- gsub(" ", "", freqCat$Var1)
    freqCat$Var1 <- gsub("\n", " ", freqCat$Var1)
    freqCat$Frac <- round(freqCat$Freq*100 / uniq_sj_count, 2)
    table2 <- tableGrob(freqCat, rows = NULL, cols = c("Category","SJs, count","Percent"))
    # title2 <- textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17), vjust = -5)
    title2 <- textGrob(sample.names[i], gp=gpar(fontface="italic", fontsize=17), vjust = -5)
    gt3.list[[i]] <- gTree(children=gList(table2, title2))
    
    # TABLE 4: Summary number of Unique Isoforms and Unique Genes
    nGenes = nrow(isoPerGene.list[[i]])
    nIso = nrow(data.class.list[[i]])
    sn = paste(sample.names[i], "\n", "Unique Genes: ", nGenes, "\n", "Unique Isoforms: ", nIso)
    gt4.list[[i]] <- textGrob(sn, gp=gpar(fontface="italic", fontsize=17), vjust = 0)
    
    # Plot Table 1 and Table 2
#    grid.arrange(gt4,gt2,gt1, layout_matrix = cbind(c(1,2),c(1,4)))
#    grid.arrange(gt3)
}
page_title = textGrob("Transcript Classification", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, gt1.list)
page_title = textGrob("Gene Classification", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, gt2.list)
page_title = textGrob("Splice Junction Classification", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, gt3.list)
page_title = textGrob("Unique Counts Summary", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, gt4.list)
rm(gt1.list)
rm(gt2.list)
rm(gt3.list)
rm(gt4.list)



# title slide
s <- textGrob("Gene Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)


p.list = vector("list", length(class.files))
#**** PLOT 0: Distribution of Number of Isoforms per Gene
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(isoPerGene.list[[i]], aes(x=nIsoCat, fill=nIsoCat)) +
                   geom_bar(stat="count", aes(y= (..count..)/sum(..count..)*100), color="black", size=0.3, width=0.7) +
                   guides(fill="none") +
                   scale_y_continuous(labels = function(x) format(x), expand = c(0,0)) +
                   scale_fill_manual(values = myPalette[c(2:5)]) +
                   labs(x ="Isoforms per gene", title=paste0(sample.names[i], "\n\n\n"), y = "Genes, %") +
                   mytheme
}
page_title = textGrob("Number of Isoforms per Gene", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)


# p7: Distribution of Number of Isoforms, separated by Novel vs Annotated Genes
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=isoPerGene.list[[i]], aes(x=novelGene)) +
                   geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=nIsoCat), color="black", size=0.3, width=0.5) +
                   scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                                      labels=c("0","25","50","75","100"), expand = c(0,0)) +
                   scale_x_discrete(drop=FALSE) +
                   scale_fill_manual(name = "Isoforms Per Gene",
                                     values = myPalette[c(2:5)]) +
                   ylab("Genes, %") +
                   xlab("Gene Type") +
                   mytheme +
                   theme(axis.title.x=element_blank()) +
                   theme(legend.position="bottom") +
                   guides(fill = guide_legend(keywidth = 0.9, keyheight = 0.9)) +
                   labs(title=paste0(sample.names[i], "\n\n\n"),
                        subtitle="Known vs Novel Genes\n\n")
}
page_title = textGrob("Number of Isoforms per Gene", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)


##**** PLOT 6: Mono vs Multi-exon distribution for Known vs Novel Genes
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(x=novelGene)) +
                   geom_bar(position="fill",aes(y = (..count..)/sum(..count..), fill=exonCat), color="black", size=0.3, width=0.5) +
                   scale_x_discrete(drop=FALSE) +
                   scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                                      labels=c("0","25","50","75","100"), expand = c(0,0)) +
                   scale_fill_manual(name = "Transcript type",
                                     values = myPalette[c(2:5)]) +
                   ylab("Transcripts, %") +
                   mytheme +
                   theme(axis.title.x=element_blank()) +
                   theme(legend.position="bottom") +
                   ggtitle(paste0(sample.names[i],"\n\n" ))
}
page_title = textGrob("Distribution of Mono- vs Multi-Exon Transcripts", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)



p.list2 = vector("list", length(class.files))
##**** PLOT  absolute and normalized % of different categories with increasing transcript length
for (i in seq_along(class.files)) {
    data.class.byLen <- data.class.list[[i]] %>% dplyr::group_by(lenCat, structural_category) %>% dplyr::summarise(count=dplyr::n() ) %>% mutate(perc=count/sum(count))
    data.class.byLen$structural_category <- factor(data.class.byLen$structural_category, levels=(xaxislabelsF1), order=TRUE)

    p.list[[i]] <- ggplot(data.class.byLen, aes(x=lenCat, y=count, fill=factor(structural_category))) +
                   geom_bar(stat='identity', color="black", size=0.3, width=0.85) +
                   scale_fill_manual(values = cat.palette, guide='none', name="Structural Category") +
                   mytheme+
                   theme(legend.position="right") +
                   guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
                   scale_y_continuous(expand=c(0,0))+
                   theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
                   labs(x="Transcript length, kb", y="Counts", title=sample.names[i])


    p.list2[[i]] <- ggplot(data.class.byLen, aes(x=lenCat, y=perc*100, fill=factor(structural_category))) +
                    geom_bar(stat='identity', color ="black", size=0.3, width=0.85) +
                    scale_fill_manual(values = cat.palette, guide='none', name="Structural Category") +
                    mytheme+
                    theme(legend.position="right")  +
                    guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
                    scale_y_continuous(expand=c(0,0))+
                    theme(axis.text.x = element_text(angle = 60, hjust = 1, size=10))+
                    labs(x="Transcript length, kb", y="%", title=paste0(sample.names[i], "\n\n\n"))
}
page_title = textGrob("Structural Categories by Transcript Length", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
page_title = textGrob("Structural Categories by Transcript Length", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list2)
rm(p.list2)
rm(data.class.byLen)


# PLOT 10: Gene Expression, if expression provided
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$iso_exp))) {
        if (!all(is.na(data.class.list[[i]]$iso_exp))) {
            p.list[[i]] <- ggplot(data=isoPerGene.list[[i]], aes(x=novelGene, y=log2(geneExp+1), fill=novelGene)) +
                           geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                           scale_x_discrete(drop=FALSE) +
                           xlab("Structural classification") +
                           ylab("log2(Gene_TPM+1)") +
                           scale_fill_manual(values = myPalette[c(3:4)]) +
                           guides(fill="none") +
                           mytheme +
                           theme(axis.title.x=element_blank()) +
                           ggtitle(paste0(sample.names[i], "\n\n"))
        }
    }
}
page_title = textGrob("Annotated vs Novel Gene Expression", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)


# PLOT 11: Gene FL number, if FL count provided
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$FL))){
        if (!all(is.na(data.class.list[[i]]$FL))){
            FL_gene <- aggregate(as.integer(data.class.list[[i]]$FL), by = list("associatedGene" = data.class.list[[i]]$associated_gene), sum)
            colnames(FL_gene)[ncol(FL_gene)] <- "FL_gene"
            isoPerGene <- merge(isoPerGene.list[[i]], FL_gene, by="associatedGene")
            total_fl <- sum(data.class.list[[i]]$FL, na.rm=T)
            isoPerGene$FL_gene_TPM <- isoPerGene.list[[i]]$FL_gene*(10**6)/total_fl
            
            p.list[[i]] <- ggplot(data=isoPerGene.list[[i]], aes(x=novelGene, y=log2(FL_gene_TPM+1), fill=novelGene)) +
                           geom_boxplot(color="black", size=0.3,outlier.size = 0.2) +
                           scale_x_discrete(drop=FALSE) +
                           ylab("log2(FL_TPM+1)") +
                           scale_fill_manual(values = myPalette[c(3:4)]) +
                           guides(fill="none") +
                           mytheme +
                           theme(axis.title.x=element_blank()) +
                           ggtitle(paste0(sample.names[i], "\n\n"))
            
        }
    }
}
page_title = textGrob("Number of FL reads per Gene by Type of Gene Annotation", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)



########################################
######### LENGTH PLOTS  ################
########################################

# PLOT length of isoforms
# p.length.all: length of all isoforms, regardless of category
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data.class.list[[i]], aes(x=length)) +
                   geom_histogram(binwidth=100) +
                   labs(x="Transcript length", y="Count", title=sample.names[i]) +
                   theme(legend.position="top") +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                   mytheme
}
page_title = textGrob("All Transcript Lengths Distribution", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)


#### ADD stacked version of this data from plot above in the format of plot below where instead of structural categories we use samples


# tmp.data.class = do.call(rbind, data.class.list)
concat.data.class = bind_rows(data.class.list)
# concat.data.class$df = as.factor(concat.data.class$df)
p.tmp <- ggplot(concat.data.class, aes(x=length, color=sample_name)) +
                geom_freqpoly(binwidth=100, size=1) +
                scale_color_manual(values = sample.palette, name="Sample") +
                labs(x="Transcript length", y="Count", title="Transcript Lengths Distribution by Sample") +
                scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                mytheme+
                theme(legend.position="bottom", legend.title=element_blank())
#grid.arrange(p.tmp)

p.tmp2 <- ggplot(concat.data.class, aes(x=length, color=sample_name, y=..density..)) +
                 geom_density(adjust = 1/5, size=1) +  # the adjust argument controls the bandwidth of the density estimate
                 scale_color_manual(values = sample.palette, name="Sample") +
                 labs(x="Transcript length", y="Density", title="Normalized Transcript Lengths Distribution by Sample") +
                 scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                 mytheme +
                 theme(legend.position="bottom", legend.title=element_blank())
grid.arrange(p.tmp, p.tmp2, ncol=2)

x_limit <- quantile(concat.data.class$length, 0.99)
p.tmp <- p.tmp + 
        coord_cartesian(xlim = c(0, x_limit)) +
        labs(title = paste("Transcript Lengths Distribution by Sample (Top 1% at >", round(x_limit), "excluded)"))
p.tmp2 <- p.tmp2 + 
         coord_cartesian(xlim = c(0, x_limit)) +
         labs(title = paste("Normalized Transcript Lengths Distribution by Sample (Top 1% at >", round(x_limit), "excluded)"))
grid.arrange(p.tmp, p.tmp2, ncol=2)
rm(p.tmp)
rm(p.tmp2)



# p.length.cat: length of isoforms, by category
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data.class.list[[i]], aes(x=length, color=structural_category)) +
                   geom_freqpoly(binwidth=100, size=1) +
                   scale_color_manual(values = cat.palette, name="Structural Category") +
                   labs(x="Transcript length", y="Count", title=sample.names[i]) +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                   mytheme+
                   theme(legend.position="bottom", legend.title=element_blank())
}
page_title = textGrob("Transcript Lengths Distribution by Structural Category", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

# p.length.exon: length of isoforms, mono- vs mult-exon/ufrc/conesa/fpardopalacios/SQANTI_QDE/SQANTI3/melanoma_example/melanoma_chr13_tappAS_annot_from_SQANTI3.gff3
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data.class.list[[i]], aes(x=length, color=exonCat)) +
                   geom_freqpoly(binwidth=100, size=1) +
                   labs(x="Transcript length", y="Count", title=sample.names[i]) +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                   mytheme +
                   theme(legend.position="bottom", legend.title=element_blank()) 
}
page_title = textGrob("Mono- vs Multi- Exon Transcript Lengths Distribution", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)



# (optional) p.length.all.sample: length of all isoforms by sample
if (length(FL_multisample_indices.list[[1]])>0) {  # has multiple samples
    p.list2 = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        df.length_by_sample <- data.frame();
        for (j in seq_along(FL_multisample_indices.list[[i]])) {
            k <- FL_multisample_indices.list[[i]][j];
            df_new <- data.class.list[[i]][data.class[k]>0,c("isoform", "length", "exonCat")];
            df_new$sample <- FL_multisample_names.list[[i]][j];
            df.length_by_sample <- rbind(df.length_by_sample, df_new);
        }

        p.list[[i]] <- ggplot(df.length_by_sample, aes(x=length, color=sample)) +
                       geom_freqpoly(binwidth=100) +
                       scale_color_manual(values = cat.palette, name="Structural Category") +
                       labs(x="Transcript length", y="Count", title=sample.names[i]) +
                       theme(legend.position="bottom") +
                       scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                       mytheme

        p.list2[[i]] <- ggplot(df.length_by_sample, aes(x=length, color=sample, lty=exonCat)) +
                        geom_freqpoly(binwidth=100) +
                        labs(x="Transcript length", y="Count", title=sample.names[i]) +
                        theme(legend.position="bottom") +
                        scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                        mytheme
    } ### else add a filler maybe? to avoid shift between pages
    page_title = textGrob("Transcript Lengths by Sample", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    page_title = textGrob("Mono- vs Multi-Exons Transcript Lengths by Sample", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list2)
    rm(p.list2)
}



# 2. general parameters by structual categories
s <- textGrob("Structural Isoform Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)



#**** PLOT 1: Structural Classification
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(x=structural_category)) +
                   geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=coding, fill=structural_category), color="black", size=0.3, width=0.7) +
                   #geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
                   scale_x_discrete(drop=FALSE) +
                   scale_alpha_manual(values=c(1,0.3),
                                      name = "Coding prediction",
                                      labels = legendLabelF1)+
                   xlab("") +
                   ylab("Transcripts, %") +
                   mytheme +
                   geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
                   theme(axis.text.x = element_text(angle = 45)) +
                   scale_fill_manual(values = cat.palette, guide='none') +
                   ggtitle(paste0(sample.names[i], "\n\n" )) +
                   theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                   theme(legend.justification=c(1,1), legend.position=c(1,1))
}
page_title = textGrob("Isoform Distribution Across Structural Categories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)



p1.s.titles = list("Isoform Distribution Across FSM\n\n",
                   "Isoform Distribution Across ISM\n\n",
                   "Isoform Distribution Across NNC\n\n",
                   "Isoform Distribution Across NIC\n\n",
                   "Isoform Distribution Across Genic Genomic\n\n",
                   "Isoform Distribution Across Antisense\n\n",
                   "Isoform Distribution Across Fusion\n\n",
                   "Isoform Distribution Across Intergenic\n\n",
                   "Isoform Distribution Across Genic Intron\n\n")

categories.list=list(data.FSM.list, data.ISM.list, data.NNC.list, data.NIC.list, data.GenicGenomic.list, data.Antisense.list, data.Fusion.list, data.Intergenic.list, data.GenicIntron.list)

for(j in 1:length(categories.list)) {
    c <- categories.list[[j]]
    for (i in seq_along(class.files)) {
        if (!(dim(c[[i]])[1]==0)) {
            p.list[[i]] <- ggplot(data=c[[i]], aes(x=subcategory)) +
                           geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=coding, fill=subcategory), color="black", size=0.3, width=0.7) +
                           scale_x_discrete(drop=TRUE) +
                           scale_alpha_manual(values=c(1,0.3), name = "Coding prediction", labels = legendLabelF1)+
                           ylab("Transcripts, %") +
                           mytheme +
                           geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
                           theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                           scale_fill_manual(values = subcat.palette, guide='none') +
                           ggtitle(paste0(sample.names[i], "\n\n"))+
                           scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                           theme(axis.title.x=element_blank())  
        }
    }
    page_title = textGrob(p1.s.titles[j] , gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(categories.list)


#****  PLOT 4: Transcript lengths by category
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(x=structural_category, y=length, fill=structural_category)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   scale_x_discrete(drop=FALSE) +
                   ylab("Transcript Length (bp)") +
                   scale_fill_manual(values = cat.palette) +
                   guides(fill="none") +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   ggtitle(paste0(sample.names[i], "\n\n")) +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Transcript Lengths by Structural Classification", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=structural_category, y=length, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   scale_x_discrete(drop=TRUE) +
                   ylab("Transcript Length (bp)") +
                   scale_fill_manual(values = subcat.palette, drop=TRUE) +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(legend.position="right", legend.title=element_blank()) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   ggtitle(paste0(sample.names[i], "\n\n")) +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Transcript Lengths by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.NICNNC.list[[i]], aes(x=structural_category, y=length, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   scale_x_discrete(drop=TRUE) +
                   ylab("Transcript Length (bp)") +
                   scale_fill_manual(values = subcat.palette, drop=TRUE) +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(legend.position="right", legend.title=element_blank())+
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   ggtitle(paste0(sample.names[i], "\n\n")) +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Transcript Lengths by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.other.list[[i]], aes(x=structural_category, y=length, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   scale_x_discrete(drop=TRUE) +
                   ylab("Transcript Length (bp)") +
                   scale_fill_manual(values = subcat.palette, drop=TRUE) +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(legend.position="right", legend.title=element_blank())+
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   ggtitle(paste0(sample.names[i], "\n\n")) +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Transcript Lengths by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)




##**** PLOT 5: Exon counts by category
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(x=structural_category, y=exons, fill=structural_category)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   ylab("Number of exons") +
                   scale_x_discrete(drop=FALSE) +
                   scale_fill_manual(values = cat.palette) +
                   guides(fill="none") +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   ggtitle(paste0(sample.names[i], "\n\n")) +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Exon Counts by Structural Classification", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

###Exon counts by subcategory
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=structural_category, y=exons, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   theme(legend.position="right", legend.title=element_blank()) +
                   theme(axis.title.x=element_blank())+
                   ylab("Number of exons") +
                   scale_x_discrete(drop=TRUE) +
                   scale_fill_manual(values = subcat.palette) +
                   ggtitle(paste0(sample.names[i], "\n\n"))
}
page_title = textGrob("Exon Counts by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.NICNNC.list[[i]], aes(x=structural_category, y=exons, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   ylab("Number of exons") +
                   scale_x_discrete(drop=TRUE) +
                   scale_fill_manual(values = subcat.palette) +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   theme(legend.position="right", legend.title=element_blank())+
                   ggtitle(paste0(sample.names[i]), "\n\n") +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Exon Counts by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.other.list[[i]], aes(x=structural_category, y=exons, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   ylab("Number of exons") +
                   scale_x_discrete(drop=TRUE) +
                   scale_fill_manual(values = subcat.palette) +
                   mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   theme(legend.position="right", legend.title=element_blank())+
                   ggtitle(paste0(sample.names[i]), "\n\n") +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Exon Counts by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)


##### STM plots
data.FSMISMNICNNC.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    data.FSMISMNICNNC.list[[i]] <- rbind(data.FSMISM.list[[i]], data.NICNNC.list[[i]])
    p.list[[i]] <- ggplot(data=data.FSMISMNICNNC.list[[i]], aes(x=structural_category)) +
                   geom_bar(aes(y = (..count..), alpha=STM, fill=structural_category), color="black", size=0.3, width=0.7) +
                   scale_x_discrete(drop=TRUE) +
                   scale_alpha_manual(values=c(1,0.3),
                                      name = "Supported Transcript Model",
                                      guide = "legend")+
                   xlab("") +
                   ylab("Transcripts, count") +
                   mytheme +
                   geom_blank(aes(y=(..count..)), stat = "count") +
                   theme(axis.text.x = element_text(angle = 45)) +
                   scale_fill_manual(values = cat.palette, guide='none') +
                   ggtitle(paste0(sample.names[i], "\n\n")) +
                   theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                   theme(legend.position = "right")
}
page_title = textGrob("Isoform Distribution Across Structural Categories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.FSMISMNICNNC.list[[i]], aes(x=structural_category)) +
                   geom_bar(aes(y = (..count..), alpha=STM, fill=structural_category), position="fill", color="black", size=0.3, width=0.7) +
                   scale_x_discrete(drop=TRUE) +
                   scale_alpha_manual(values=c(1,0.3),
                                      name = "Supported Transcript Model",
                                      guide = "legend")+
                   xlab("") +
                   ylab("Transcripts, %") +
                   mytheme +
                   theme(axis.text.x = element_text(angle = 45)) +
                   scale_fill_manual(values = cat.palette, guide='none') +
                   ggtitle(paste0(sample.names[i], "\n\n")) +
                   theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1)), labels = scales::percent) +
                   theme(legend.position = "right")
}
page_title = textGrob("Isoform Distribution Across Structural Categories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=subcategory)) +
                   geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), color="black", size=0.3, width=0.7) +
                   scale_x_discrete(drop=TRUE) +
                   scale_alpha_manual(values=c(1,0.3),
                                      name = "Supported Transcript Model",
                                      guide = "legend")+
                   xlab("") +
                   ylab("Transcripts, count") +
                   mytheme +
                   facet_grid(.~ structural_category, scales = "free_x") +
                   theme(axis.text.x = element_text(angle = 90)) +
                   scale_fill_manual(values = subcat.palette, guide = "none") +
                   ggtitle(paste0(sample.names[i], "\n\n"),
                           subtitle = "FSM and ISM" ) +
                   theme(axis.title.x=element_blank()) +  theme(axis.text.x=element_text(size=10)) +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                   theme(legend.position = "right")
}
page_title = textGrob("Isoform Distribution Across Structural Subcategories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=subcategory)) +
                   geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), position="fill", color="black", size=0.3, width=0.7) +
                   scale_x_discrete(drop=TRUE) +
                   scale_alpha_manual(values=c(1,0.3),
                                      name = "Supported Transcript Model",
                                      guide = "legend")+
                   xlab("") +
                   ylab("Transcripts, %") +
                   mytheme +
                   facet_grid(.~ structural_category, scales = "free_x") +
                   theme(axis.text.x = element_text(angle = 90)) +
                   scale_fill_manual(values = subcat.palette, guide = "none") +
                   ggtitle(paste0(sample.names[i], "\n\n"),
                           subtitle = "FSM and ISM" ) +
                   theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
                   scale_y_continuous(expand=expansion(mult = c(0,0)), labels = scales::percent) +
                   theme(legend.position = "right")
}
page_title = textGrob("Isoform Distribution Across Structural Subcategories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.NICNNC.list[[i]], aes(x=subcategory)) +
                   geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), color="black", size=0.3, width=0.7) +
                   scale_x_discrete(drop=TRUE) +
                   scale_alpha_manual(values=c(1,0.3),
                                      name = "Supported Transcript Model",
                                      guide = "legend")+
                   xlab("") +
                   ylab("Transcripts, count") +
                   mytheme +
                   facet_grid(.~ structural_category, scales = "free_x") +
                   theme(axis.text.x = element_text(angle = 90)) +
                   scale_fill_manual(values = subcat.palette, guide='none') +
                   ggtitle(paste0(sample.names[i], "\n\n"),
                                   subtitle = "NIC and NNC") +
                   theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
                   scale_y_continuous(expand=expansion(mult = c(0,0.1))) +
                   theme(legend.position = "right")
}
page_title = textGrob("Isoform Distribution Across Structural Subcategories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.NICNNC.list[[i]], aes(x=subcategory)) +
                   geom_bar(aes(y = (..count..), alpha=STM, fill=subcategory), position="fill", color="black", size=0.3, width=0.7) +
                   scale_x_discrete(drop=TRUE) +
                   scale_alpha_manual(values=c(1,0.3),
                                      name = "Supported Transcript Model",
                                      guide = "legend")+
                   xlab("") +
                   ylab("Transcripts, %") +
                   mytheme +
                   facet_grid(.~ structural_category, scales = "free_x") +
                   theme(axis.text.x = element_text(angle = 90)) +
                   scale_fill_manual(values = subcat.palette, guide='none') +
                   ggtitle(paste0(sample.names[i], "\n\n"),
                                   subtitle = "NIC and NNC" ) +
                   theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(size=10)) +
                   scale_y_continuous(expand=expansion(mult = c(0,0)), labels = scales::percent) +
                   theme(legend.position = "right")
}
page_title = textGrob("Isoform Distribution Across Structural Subcategories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)




# requires short reads
##**** PLOT 8: Expression, if isoform expression provided (iso_exp is in TPM)
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$iso_exp))) {
        p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(x=structural_category, y=log2(iso_exp+1), fill=structural_category)) +
                       geom_boxplot(color="black", size=0.3,  outlier.size = 0.2) +
                       scale_x_discrete(drop=FALSE) +
                       ylab("log2(TPM+1)") +
                       scale_fill_manual(values = cat.palette) +
                       guides(fill="none") +
                       mytheme  + theme(axis.text.x = element_text(angle = 45)) +
                       theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                       theme(axis.title.x=element_blank()) +
                       ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Transcript Expression by Structural Category", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    if (!all(is.na(data.FSMISM.list[[i]]$iso_exp))) {
    ###Expression, if isoform expression provided (iso_exp is in TPM) by subcategory
        p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
                      geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
                      scale_x_discrete(drop=TRUE) +
                      facet_grid(.~ structural_category, scales = "free_x") +
                      ylab("log2(TPM+1)") +
                      scale_fill_manual(values = subcat.palette, guide="none") +
                      mytheme  + theme(axis.text.x = element_text(angle = 90)) +
                      theme(legend.position="right", legend.title=element_blank()) +
                      theme(axis.text.x  = element_text(size=10))+
                      theme(axis.title.x=element_blank()) +
                      ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Transcript Expression by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    if (!all(is.na(data.NICNNC.list[[i]]$iso_exp))) {
        p.list[[i]] <- ggplot(data=data.NICNNC.list[[i]], aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
                       geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
                       scale_x_discrete(drop=TRUE) +
                       facet_grid(.~ structural_category, scales = "free_x") +
                       ylab("log2(TPM+1)") +
                       scale_fill_manual(values = subcat.palette, guide="none") +
                       mytheme  + theme(axis.text.x = element_text(angle = 90)) +
                       theme(legend.position="right", legend.title=element_blank()) +
                       theme(axis.text.x  = element_text(size=10))+
                       theme(axis.title.x=element_blank()) +
                       ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Transcript Expression by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    if (!all(is.na(data.other.list[[i]]$iso_exp))) {
        p.list[[i]] <- ggplot(data=data.other.list[[i]], aes(x=subcategory, y=log2(iso_exp+1), fill=subcategory)) +
                       geom_boxplot(color="black", size=0.3,  outlier.size = 0.2, position="dodge") +
                       scale_x_discrete(drop=TRUE) +
                       facet_grid(.~ structural_category, scales = "free_x") +
                       ylab("log2(TPM+1)") +
                       scale_fill_manual(values = subcat.palette, guide="none") +
                       mytheme  + theme(axis.text.x = element_text(angle = 90)) +
                       theme(legend.position="right", legend.title=element_blank()) +
                       theme(axis.text.x  = element_text(size=10))+
                       theme(axis.title.x=element_blank()) +
                       ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Transcript Expression by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

# requires --fl_count
# PLOT 9: FL number, if FL count provided
# convert FL count to TPM
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$FL))) {
        p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(x=structural_category, y=log2(FL_TPM+1), fill=structural_category)) +
                       geom_boxplot(color="black", size=0.3, outlier.size=0.2) +
                       ylab("log2(FL_TPM+1)") +
                       scale_x_discrete(drop=FALSE) +
                       scale_fill_manual(values = cat.palette) +
                       guides(fill="none") +
                       mytheme +
                       theme(axis.text.x = element_text(angle = 45)) +
                       theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                       theme(axis.title.x=element_blank()) +
                       ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Long Reads Count by Structural Category", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    if (!all(is.na(data.FSMISM.list[[i]]$FL))) {
        p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
                       geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
                       facet_grid(.~ structural_category, scales = "free_x") +
                       ylab("log2(FL_TPM+1)") +
                       scale_x_discrete(drop=TRUE) +
                       scale_fill_manual(values = subcat.palette, guide="none") +
                       mytheme +
                       theme(legend.position="right", legend.title=element_blank()) +
                       theme(axis.text.x = element_text(angle = 90)) +
                       theme(axis.text.x  = element_text(size=10))+
                       theme(axis.title.x=element_blank()) +
                       ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Long Reads Count by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    if (!all(is.na(data.NICNNC.list[[i]]$FL))) {
        p.list[[i]] <- ggplot(data=data.NICNNC.list[[i]], aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
                       geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
                       facet_grid(.~ structural_category, scales = "free_x") +
                       ylab("log2(FL_TPM+1)") +
                       scale_x_discrete(drop=TRUE) +
                       scale_fill_manual(values = subcat.palette, guide="none") +
                       mytheme +
                       theme(legend.position="right", legend.title=element_blank()) +
                       theme(axis.text.x = element_text(angle = 90)) +
                       theme(axis.text.x  = element_text(size=10))+
                       theme(axis.title.x=element_blank()) +
                       ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Long Reads Count by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    if (!all(is.na(data.other.list[[i]]$FL))) {
        p.list[[i]] <- ggplot(data=data.other.list[[i]], aes(x=subcategory, y=log2(FL_TPM+1), fill=subcategory)) +
                       geom_boxplot(color="black", size=0.3, outlier.size=0.1) +
                       facet_grid(.~ structural_category, scales = "free_x") +
                       ylab("log2(FL_TPM+1)") +
                       scale_x_discrete(drop=TRUE) +
                       scale_fill_manual(values = subcat.palette, guide="none") +
                       mytheme +
                       theme(legend.position="right", legend.title=element_blank()) +
                       theme(axis.text.x = element_text(angle = 90)) +
                       theme(axis.text.x  = element_text(size=10))+
                       theme(axis.title.x=element_blank()) +
                       ggtitle(paste0(sample.names[i], "\n\n"))
    }
}
page_title = textGrob("Long Reads Count by Subcategory", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)






#**** PLOTS 2-3: refLength and refExons for ISM and FSM transcripts. Plot if any ISM or FSM transcript

for (i in seq_along(class.files)) {
    if (nrow(data.FSMISM.list[[i]]) > 0) {
        p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=structural_category, y=ref_length/1000, fill=structural_category)) +
                       geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
                       scale_fill_manual(values = myPalette) +
                       scale_x_discrete(drop = TRUE) +
                       guides(fill="none") +
                       xlab("") +
                       ylab("Matched reference length, kb") +
                       labs(title=paste0(sample.names[i], "\n\n\n"),
                            subtitle="Applicable Only to FSM and ISM Categories\n\n")
    }
}
page_title = textGrob("Length Distribution of Matched Reference Transcripts", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)

for (i in seq_along(class.files)) {
    if (nrow(data.FSMISM.list[[i]]) > 0) {
        p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(x=structural_category, y=ref_exons, fill=structural_category)) +
                       geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                       scale_x_discrete(drop = TRUE) +
                       xlab("") +
                       ylab("Matched reference exon count") +
                       scale_fill_manual(values = myPalette) +
                       guides(fill="none") +
                       mytheme +
                       labs(title=paste0(sample.names[i], "\n\n\n"),
                            subtitle="Applicable Only to FSM and ISM Categories\n\n")
    }
}
page_title = textGrob("Exon Count Distribution of Matched Reference Transcripts", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)





# PLOT 12: NNC expression genes vs not NNC expression genes
# NNC expression genes vs not NNC expression genes
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$gene_exp))) {
        if (nrow(data.class.list[[i]][data.class.list[[i]]$structural_category=="NNC",])!=0) {
            
            NNC_genes <- unique(data.class.list[[i]][data.class.list[[i]]$structural_category=="NNC","associated_gene"])
            notNNC_genes <- unique(data.class.list[[i]][!data.class.list[[i]]$associated_gene%in%NNC_genes,"associated_gene"])
            isoPerGene.list[[i]][isoPerGene.list[[i]]$associatedGene %in% notNNC_genes, "NNC_class"] <- "Genes without\n NNC isoforms"
            isoPerGene.list[[i]][isoPerGene.list[[i]]$associatedGene %in% NNC_genes, "NNC_class"] <- "Genes with\n NNC isoforms"
            
            isoPerGene.list[[i]]$NNC_class <- factor(isoPerGene.list[[i]]$NNC_class, levels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"),
                                                                 labels=c("Genes with\n NNC isoforms","Genes without\n NNC isoforms"), order=T)
            
            p.list[[i]] <- ggplot(data=isoPerGene.list[[i]][!is.na(isoPerGene.list[[i]]$NNC_class),], aes(x=NNC_class, y=log2(geneExp+1), fill=NNC_class)) +
                           geom_boxplot(color="black", size=0.3, outlier.size=0.2) +
                           xlab("") +  
                           ylab("log2(Gene_TPM+1)") +
                           scale_x_discrete(drop=FALSE) +
                           scale_fill_manual(values = c(myPalette[4],"grey38")) +
                           guides(fill="none") +
                           mytheme +
                           theme(axis.title.x=element_blank()) + 
                           ggtitle(paste0(sample.names[i], "\n\n"))
        }
    }
}
page_title = textGrob("Gene Expression of NNC And Not NNC Containing Genes", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)



# commented out plots in generatePDFreport.R
#    if (!all(is.na(data.class$gene_exp))){
#        if (nrow(data.class[data.class$structural_category=="NNC",])!=0 & nrow(data.class[data.class$structural_category=="FSM",])!=0 ){
#            
#            FSM_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="FSM","associated_gene"])
#            NNC_just_genes = unique(data.class[data.class$FSM_class=="A" & data.class$structural_category=="NNC","associated_gene"])
#            FSMandNNCgenes = unique(data.class[data.class$FSM_class=="C" & data.class$structural_category=="NNC","associated_gene"])
#            isoPerGene[isoPerGene$associatedGene %in% FSMandNNCgenes, "FSM_NNC_class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
#            isoPerGene[isoPerGene$associatedGene %in% NNC_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only NNC isoforms"
#            isoPerGene[isoPerGene$associatedGene %in% FSM_just_genes, "FSM_NNC_class"] <- "Genes expressing\n only FSM isoforms"
#            data.class[data.class$associated_gene %in% FSMandNNCgenes, "class"] <- "Genes expressing\nboth NNC and\n FSM isoforms"
#            data.class[data.class$associated_gene %in% NNC_just_genes, "class"] <- "Genes expressing\n only NNC isoforms"
#            data.class[data.class$associated_gene %in% FSM_just_genes, "class"] <- "Genes expressing\n only FSM isoforms"
#            
#            isoPerGene$FSM_NNC_class = factor(isoPerGene$FSM_NNC_class, levels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"),
#                                                                        labels=c("Genes expressing\nboth NNC and\n FSM isoforms","Genes expressing\n only NNC isoforms","Genes expressing\n only FSM isoforms"), order=T)
#            
#            p13 <- ggplot(data=isoPerGene[!is.na(isoPerGene$FSM_NNC_class),], aes(x=FSM_NNC_class, y=log2(geneExp+1), fill=FSM_NNC_class)) +
#                geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#                ylab("log2( # Short reads per gene + 1)") +
#                theme(axis.title.x=element_blank()) +
#                #theme(plot.margin = unit(c(1.5,1,0.5,1), "cm")) +
#                scale_fill_manual(values = c("grey38",myPalette[[4]],myPalette[[1]])) +
#                guides(fill="none") +
#                mytheme +
#                theme(axis.title.x=element_blank()) +
#                ggtitle("Gene Expression Level in NNC/FSM Containing Genes\n\n" ) +
#                scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
#                                          "Genes expressing\n only FSM isoforms",
#                                          "Genes expressing\n only NNC isoforms"),
#                                 labels=c("NNC/FSM genes",
#                                          "FSM genes",
#                                          "NNC genes"), drop=FALSE) 
#            
#            
#            p13.c <- ggplot(data=data.class[!is.na(data.class$class),], aes(x=class, y=log2(iso_exp+1), fill=structural_category)) +
#                geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
#                ylab("log2( # Short reads per transcript + 1)") +
#                theme(axis.title.x=element_blank()) +
#                scale_fill_manual(values = myPalette) +
#                guides(fill="none") +
#                mytheme +
#                theme(axis.title.x=element_blank()) +
#                ggtitle("Transcript Expression Level in NNC/FSM Containing Genes\n\n" ) +
#                scale_x_discrete(breaks=c("Genes expressing\nboth NNC and\n FSM isoforms",
#                                          "Genes expressing\n only FSM isoforms",
#                                          "Genes expressing\n only NNC isoforms"),
#                                 labels=c("NNC/FSM genes",
#                                          "FSM genes",
#                                          "NNC genes"), drop=F) 
#            
#        }
#    }


s <- textGrob("Splice Junction Characterization", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)

# PLOT 23: Junction categories
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.junction.list[[i]]) > 0){
        data.junction.list[[i]]$junctionLabel = with(data.junction.list[[i]], paste(chrom, strand,genomic_start_coord, genomic_end_coord, sep="_"))
        
        data.junction.list[[i]]$canonical_known = with(data.junction.list[[i]], paste(junction_category,canonical,"SJ", sep="_"))
        data.junction.list[[i]]$canonical_known = as.factor(data.junction.list[[i]]$canonical_known)
        data.junction.list[[i]]$canonical_known = factor(data.junction.list[[i]]$canonical_known, levels=c("known_canonical_SJ", "known_non_canonical_SJ", "novel_canonical_SJ", "novel_non_canonical_SJ"),
                                                                              labels=c("Known\ncanonical ", "Known\nNon-canonical ", "Novel\ncanonical ", "Novel\nNon-canonical "), order=T) 
        data.junction.list[[i]]$structural_category = data.class.list[[i]][data.junction.list[[i]]$isoform,"structural_category"]
        ##    data.junction.list[[i]]$TSSrange =cut(data.junction.list[[i]]$transcript_coord, breaks = c(0, 40, 80, 120, 160, 200, 10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))
        
        p.list[[i]] <- ggplot(data.junction.list[[i]], aes(x=structural_category)) +
                       geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=SJ_type), color="black",  size=0.3, width = 0.7) +
                       scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                                          labels=c("0","25","50","75","100"), expand = c(0,0)) +
                       scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
                       ylab("Splice junctions, %") +
                       mytheme +
                       guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
                       theme(legend.position="bottom", legend.title=element_blank())  +
                       theme(axis.text.x = element_text(angle = 45)) +
                       theme(axis.text.x = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                       theme(axis.title.x = element_blank()) +
                       ggtitle(paste0(sample.names[i],"\n\n\n"))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distribution of Splice Junctions by Structural Classification", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    t <- subset(data.class.list[[i]], exons > 1)  # select only multi-exon isoforms
    p.list[[i]] <- ggplot(data=t, aes(x=structural_category)) +
                   geom_bar(position="fill", aes(y = (..count..)/sum(..count..), fill=all_canonical), color="black", size=0.3, width = 0.7) +
                   scale_y_continuous(breaks=c(0.0,0.25,0.5,0.75,1),
                                      labels=c("0","25","50","75","100"), expand = c(0,0)) +
                   scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
                   xlab("") +
                   ylab("Transcripts, %") +
                   mytheme +
                   guides(fill = guide_legend(keywidth = 0.7, keyheight = 0.3))+
                   theme(legend.position="bottom", legend.title=element_blank())  +
                   theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   theme(axis.title.x=element_blank()) +
                   ggtitle(paste0(sample.names[i], "\n\n\n"))
}
page_title = textGrob("Distribution of Transcripts by Splice Junctions", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)


# PLOT pn4-5: Splice Junction Coverage (if coverage provided)
p.list = vector("list", length(class.files))
p.list2 = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (!all(is.na(data.junction.list[[i]]$total_coverage_unique))) {
        
        uniqJuncCov <- unique(data.junction.list[[i]][,c("junctionLabel","SJ_type", "total_coverage_unique")])
        
        e <- data.frame(table(uniqJuncCov$SJ_type))
        f <- data.frame(table(uniqJuncCov[which(uniqJuncCov$total_coverage_unique>0),"SJ_type"]))
        
        
        df.juncSupport <- data.frame(type=e$Var1, count=e$Freq-f$Freq, name='Unsupported')
        df.juncSupport <- rbind(df.juncSupport, data.frame(type=f$Var1, count=f$Freq, name='Supported'))
        
        p.list[[i]] <- ggplot(df.juncSupport, aes(x=type, y=count, fill=name)) +
                       geom_bar(stat='identity', color="black", size=0.3, width=0.7) +
                       scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
                       scale_y_continuous( expand = c(0,0)) +
                       labs(x='', y='Junctions, count', title=paste0(sample.names[i], '\n\n\n')) +
                       mytheme +
                       theme(legend.position="bottom", legend.title=element_blank()) +
                       guides(fill = guide_legend(title = "") )
        
        
        df.SJcov <- merge(e, f, by="Var1")
        # calculate the percentage of junctions that have zero short read junction coverage
        df.SJcov$perc <- 100-df.SJcov$Freq.y / df.SJcov$Freq.x * 100 ;
        df.SJcov[is.na(df.SJcov$perc), "perc"] <- 0
        maxP=min(100, (max(df.SJcov$perc) %/% 5) * 5) + 5 ;
        p.list2[[i]] <- ggplot(df.SJcov, aes(x=Var1,fill=Var1, y=perc)) +
                        geom_bar(stat="identity", position = position_dodge(), color="black", size=0.3, width=0.7) +
                        geom_text(label=paste(round(df.SJcov$perc,1.05),"%",sep=''), nudge_y=1.5) +
                        scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
                        scale_y_continuous( expand = c(0,0), limits = c(0,maxP)) +
                        labs(x='', y="Junctions, %", title=paste0(sample.names[i], '\n\n\n')) +
                        mytheme +
                        guides(fill="none")
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Unique Junctions w/ or w/out Short Read Coverage", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    page_title = textGrob("Unique Junctions w/out Short Read Coverage", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list2)
}
rm(p.list)
rm(p.list2)
# rm(uniqJuncCov)
# rm(dj.SJcov)
# rm(e)
# rm(f)



    
# PLOT 29: RT-switching
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (sum(data.junction.list[[i]]$RTS_junction=='TRUE') > 0) {
        
        a <- data.frame(table(data.junction.list[[i]]$SJ_type));
        b <- data.frame(table(subset(data.junction.list[[i]], RTS_junction=="TRUE")$SJ_type));
        
        df.RTS <- merge(a, b, by="Var1");
        df.RTS$perc <- df.RTS$Freq.y/df.RTS$Freq.x *100
        df.RTS[is.na(df.RTS$perc),"perc"] <- 0
        
        maxH <- min(100, (max(df.RTS$perc) %/% 5) * 5 + 5);
        
        p.list[[i]] <- ggplot(data=df.RTS, aes(x=Var1, y=perc, fill=Var1)) +
                       geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
                       geom_text(label=paste(round(df.RTS$perc, 2),"%",sep=''), nudge_y=0.3) +
                       scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
                       labs(x="", y="RT-switching junctions, %") +
                       ggtitle(paste0(sample.names[i],"\n\n" )) +
                       mytheme +
                       guides(fill="none") +
                       scale_y_continuous(expand = c(0,0), limits = c(0,maxH)) +
                       theme(axis.text.x = element_text(size=11))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("RT-Switching All Junctions", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    rm(df.RTS)
    rm(a)
    rm(b)
}
rm(p.list)



p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (sum(data.junction.list[[i]]$RTS_junction=='TRUE') > 0) {
        c <- data.frame(table(uniqJuncRTS.list[[i]]$SJ_type));
        d <- data.frame(table(subset(uniqJuncRTS.list[[i]], RTS_junction=='TRUE')$SJ_type));

        df.uniqRTS <- merge(c, d, by="Var1");
        df.uniqRTS$perc <- df.uniqRTS$Freq.y/df.uniqRTS$Freq.x *100
        df.uniqRTS[is.na(df.uniqRTS$perc),"perc"] <- 0
        
        maxH <- min(100, (max(df.uniqRTS$perc) %/% 5) * 5 + 5);
        
        p.list[[i]] <- ggplot(data=df.uniqRTS, aes(x=Var1, y=perc, fill=Var1)) +
                       geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
                       geom_text(label=paste(round(df.uniqRTS$perc, 2),"%",sep=''), nudge_y=0.3) +
                       scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=F) +
                       labs(x="", y="RT-switching junctions, %") +
                       ggtitle(paste0(sample.names[i],"\n\n" )) +
                       mytheme +
                       guides(fill="none") +
                       scale_y_continuous(expand = c(0,0), limits = c(0,maxH)) +
                       theme(axis.text.x = element_text(size=11))
        
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Unique Junctions RT-switching", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    rm(df.uniqRTS)
    rm(c)
    rm(d)
    rm(uniqJuncRTS.list)
    n.data.junction.list = lapply(data.junction.list, nrow)
    rm(data.junction.list)
}
rm(p.list)



s <- textGrob("Comparison With Annotated TSS and TTS", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)

alpha_TTS.list = vector("character", length(class.files))
alpha_TTS_labs.list = vector("character", length(class.files))
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.FSM.list[[i]]) > 0) {
        diff_max <- max(max(abs(data.FSM.list[[i]]$diff_to_TSS)), max(abs(data.FSM.list[[i]]$diff_to_TTS)));
        diff_breaks <- c(-(diff_max+1), seq(-1000, 1000, by = 100), diff_max+1);
        breaks_labels <- c("Larger than -1 kb", "-1 to -0.9 kb", "-0.9 to -0.8 kb", "-0.8 to -0.7 kb", "-0.7 to -0.6 kb", "-0.6 to -0.5 kb",
                           "-0.5 to -0.4 kb", "-0.4 to -0.3 kb", "-0.3 to -0.2 kb", "-0.2 to -0.1 kb", "-0.1 to 0 kb", "0 to 0.1 kb", "0.1 to 0.2 kb",
                           "0.2 to 0.3 kb", "0.3 to 0.4 kb", "0.4 to 0.5 kb", "0.5 to 0.6 kb", "0.6 to 0.7 kb", "0.7 to 0.8 kb", "0.8 to 0.9 kb", "0.9 to 1 kb",
                           "Larger than 1 kb")
        data.FSM.list[[i]]$diffTTSCat = cut(-(data.FSM.list[[i]]$diff_to_TTS), breaks = diff_breaks);
        data.FSM.list[[i]]$diffTSSCat = cut(-(data.FSM.list[[i]]$diff_to_TSS), breaks = diff_breaks);

        # plot histogram of distance to polyA site, Y-axis absolute count
        
        if (!all(is.na(data.FSM.list[[i]]$polyA_motif))) {
            alpha_TTS.list[[i]] = '!is.na(polyA_motif)'
            alpha_TTS_labs.list[[i]] = "polyA motif found"
        } else {
            alpha_TTS.list[[i]] = NULL
            alpha_TTS_labs.list[[i]] = NULL
        }

        max_height <- max(table(data.FSM.list[[i]]$diffTTSCat));
        max_height <- (max_height %/% 10+1) * 10;
        p.list[[i]] <- ggplot(data=data.FSM.list[[i]], aes(x=diffTTSCat)) +
                       geom_bar(fill=myPalette[4], color="black", size=0.3, aes(alpha=eval(parse(text = alpha_TTS.list[[i]])))) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                       mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       ylab("Transcripts, count")+
                       xlab("Distance to annotated Transcription Termination Site (TTS), bp")+
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate upstream of annotated termination site\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to annotated Transcription Termination Site (TTS)\nFSM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.FSM.list[[i]]) > 0) {
        p.list[[i]] <- ggplot(data=data.FSM.list[[i]], aes(x=diffTTSCat)) +
                       geom_bar(aes(y = (..count..)/sum(..count..) , alpha= eval(parse(text = alpha_TTS.list[[i]]))), fill=myPalette[4], color="black", size=0.3)+
                       scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                          labels=c("0","25","50","75","100"), expand = c(0,0)) +
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
                       ylab("Transcripts, %")+
                       xlab("Distance to annotated Transcription Termination Site (TTS), bp")+
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate upstream of annotated termination site\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to annotated Transcription Termination Site (TTS)\nFSM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)



alpha_TSS.list = vector("character", length(class.files))
alpha_TSS_labs.list = vector("character", length(class.files))
p.list = vector("list", length(class.files))
# plot histogram of distance to start site, Y-axis absolute count
for (i in seq_along(class.files)) {
    if (nrow(data.FSM.list[[i]]) > 0) {
        max_height <- max(table(data.FSM.list[[i]]$diffTSSCat));
        max_height <- (max_height %/% 10+1) * 10;
        if (!all(is.na(data.FSM.list[[i]]$within_CAGE_peak))) {
            alpha_TSS.list[[i]] = 'within_CAGE_peak'
            alpha_TSS_labs.list[[i]] = "TSS within a CAGE peak"
        } else {
            alpha_TSS.list[[i]] = NULL
            alpha_TSS_labs.list[[i]] = NULL
        }
        
        p.list[[i]] <- ggplot(data=data.FSM.list[[i]], aes(x=diffTSSCat)) +
                       geom_bar(fill=myPalette[6], color="black", size=0.3 , aes(alpha=eval(parse(text = alpha_TSS.list[[i]]))))+
                       scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       mytheme + labs(alpha=alpha_TSS_labs.list[[i]]) +
                       ylab("Transcripts, count") +
                       xlab("Distance to Annotated Transcription Start Site (TSS), bp") +
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to Annotated Transcription Start Site for FSM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.FSM.list[[i]]) > 0) {
        p.list[[i]] <- ggplot(data=data.FSM.list[[i]], aes(x=diffTSSCat)) +
                       geom_bar(aes(y = (..count..)/sum(..count..) , alpha=eval(parse(text = alpha_TSS.list[[i]]))), fill=myPalette[6], color="black", size=0.3)+
                       scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                          labels=c("0","25","50","75","100"), expand = c(0,0)) +
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       mytheme + labs(alpha=alpha_TSS_labs.list[[i]]) +
                       ylab("Transcripts, %")+
                       xlab("Distance to Annotated Transcription Start Site (TSS), bp")+
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to Annotated Transcription Start Site for FSM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)



breaks_labels <- c("Larger than -5 kb", "-5 to -4.5 kb", "-4.5 to -4 kb", "-4 to -3.5 kb", "-3.5 to -3 kb", "-3 to -2.5 kb",
                   "-2.5 to -2 kb", "-2 to -1.5 kb", "-1.5 to -1 kb", "-1 to -0.5 kb", "-0.5 to 0 kb", "0 to 0.5 kb", "0.5 to 1 kb",
                   "1 to 1.5 kb", "1.5 to 2 kb", "2 to 2.5 kb", "2.5 to 3 kb", "3 to 3.5 kb", "3.5 to 4 kb", "4 to 4.5 kb", "4.5 to 5 kb", "Larger than 5 kb")
max_height.list = vector("double", length(class.files))
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.ISM.list[[i]]) > 0) {
        diff_max <- max(max(abs(data.ISM.list[[i]]$diff_to_TSS)), max(abs(data.ISM.list[[i]]$diff_to_TTS)));
        diff_breaks <- c(-(diff_max+1), seq(-5000, 5000, by = 500), diff_max+1);
        
        data.ISM.list[[i]]$diffTTSCat = cut(-(data.ISM.list[[i]]$diff_to_TTS), breaks = diff_breaks);
        data.ISM.list[[i]]$diffTSSCat = cut(-(data.ISM.list[[i]]$diff_to_TSS), breaks = diff_breaks);
        
        max_height <- max(table(data.ISM.list[[i]]$diffTTSCat));
        max_height.list[[i]] <- (max_height %/% 10+1) * 10;
        
        p.list[[i]] <- ggplot(data=data.ISM.list[[i]], aes(x=diffTTSCat)) +
                       geom_bar(fill=myPalette[4], color="black", size=0.3,  aes( alpha= !is.na(polyA_motif))) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,max_height.list[[i]]))+
                       mytheme + labs(alpha = "polyA motif found") +
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       ylab("Transcripts, count")+
                       xlab("Distance to annotated polyadenylation site, bp")+
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to Annotated Polyadenylation Site for ISM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


# plot histogram of distance to polyA site, Y-axis percentages
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.ISM.list[[i]]) > 0) {
        p.list[[i]] <- ggplot(data=data.ISM.list[[i]], aes(x=diffTTSCat)) +
                       geom_bar(aes(y = (..count..)/sum(..count..), alpha= !is.na(polyA_motif)), fill=myPalette[4], color="black", size=0.3)+
                       scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                          labels=c("0","25","50","75","100"), expand = expansion(mult=c(0,0.1))) +
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       mytheme + labs(alpha = "polyA motif found") +
                       ylab("Transcripts, %")+
                       xlab("Distance to annotated polyadenylation site, bp")+
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to Annotated Polyadenylation Site for ISM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)



# plot histogram of distance to start site, Y-axis absolute count
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.ISM.list[[i]]) > 0) {

        max_height <- max(table(data.ISM.list[[i]]$diffTSSCat));
        max_height.list[[i]] <- (max_height %/% 10+1) * 10;

        p.list[[i]] <- ggplot(data=data.ISM.list[[i]], aes(x=diffTSSCat)) +
                       geom_bar(fill=myPalette[6], color="black", size=0.3, aes(alpha=within_CAGE_peak))+
                       scale_y_continuous(expand = c(0,0), limits = c(0,max_height.list[[i]]))+
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       mytheme + labs(alpha="TSS within a CAGE peak") +
                       ylab("Transcripts, count")+
                       xlab("Distance to annotated transcription start site, bp")+
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to Annotated Transcription Start Site for ISM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


# plot histogram of distance to start site, Y-axis absolute count
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.ISM.list[[i]]) > 0) {
        p.list[[i]] <- ggplot(data=data.ISM.list[[i]], aes(x=diffTSSCat)) +
                       geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), fill=myPalette[6], color="black", size=0.3)+
                       scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                          labels=c("0","25","50","75","100"), expand = c(0,0)) +
                       scale_x_discrete(drop=F, labels=breaks_labels) +
                       mytheme + labs(alpha="TSS within a CAGE peak") +
                       ylab("Transcripts, %")+
                       xlab("Distance to annotated transcription start site, bp")+
                       labs(title=paste0(sample.names[[i]], "\n\n"),
                            subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                       theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                       theme(legend.justification=c(1,1), legend.position=c(1,1))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance to Annotated Transcription Start Site for ISM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)




s <- textGrob("Comparison With Annotated TSS and TTS \nby Subcategories", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)

p21.stitles.FSM = list(textGrob("Distance to annotated Transcription Termination Site (TTS)\nFSM Alternative 3'End", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to annotated Transcription Termination Site (TTS)\nFSM Alternative 3'5'End", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to annotated Transcription Termination Site (TTS)\nFSM Alternative 5'End", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to annotated Transcription Termination Site (TTS)\nFSM Reference Match", gp=gpar(fontface="italic", fontsize=17)))
# p.list = vector("list", length(p21.stitles.FSM))
# p.list2 = vector("list", length(p21.stitles.FSM))
# for (j in seq_along(p21.stitles.FSM)) {
#     p.list[[j]] = vector("list", length(class.files))
#     p.list2[[j]] = vector("list", length(class.files))
# }
# for (i in seq_along(class.files)) {
#     if (nrow(data.FSM.list[[i]]) > 0) {
#         if (!all(is.na(data.FSM.list[[i]]$polyA_motif))) {
#             for(j in 1:length(subcategories.FSM.list[[i]])) {
#                 c <- data.frame(subcategories.FSM.list[[i]][[j]])
#                 if (!(dim(c))[1]==0 & !all(is.na(c$polyA_motif))) {
#                     diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
#                     diff_breaks <- c(-(diff_max+1), seq(-1000, 1000, by = 100), (diff_max+1));
#                     c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
#                     max_height <- max(table(c$diffTTSCat));
#                     max_height <- (max_height %/% 10+1) * 10;
#                     #p.list[[j]][[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
#                     #                    geom_bar(fill=myPalette[4], color="black", size=0.3, aes( alpha=eval(parse(text = alpha_TTS.list[[i]])))) +
#                     #                    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
#                     #                    mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
#                     #                    scale_x_discrete(drop=F, labels=breaks_labels) +
#                     #                    ylab("Transcripts, count")+
#                     #                    xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
#                     #                    labs(title=sample.names[[i]],
#                     #                         subtitle="Negative values indicate upstream of annotated termination site\n\n") +
#                     #                    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                     #                    theme(legend.justification=c(1,1), legend.position=c(1,1))
#                     #p.list2[[j]][[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
#                     #                     geom_bar(aes(alpha=eval(parse(text = alpha_TTS.list[[i]])), y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
#                     #                     scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
#                     #                                        labels=c("0","25","50","75","100"), expand=c(0,0)) +
#                     #                     scale_x_discrete(drop=F, labels=breaks_labels) +
#                     #                     mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
#                     #                     ylab("Transcripts, %")+
#                     #                     xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
#                     #                     labs(title=sample.names[[i]],
#                     #                          subtitle="Negative values indicate upstream of annotated termination site\n\n") +
#                     #                     theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                     #                     theme(legend.justification=c(1,1), legend.position=c(1,1))
#                     p.list[[j]][[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
#                                         geom_bar(fill=myPalette[4], color="black", size=0.3, aes( alpha=polyA_motif_found)) +
#                                         scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
#                                         mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
#                                         scale_x_discrete(drop=F, labels=breaks_labels) +
#                                         ylab("Transcripts, count")+
#                                         xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
#                                         labs(title=sample.names[[i]],
#                                              subtitle="Negative values indicate upstream of annotated termination site\n\n") +
#                                         theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                         theme(legend.justification=c(1,1), legend.position=c(1,1))
#                     p.list2[[j]][[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
#                                          geom_bar(aes(alpha=polyA_motif_found, y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
#                                          scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
#                                                             labels=c("0","25","50","75","100"), expand=c(0,0)) +
#                                          scale_x_discrete(drop=F, labels=breaks_labels) +
#                                          mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
#                                          ylab("Transcripts, %")+
#                                          xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
#                                          labs(title=sample.names[[i]],
#                                               subtitle="Negative values indicate upstream of annotated termination site\n\n") +
#                                          theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                          theme(legend.justification=c(1,1), legend.position=c(1,1))
#                 }
#             }
#         }
#     }
# }
# if (!all(sapply(unlist(p.list), is.null))) {
#     for (j in seq_along(p21.stitles.FSM)) {
#         organize_in_grid_with_title(p21.stitles.FSM[[j]], p.list[[j]])
#         organize_in_grid_with_title(p21.stitles.FSM[[j]], p.list2[[j]])
#     }
# }
# rm(p.list)
# rm(p.list2)




for (j in seq_along(p21.stitles.FSM)) {
    p.list = vector("list", length(class.files))
    p.list2 = vector("list", length(class.files))

    for (i in seq_along(class.files)) {
        if ((nrow(data.FSM.list[[i]]) > 0) & (!all(is.na(data.FSM.list[[i]]$polyA_motif)))) {
            c <- data.frame(subcategories.FSM.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$polyA_motif))) {
                diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
                diff_breaks <- c(-(diff_max+1), seq(-1000, 1000, by = 100), (diff_max+1));
                c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
                max_height <- max(table(c$diffTTSCat));
                max_height <- (max_height %/% 10+1) * 10;
                p.list[[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
                                    geom_bar(fill=myPalette[4], color="black", size=0.3, aes( alpha=polyA_motif_found)) +
                                    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                                    mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
                                    scale_x_discrete(drop=F, labels=breaks_labels) +
                                    ylab("Transcripts, count")+
                                    xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
                                    labs(title=sample.names[[i]],
                                         subtitle="Negative values indicate upstream of annotated termination site\n\n") +
                                    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                                    theme(legend.justification=c(1,1), legend.position=c(1,1))
                p.list2[[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
                                     geom_bar(aes(alpha=polyA_motif_found, y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3)+
                                     scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                        labels=c("0","25","50","75","100"), expand=c(0,0)) +
                                     scale_x_discrete(drop=F, labels=breaks_labels) +
                                     mytheme + labs(alpha = alpha_TTS_labs.list[[i]]) +
                                     ylab("Transcripts, %")+
                                     xlab("Distance to Annotated Transcription Termination Site (TTS), bp")+
                                     labs(title=sample.names[[i]],
                                          subtitle="Negative values indicate upstream of annotated termination site\n\n") +
                                     theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                                     theme(legend.justification=c(1,1), legend.position=c(1,1))

            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        organize_in_grid_with_title(p21.stitles.FSM[[j]], p.list)
        organize_in_grid_with_title(p21.stitles.FSM[[j]], p.list2)
    }
    rm(p.list)
    rm(p.list2)
    gc()
}





p21.stitles.ISM = list(textGrob("Distance to Annotated Polyadenylation Site for ISM\n3' Fragment", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Polyadenylation Site for ISM\nInternal Fragment", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Polyadenylation Site for ISM\nA5' Fragment", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Polyadenylation Site for ISM\nIntron Retention", gp=gpar(fontface="italic", fontsize=17)))
# p.list = vector("list", length(p21.stitles.ISM))
# p.list2 = vector("list", length(p21.stitles.ISM))
# for (j in seq_along(p21.stitles.ISM)) {
#     p.list[[j]] = vector("list", length(class.files))
#     p.list2[[j]] = vector("list", length(class.files))
# }
# for (i in seq_along(class.files)) {
#     if (nrow(data.ISM.list[[i]]) > 0) {
#         for (j in 1:length(subcategories.ISM.list[[i]])) {
#             c <- data.frame(subcategories.ISM.list[[i]][[j]])
#             if (!(dim(c))[1]==0 & !all(is.na(c$polyA_motif))) {
#                 diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
#                 diff_breaks <- c(-(diff_max+1), seq(-5000, 5000, by = 500), diff_max+1);
#                 c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
#                 max_height <-  max(table(c$diffTTSCat));
#                 max_height <- (max_height %/% 10+1) * 10;
#                 p.list[[j]][[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
#                                     geom_bar(fill=myPalette[4], color="black", size=0.3, aes( alpha= !is.na(polyA_motif))) +
#                                     scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
#                                     scale_x_discrete(drop=F, labels=breaks_labels) +
#                                     mytheme + labs(alpha = "polyA motif found") +
#                                     ylab("Transcripts, count") +
#                                     xlab("Distance to annotated polyadenylation site, bp") +
#                                     labs(title=sample.names[[i]],
#                                          subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
#                                     theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                     theme(legend.justification=c(1,1), legend.position=c(1,1))
#                 p.list2[[j]][[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
#                                      geom_bar(aes(alpha= !is.na(polyA_motif), y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3) +
#                                      scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
#                                                         labels=c("0","25","50","75","100"), expand = c(0,0)) +
#                                      scale_x_discrete(drop=F, labels=breaks_labels) +
#                                      mytheme + labs(alpha = "polyA motif found") +
#                                      ylab("Transcripts, %") +
#                                      xlab("Distance to annotated polyadenylation site, bp") +
#                                      labs(title=sample.names[[i]],
#                                           subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
#                                      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                      theme(legend.justification=c(1,1), legend.position=c(1,1))
#             }
#         }
#     }
# }
# if (!all(sapply(unlist(p.list), is.null))) {
#     for (j in seq_along(p21.stitles.ISM)) {
#         organize_in_grid_with_title(p21.stitles.ISM[[j]], p.list[[j]])
#         organize_in_grid_with_title(p21.stitles.ISM[[j]], p.list2[[j]])
#     }
# }
# rm(p.list)
# rm(p.list2)




for (j in seq_along(p21.stitles.ISM)) {
    p.list = vector("list", length(class.files))
    p.list2 = vector("list", length(class.files))

    for (i in seq_along(class.files)) {
        if (nrow(data.ISM.list[[i]]) > 0) {
            c <- data.frame(subcategories.ISM.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$polyA_motif))) {
                diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
                diff_breaks <- c(-(diff_max+1), seq(-5000, 5000, by = 500), diff_max+1);
                c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
                max_height <-  max(table(c$diffTTSCat));
                max_height <- (max_height %/% 10+1) * 10;
                p.list[[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
                                    geom_bar(fill=myPalette[4], color="black", size=0.3, aes( alpha= polyA_motif_found)) +
                                    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                                    scale_x_discrete(drop=F, labels=breaks_labels) +
                                    mytheme + labs(alpha = "polyA motif found") +
                                    ylab("Transcripts, count") +
                                    xlab("Distance to annotated polyadenylation site, bp") +
                                    labs(title=sample.names[[i]],
                                         subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
                                    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                                    theme(legend.justification=c(1,1), legend.position=c(1,1))
                p.list2[[i]] <- ggplot(data=c, aes(x=diffTTSCat)) +
                                     geom_bar(aes(alpha= polyA_motif_found, y = (..count..)/sum(..count..)), fill=myPalette[4], color="black", size=0.3) +
                                     scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                        labels=c("0","25","50","75","100"), expand = c(0,0)) +
                                     scale_x_discrete(drop=F, labels=breaks_labels) +
                                     mytheme + labs(alpha = "polyA motif found") +
                                     ylab("Transcripts, %") +
                                     xlab("Distance to annotated polyadenylation site, bp") +
                                     labs(title=sample.names[[i]],
                                          subtitle="Negative values indicate upstream of annotated polyA site\n\n") +
                                     theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                                     theme(legend.justification=c(1,1), legend.position=c(1,1))
            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        organize_in_grid_with_title(p21.stitles.ISM[[j]], p.list)
        organize_in_grid_with_title(p21.stitles.ISM[[j]], p.list2)
    }
    rm(p.list)
    rm(p.list2)
    gc()
}




#FSM_TSS
p22.stitles.FSM = list(textGrob("Distance to Annotated Transcription Start Site\nFSM Alternative 3' End", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Transcription Start Site\nFSM Alternative 3'5' End", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Transcription Start Site\nFSM Alternative 5' End", gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Transcription Start Site\nFSM Reference Match", gp=gpar(fontface="italic", fontsize=17)))
# p.list = vector("list", length(p22.stitles.FSM))
# p.list2 = vector("list", length(p22.stitles.FSM))
# for (j in seq_along(p22.stitles.FSM)) {
#     p.list[[j]] = vector("list", length(class.files))
#     p.list2[[j]] = vector("list", length(class.files))
# }
# 
# for (i in seq_along(class.files)) {
#     if (nrow(data.FSM.list[[i]]) > 0 && !all(is.na(data.FSM.list[[i]]$within_cage_peak))) {
#         for(j in 1:length(subcategories.FSM.list[[i]])) {
#             c <- data.frame(subcategories.FSM.list[[i]][[j]])
#             if (!(dim(c))[1]==0 & !all(is.na(c$within_CAGE_peak))) {
#                 diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
#                 diff_breaks <- c(-(diff_max+1), seq(-1000, 1000, by = 100), (diff_max+1));
#                 c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
#                 c$diffTSSCat = cut(-(c$diff_to_TSS), breaks = diff_breaks);
#                 max_height <- max(table(c$diffTSSCat));
#                 max_height <- (max_height %/% 10+1) * 10;
#                 p.list[[j]][[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
#                                     geom_bar(fill=myPalette[6], color="black", size=0.3, aes(alpha=eval(parse(text = alpha_TSS.list[[i]])))) +
#                                     scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
#                                     scale_x_discrete(drop=F, labels=breaks_labels) +
#                                     mytheme + labs(alpha=alpha_TSS_labs.list[[i]]) +
#                                     ylab("Transcripts, count")+
#                                     xlab("Distance to annotated Transcription Start Site (TSS), bp")+
#                                     labs(title=sample.names[[i]],
#                                          subtitle="Negative values indicate downstream of annotated TSS\n\n") +
#                                     theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                     theme(legend.justification=c(1,1), legend.position=c(1,1))
#                 p.list2[[j]][[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
#                                      geom_bar(aes(alpha=eval(parse(text = alpha_TSS.list[[i]])), y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
#                                      scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
#                                                                              labels=c("0","25","50","75","100"), expand = c(0,0)) +
#                                      scale_x_discrete(drop=F, labels=breaks_labels) +
#                                      mytheme + labs(alpha=alpha_TSS_labs.list[[i]]) +
#                                      ylab("Transcripts, %")+
#                                      xlab("Distance to annotated Transcription Start Site (TSS), bp")+
#                                      labs(title=sample.names[[i]],
#                                           subtitle="Negative values indicate downstream of annotated TSS\n\n") +
#                                      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                      theme(legend.justification=c(1,1), legend.position=c(1,1))
#             }
#         }
#     }
# }
# if (!all(sapply(unlist(p.list), is.null))) {
#     for (j in seq_along(p22.stitles.FSM)) {
#         organize_in_grid_with_title(p22.stitles.FSM[[j]], p.list[[j]])
#         organize_in_grid_with_title(p22.stitles.FSM[[j]], p.list2[[j]])
#     }
# }
# rm(p.list)
# rm(p.list2)




for (j in seq_along(p22.stitles.FSM)) {
    p.list = vector("list", length(class.files))
    p.list2 = vector("list", length(class.files))

    for (i in seq_along(class.files)) {
        if (nrow(data.FSM.list[[i]]) > 0 && !all(is.na(data.FSM.list[[i]]$within_cage_peak))) {
            c <- data.frame(subcategories.FSM.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$within_CAGE_peak))) {
                diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
                diff_breaks <- c(-(diff_max+1), seq(-1000, 1000, by = 100), (diff_max+1));
                c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
                c$diffTSSCat = cut(-(c$diff_to_TSS), breaks = diff_breaks);
                max_height <- max(table(c$diffTSSCat));
                max_height <- (max_height %/% 10+1) * 10;
                p.list[[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
                               geom_bar(fill=myPalette[6], color="black", size=0.3, aes(alpha="within_CAGE_peak")) +
                               scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                               scale_x_discrete(drop=F, labels=breaks_labels) +
                               mytheme + labs(alpha=alpha_TSS_labs.list[[i]]) +
                               ylab("Transcripts, count")+
                               xlab("Distance to annotated Transcription Start Site (TSS), bp")+
                               labs(title=sample.names[[i]],
                                    subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                               theme(legend.justification=c(1,1), legend.position=c(1,1))
                p.list2[[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
                                geom_bar(aes(alpha="within_CAGE_peak", y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3)+
                                scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                                        labels=c("0","25","50","75","100"), expand = c(0,0)) +
                                scale_x_discrete(drop=F, labels=breaks_labels) +
                                mytheme + labs(alpha=alpha_TSS_labs.list[[i]]) +
                                ylab("Transcripts, %")+
                                xlab("Distance to annotated Transcription Start Site (TSS), bp")+
                                labs(title=sample.names[[i]],
                                     subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                                theme(legend.justification=c(1,1), legend.position=c(1,1))

            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        organize_in_grid_with_title(p22.stitles.FSM[[j]], p.list)
        organize_in_grid_with_title(p22.stitles.FSM[[j]], p.list2)
    }
    rm(p.list)
    rm(p.list2)
    gc()
}





#ISM_TSS

p22.stitles.ISM = list(textGrob("Distance to Annotated Transcription Start Site for ISM\n3' Fragment" , gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Transcription Start Site for ISM\nInternal Fragment" , gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Transcription Start Site for ISM\nA5' Fragment" , gp=gpar(fontface="italic", fontsize=17)),
                       textGrob("Distance to Annotated Transcription Start Site for ISM\nIntron Retention" , gp=gpar(fontface="italic", fontsize=17)))
# p.list = vector("list", length(p22.stitles.ISM))
# p.list2 = vector("list", length(p22.stitles.ISM))
# for (j in seq_along(p22.stitles.ISM)) {
#     p.list[[j]] = vector("list", length(class.files))
#     p.list2[[j]] = vector("list", length(class.files))
# }
# for (i in seq_along(class.files)) {
#     if (nrow(data.ISM.list[[i]]) > 0 && !all(is.na(data.FSM.list[[i]]$within_cage_peak))) {
#         for(j in 1:length(subcategories.ISM.list[[i]])) {
#             c <- data.frame(subcategories.ISM.list[[i]][[j]])
#             if (!(dim(c))[1]==0 & !all(is.na(c$within_CAGE_peak))) {
#                 diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
#                 diff_breaks <- c(-(diff_max+1), seq(-5000, 5000, by = 500), diff_max+1)
#                 c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
#                 c$diffTSSCat = cut(-(c$diff_to_TSS), breaks = diff_breaks);
#                 max_height <- max(max(table(c$diffTSSCat)), max(table(c$diffTTSCat)));
#                 max_height <- (max_height %/% 10+1) * 10;
#                 p.list[[j]][[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
#                                     geom_bar(fill=myPalette[6], color="black", size=0.3, aes( alpha= within_CAGE_peak)) +
#                                     scale_y_continuous(expand = c(0,0), limits = c(0, max_height)) +
#                                     scale_x_discrete(drop=F, labels=breaks_labels) +
#                                     mytheme + labs(alpha = "TSS within a CAGE peak") +
#                                     ylab("Transcripts, count") +
#                                     xlab("Distance to annotated transcription start site, bp") +
#                                     labs(title=sample.names[[i]],
#                                          subtitle="Negative values indicate downstream of annotated TSS\n\n") +
#                                     theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                     theme(legend.justification=c(1,1), legend.position=c(1,1))
#                 p.list2[[j]][[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
#                                      geom_bar(aes( alpha= within_CAGE_peak, y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3) +
#                                      scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
#                                                         labels=c("0","25","50","75","100"), expand = c(0,0)) +
#                                      scale_x_discrete(drop=F, labels=breaks_labels) +
#                                      mytheme + labs(alpha = "TSS within a CAGE peak") +
#                                      ylab("Transcripts, %") +
#                                      xlab("Distance to annotated transcription start site, bp") +
#                                      labs(title=sample.names[[i]],
#                                           subtitle="Negative values indicate downstream of annotated TSS\n\n") +
#                                      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#                                      theme(legend.justification=c(1,1), legend.position=c(1,1))
#             }
#         }
#     }
# }
# if (!all(sapply(unlist(p.list), is.null))) {
#     for (j in seq_along(p22.stitles.ISM)) {
#         organize_in_grid_with_title(p22.stitles.ISM[[j]], p.list[[j]])
#         organize_in_grid_with_title(p22.stitles.ISM[[j]], p.list2[[j]])
#     }
# }
# rm(p.list)
# rm(p.list2)




for (j in seq_along(p22.stitles.ISM)) {
    p.list = vector("list", length(class.files))
    p.list2 = vector("list", length(class.files))

    for (i in seq_along(class.files)) {
        if (nrow(data.ISM.list[[i]]) > 0 && !all(is.na(data.ISM.list[[i]]$within_cage_peak))) {
            c <- data.frame(subcategories.ISM.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$within_CAGE_peak))) {
                diff_max <- max(max(abs(c$diff_to_TSS)), max(abs(c$diff_to_TTS)));
                diff_breaks <- c(-(diff_max+1), seq(-5000, 5000, by = 500), diff_max+1)
                c$diffTTSCat = cut(-(c$diff_to_TTS), breaks = diff_breaks);
                c$diffTSSCat = cut(-(c$diff_to_TSS), breaks = diff_breaks);
                max_height <- max(max(table(c$diffTSSCat)), max(table(c$diffTTSCat)));
                max_height <- (max_height %/% 10+1) * 10;
                p.list[[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
                               geom_bar(fill=myPalette[6], color="black", size=0.3, aes( alpha= within_CAGE_peak)) +
                               scale_y_continuous(expand = c(0,0), limits = c(0, max_height)) +
                               scale_x_discrete(drop=F, labels=breaks_labels) +
                               mytheme + labs(alpha = "TSS within a CAGE peak") +
                               ylab("Transcripts, count") +
                               xlab("Distance to annotated transcription start site, bp") +
                               labs(title=sample.names[[i]],
                                    subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                               theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                               theme(legend.justification=c(1,1), legend.position=c(1,1))
                p.list2[[i]] <- ggplot(data=c, aes(x=diffTSSCat)) +
                                geom_bar(aes( alpha= within_CAGE_peak, y = (..count..)/sum(..count..)), fill=myPalette[6], color="black", size=0.3) +
                                scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                   labels=c("0","25","50","75","100"), expand = c(0,0)) +
                                scale_x_discrete(drop=F, labels=breaks_labels) +
                                mytheme + labs(alpha = "TSS within a CAGE peak") +
                                ylab("Transcripts, %") +
                                xlab("Distance to annotated transcription start site, bp") +
                                labs(title=sample.names[[i]],
                                     subtitle="Negative values indicate downstream of annotated TSS\n\n") +
                                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                                theme(legend.justification=c(1,1), legend.position=c(1,1))
            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        organize_in_grid_with_title(p22.stitles.ISM[[j]], p.list)
        organize_in_grid_with_title(p22.stitles.ISM[[j]], p.list2)
    }
    rm(p.list)
    rm(p.list2)
    gc()
}








#PolyA Distance Analysis by categories
s <- textGrob("PolyA Distance Analysis", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)

p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (sum(!is.na(data.class.list[[i]]$polyA_dist)) > 10) {
    p.list[[i]] <- ggplot(data.class.list[[i]], aes(x=polyA_dist, color=structural_category)) +
                   geom_freqpoly(binwidth=1, size=1) +
                   scale_color_manual(values = cat.palette)+
                   xlab("Distance of polyA motif from 3' end, bp") +
                   ylab("Count") +
                   labs(title=sample.names[[i]]) +
                   mytheme+
                   theme(legend.title=element_blank())
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance of Detected PolyA Motif From 3' end", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


# PLOT polyA motif ranking, distance from 3' end
p.list = vector("list", length(class.files))
df.polyA_freq.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    df.polyA <- as.data.frame(group_by(data.class.list[[i]], by=structural_category) %>%
                              dplyr::summarise(count=dplyr::n(),
                                               polyA_detected=sum(!is.na(polyA_motif)),
                                               polyA_detected_perc=round(polyA_detected*100/count) ,
                                               .groups = 'keep'))
    df.polyA_freq <- as.data.frame(sort(table(data.class.list[[i]]$polyA_motif),decreasing=T))
    df.polyA_freq$perc <- round(df.polyA_freq$Freq*100/sum(df.polyA_freq$Freq),1)

    
    table.polyA <- tableGrob(df.polyA, rows = NULL, cols = c("Category","Count","polyA\nDetected","%"))
    title.polyA <- textGrob("Number of polyA Motifs Detected", gp=gpar(fontface="italic", fontsize=15), vjust=-12)
    gt.polyA <- gTree(children=gList(table.polyA, title.polyA))
    
    
    table.polyA_freq <- tableGrob(df.polyA_freq, rows = NULL, cols = c("Motif", "Count", "%"))
    title.polyA_freq <- textGrob("Frequency of PolyA Motifs", gp=gpar(fontface="italic", fontsize=15), vjust=-18)
    gt.polyA_freq <- gTree(children=gList(title.polyA_freq, table.polyA_freq))
    
    df.polyA_freq.list[[i]] = df.polyA_freq
    p.list[[i]] = arrangeGrob(gt.polyA, gt.polyA_freq, ncol=2, top=textGrob(sample.names[[i]], gp=gpar(fontsize=15)))
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance of Detected PolyA Motif From 3' end", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)



p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data.FSMISM.list[[i]], aes(x=polyA_dist, color=subcategory)) +
                   geom_freqpoly(binwidth=1, size=1) +
                   scale_color_manual(values = subcat.palette) +
                   xlab("Distance of polyA motif from 3' end, bp") +
                   ylab("Count") +
                   labs(title=sample.names[[i]]) +
                   mytheme +
                   theme(legend.title=element_blank())

}
page_title = textGrob("Distance of Detected PolyA Motif From 3'End\nby FSM and ISM Subcategories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    data.s2 = rbind(data.other.list[[i]], data.NICNNC.list[[i]])
    p.list[[i]] <- ggplot(data.s2, aes(x=polyA_dist, color=subcategory)) +
                   geom_freqpoly(binwidth=1, size=1) +
                   scale_color_manual(values = subcat.palette) +
                   xlab("Distance of polyA motif from 3' end, bp") +
                   ylab("Count") +
                   labs(title=sample.names[[i]]) +
                   mytheme +
                   theme(legend.title=element_blank())
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Distance of Detected PolyA Motif From 3'End\nby Non-FSM/ISM  Subcategories", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
# rm(data.s2)
rm(p.list)




# PLOT polyA motif ranking, distance from 3' end by subcategory
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (sum(!is.na(data.class.list[[i]]$polyA_dist)) > 10) {
 
        df.polyA_subcat <- as.data.frame(group_by(data.class.list[[i]], by=subcategory) %>%
                                         dplyr::summarise(count=dplyr::n(),
                                                          polyA_detected=sum(!is.na(polyA_motif)),
                                                          polyA_detected_perc=round(polyA_detected*100/count) ,
                                                          .groups = 'keep'))

        df.polyA_subcat <- tableGrob(df.polyA_subcat, rows = NULL, cols = c("Subcategory", "Count", "polyA\nDetected","%"))
        title.polyA <- textGrob(sample.names[[i]], gp=gpar(fontface="italic", fontsize=15), vjust=-18)
        p.list[[i]] <- gTree(children=gList(df.polyA_subcat, title.polyA))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Number of polyA Motifs Detected", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (sum(!is.na(data.class.list[[i]]$polyA_dist)) > 10) {
        table.polyA_freq <- tableGrob(df.polyA_freq.list[[i]], rows = NULL, cols = c("Motif", "Count", "%"))
        title.polyA_freq <- textGrob(sample.names[[i]], gp=gpar(fontface="italic", fontsize=15), vjust=-18)
        p.list[[i]] <- gTree(children=gList(title.polyA_freq, table.polyA_freq))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Frequency of PolyA Motifs", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)
rm(df.polyA_freq.list)






##### Distances to CAGE peaks by FSM and ISM
if (!all(is.na(data.class.list[[1]]$dist_to_CAGE_peak))) {
    s <- textGrob("CAGE Distances Analysis", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
    grid.arrange(s)


    diff_max = 11000     ##diff_max <- max(abs(data.class$dist_to_cage_peak));
    diff_breaks <- c(-(diff_max+1), seq(-200, 100, by = 20), diff_max+1);
    
    breaks_labels <- c("Larger than -200", "-200 to -180","-180 to -160","-160 to -140","-140 to -120",
                       "-120 to -100", "-100 to -80", "-80 to -60", "-60 to -40", "-40 to -20", "-20 to 0",
                       "0 to 20", "20 to 40", "40 to 60", "60 to 80", "80 to 100", "Larger than 100")
    
    for (i in seq_along(class.files)) {
        data.FSM.list[[i]]$dist_CAGE_Cat = cut(-(data.FSM.list[[i]]$dist_to_CAGE_peak), breaks = diff_breaks);
        data.ISM.list[[i]]$dist_CAGE_Cat = cut(-(data.ISM.list[[i]]$dist_to_CAGE_peak), breaks = diff_breaks);
        data.NNC.list[[i]]$dist_CAGE_Cat = cut(-(data.NNC.list[[i]]$dist_to_CAGE_peak), breaks = diff_breaks);
        data.NIC.list[[i]]$dist_CAGE_Cat = cut(-(data.NIC.list[[i]]$dist_to_CAGE_peak), breaks = diff_breaks);
    }

    #  max_height <- max(max(table(data.class$dist_cage_Cat)), max(table(data.class$dist_cage_Cat)));
    #  max_height <- (max_height %/% 10+1) * 10;
    
    #  ggplot(data=d.fsm, aes(x=dist_to_cage_peak , fill=structural_category)) +
    #    geom_density( color="black", size=0.3) + xlim(c(-50,50))+
    #   scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
    #    mytheme  + facet_wrap(~structural_category , nrow=2)+
    #   scale_x_discrete(limits = c(-50,50)) +
    #    ylab("Number of transcripts")+
    #    xlab("Distance to CAGE peak, bp")+
    #    labs(     title="Distance to CAGE Peak",
    #            subtitle="Negative values indicate downstream of annotated CAGE peak\n\n") +
    #    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    ## FSM cage hist number of Isoforms
    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.FSM.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat, fill=structural_category)) +
                      geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
                      mytheme  +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
                      theme(legend.position="bottom") +
                      ylab("Number of transcripts")+
                      xlab("Distance to CAGE peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic FSM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    rm(p.list)

    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.FSM.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                      geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
                      scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                           labels=c("0","25","50","75","100"), 
                                                           expand = expansion(mult = c(0,0.1))) +
                      mytheme  +
                      theme(legend.position="bottom") +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      ylab("Transcripts, %")+
                      xlab("Distance to CAGE peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic FSM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    rm(p.list)


    cage.titles.FSM = list(textGrob("Distance to CAGE Peak of Multi-Exonic FSM\nAlternative 3' End", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic FSM\nAlternative 3'5' End", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic FSM\nAlternative 5' End", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic FSM\nReference Match", gp=gpar(fontface="italic", fontsize=17)))
    p.list = vector("list", length(cage.titles.FSM))
    p.list2 = vector("list", length(cage.titles.FSM))
    for (j in seq_along(cage.titles.FSM)) {
        p.list[[j]] = vector("list", length(class.files))
        p.list2[[j]] = vector("list", length(class.files))
    }
    for (i in seq_along(class.files)) {
        for(j in 1:length(subcategories.FSM.list[[i]])) {
            c <- data.frame(subcategories.FSM.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))) {
                c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
                p.list[[j]][[i]] = ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                   geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
                                   mytheme  +
                                   scale_x_discrete(drop=F, labels=breaks_labels) +
                                   scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
                                   theme(legend.position="bottom") +
                                   ylab("Number of transcripts")+
                                   xlab("Distance to CAGE peak, bp")+
                                   labs(title=sample.names[[i]],
                                        subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                        alpha = "TSS within a CAGE peak") +
                                   theme(axis.text.x = element_text(angle = 60, hjust = 1))
                p.list2[[j]][[i]] = ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                    geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[1]) +
                                    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                       labels=c("0","25","50","75","100"), 
                                                       expand = expansion(mult = c(0,0.1))) +
                                    mytheme  +
                                    scale_x_discrete(drop=F, labels=breaks_labels) +
                                    theme(legend.position="bottom") +
                                    ylab("Transcripts, %")+
                                    xlab("Distance to CAGE peak, bp")+
                                    labs(title=sample.names[[i]],
                                         subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                         alpha = "TSS within a CAGE peak") +
                                    theme(axis.text.x = element_text(angle = 60, hjust = 1))
            }
        }
    }
    for (j in seq_along(cage.titles.FSM)) {
        if (!all(sapply(p.list[[j]], is.null))) {
            organize_in_grid_with_title(cage.titles.FSM[[j]], p.list[[j]])
            organize_in_grid_with_title(cage.titles.FSM[[j]], p.list2[[j]])
        }
    }
    rm(p.list)
    rm(p.list2)




    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.ISM.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                      geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[2]) +
                      #    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                      mytheme  +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
                      theme(legend.position="bottom") +
                      #    scale_x_discrete(limits = c(-50,50)) +
                      ylab("Number of transcripts")+
                      xlab("Distance to CAGE peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic ISM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    rm(p.list)


    p.list = vector("list", length(class.files))
        for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.ISM.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                      geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[2]) +
                      scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                           labels=c("0","25","50","75","100"), 
                                                           expand = expansion(mult = c(0,0.1))) +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      mytheme  + theme(legend.position="bottom") +
                      #    scale_x_discrete(limits = c(-50,50)) +
                      ylab("Transcripts, %")+
                      xlab("Distance to CAGE Peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic ISM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    rm(p.list)






    cage.titles.ISM = list(textGrob("Distance to CAGE Peak of Multi-Exonic ISM\n3' Fragment", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic ISM\nInternal Fragment", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic ISM\n5' Fragment", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic ISM\nIntron Retention", gp=gpar(fontface="italic", fontsize=17)))
    p.list = vector("list", length(cage.titles.ISM))
    p.list2 = vector("list", length(cage.titles.ISM))
    for (j in seq_along(cage.titles.ISM)) {
        p.list[[j]] = vector("list", length(class.files))
        p.list2[[j]] = vector("list", length(class.files))
    }

    for (i in seq_along(class.files)) {
        for(j in 1:length(subcategories.ISM.list[[i]])) {
            c <- data.frame(subcategories.ISM.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))) {
                c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
                p.list[[j]][[i]] = ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                   geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[2]) +
                                   mytheme  +
                                   scale_x_discrete(drop=F, labels=breaks_labels) +
                                   scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
                                   theme(legend.position="bottom") +
                                   ylab("Number of transcripts")+
                                   xlab("Distance to CAGE peak, bp")+
                                   labs(title=sample.names[[i]],
                                        subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                        alpha = "TSS within a CAGE peak") +
                                   theme(axis.text.x = element_text(angle = 60, hjust = 1))
                p.list2[[j]][[i]] = ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                    geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[2]) +
                                    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                       labels=c("0","25","50","75","100"), 
                                                       expand = expansion(mult = c(0,0.1))) +
                                    mytheme  +
                                    scale_x_discrete(drop=F, labels=breaks_labels) +
                                    theme(legend.position="bottom") +
                                    ylab("Transcripts, %")+
                                    xlab("Distance to CAGE peak, bp")+
                                    labs(title=sample.names[[i]],
                                         subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                         alpha = "TSS within a CAGE peak") +
                                    theme(axis.text.x = element_text(angle = 60, hjust = 1))
            }
        }
    }
    for (j in seq_along(cage.titles.ISM)) {
        if (!all(sapply(p.list[[j]], is.null))) {
            organize_in_grid_with_title(cage.titles.ISM[[j]], p.list[[j]])
            organize_in_grid_with_title(cage.titles.ISM[[j]], p.list2[[j]])
        }
    }
    rm(p.list)
    rm(p.list2)



    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.NIC.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                      geom_bar( aes(alpha=within_CAGE_peak),color="black", size=0.3, fill=myPalette[3]) +
                      #    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                      mytheme  + theme(legend.position="bottom") +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
                      ylab("Number of transcripts")+
                      xlab("Distance to CAGE peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic NIC", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)


    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.NIC.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                      geom_bar( aes(y = (..count..)/sum(..count..),alpha=within_CAGE_peak),color="black", size=0.3, fill=myPalette[3]) +
                      scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                         labels=c("0","25","50","75","100"), 
                                         expand = expansion(mult = c(0,0.1))) +
                      mytheme  + theme(legend.position="bottom") +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      #    scale_x_discrete(limits = c(-50,50)) +
                      ylab("Transcripts, %")+
                      xlab("Distance to CAGE peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic NIC", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)




    cage.titles.NIC = list(textGrob("Distance to CAGE Peak of Multi-Exonic NIC\nCombination of Annotated Junctions", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic NIC\nCombination of Annotated Splice Sites", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic NIC\nIntron Retention", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic NIC\nMono-Exon by Intron Retention", gp=gpar(fontface="italic", fontsize=17)))
    p.list = vector("list", length(cage.titles.NIC))
    p.list2 = vector("list", length(cage.titles.NIC))
    for (j in seq_along(cage.titles.NIC)) {
        p.list[[j]] = vector("list", length(class.files))
        p.list2[[j]] = vector("list", length(class.files))
    }

    for (i in seq_along(class.files)) {
        for(j in 1:length(subcategories.NIC.list[[i]])){
            c <- data.frame(subcategories.NIC.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))){
                c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
                p.list[[j]][[i]] = ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                   geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[3]) +
                                   mytheme  +
                                   scale_x_discrete(drop=F, labels=breaks_labels) +
                                   scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
                                   theme(legend.position="bottom") +
                                   ylab("Number of transcripts")+
                                   xlab("Distance to CAGE peak, bp")+
                                   labs(title=sample.names[[i]],
                                        subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                        alpha = "TSS within a CAGE peak") +
                                   theme(axis.text.x = element_text(angle = 60, hjust = 1))
                p.list2[[j]][[i]] = ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                    geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[3]) +
                                    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                       labels=c("0","25","50","75","100"), 
                                                       expand = expansion(mult = c(0,0.1))) +
                                    mytheme  +
                                    scale_x_discrete(drop=F, labels=breaks_labels) +
                                    theme(legend.position="bottom") +
                                    ylab("Transcripts, %")+
                                    xlab("Distance to CAGE peak, bp")+
                                    labs(title=sample.names[[i]],
                                         subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                         alpha = "TSS within a CAGE peak") +
                                    theme(axis.text.x = element_text(angle = 60, hjust = 1))
            }
        }
    }
    for (j in seq_along(cage.titles.NIC)) {
        if (!all(sapply(p.list[[j]], is.null))) {
            organize_in_grid_with_title(cage.titles.NIC[[j]], p.list[[j]])
            organize_in_grid_with_title(cage.titles.NIC[[j]], p.list2[[j]])
        }
    }
    rm(p.list)
    rm(p.list2)


    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.NNC.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                      geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[4]) +
                      #    scale_y_continuous(expand = c(0,0), limits = c(0,max_height))+
                      mytheme  + theme(legend.position="bottom") +
                      #    scale_x_discrete(limits = c(-50,50)) +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      ylab("Number of transcripts")+
                      xlab("Distance to CAGE peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic NNC", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)


    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        p.list[[i]] = ggplot(data=subset(data.NNC.list[[i]], !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                      geom_bar(aes(y = (..count..)/sum(..count..),alpha=within_CAGE_peak), color="black", size=0.3, fill=myPalette[4]) +
                      scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                         labels=c("0","25","50","75","100"), 
                                         expand = expansion(mult = c(0,0.1))) +
                      mytheme + theme(legend.position="bottom") +
                      scale_x_discrete(drop=F, labels=breaks_labels) +
                      #    scale_x_discrete(limits = c(-50,50)) +
                      ylab("Transcripts, %")+
                      xlab("Distance to CAGE peak, bp")+
                      labs(title=sample.names[[i]],
                           subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                           alpha = "TSS within a CAGE peak") +
                      theme(axis.text.x = element_text(angle = 60, hjust = 1))
    }
    page_title = textGrob("Distance to CAGE Peak of Multi-Exonic NNC", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)


    
    cage.titles.NNC = list(textGrob("Distance to CAGE Peak of Multi-Exonic NNC\nCombination of Annotated Junctions", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic NNC\nCombination of Annotated Splice Sites", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic NNC\nIntron Retention", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic NNC\nMono-Exon by Intron Retention", gp=gpar(fontface="italic", fontsize=17)),
                           textGrob("Distance to CAGE Peak of Multi-Exonic NNC\nAt Least One Annotated Donor/Acceptor", gp=gpar(fontface="italic", fontsize=17)))
    p.list = vector("list", length(cage.titles.NNC))
    p.list2 = vector("list", length(cage.titles.NNC))
    for (j in seq_along(cage.titles.NNC)) {
        p.list[[j]] = vector("list", length(class.files))
        p.list2[[j]] = vector("list", length(class.files))
    }
    for (i in seq_along(class.files)) {
        for(j in 1:length(subcategories.NNC.list[[i]])){
            c <- data.frame(subcategories.NNC.list[[i]][[j]])
            if (!(dim(c))[1]==0 & !all(is.na(c$dist_to_CAGE_peak))){
                c$dist_CAGE_Cat = cut(-(c$dist_to_CAGE_peak), breaks = diff_breaks);
                p.list[[j]][[i]]= ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                  geom_bar(aes(alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[4]) +
                                  mytheme  +
                                  scale_x_discrete(drop=F, labels=breaks_labels) +
                                  scale_y_continuous( expand = expansion(mult = c(0,0.1)) ) +
                                  theme(legend.position="bottom") +
                                  ylab("Number of transcripts")+
                                  xlab("Distance to CAGE peak, bp")+
                                  labs(title=sample.names[[i]],
                                       subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                       alpha = "TSS within a CAGE peak") +
                                  theme(axis.text.x = element_text(angle = 60, hjust = 1))
                p.list2[[j]][[i]] = ggplot(data=subset(c, !is.na(dist_CAGE_Cat)), aes(x=dist_CAGE_Cat , fill=structural_category)) +
                                    geom_bar(aes(y = (..count..)/sum(..count..), alpha=within_CAGE_peak), color="black", size=0.3,  fill=myPalette[4]) +
                                    scale_y_continuous(breaks=c(0.0,0.25,0.50,0.75,1),
                                                       labels=c("0","25","50","75","100"), 
                                                       expand = expansion(mult = c(0,0.1))) +
                                    mytheme  +
                                    scale_x_discrete(drop=F, labels=breaks_labels) +
                                    theme(legend.position="bottom") +
                                    ylab("Transcripts, %")+
                                    xlab("Distance to CAGE peak, bp")+
                                    labs(title=sample.names[[i]],
                                         subtitle="Negative values indicate downstream of annotated CAGE peak\n\n",
                                         alpha = "TSS within a CAGE peak") +
                                    theme(axis.text.x = element_text(angle = 60, hjust = 1))
            }
        }
    }
    for (j in seq_along(cage.titles.NNC)) {
        if (!all(sapply(p.list[[j]], is.null))) {
            organize_in_grid_with_title(cage.titles.NNC[[j]], p.list[[j]])
            organize_in_grid_with_title(cage.titles.NNC[[j]], p.list2[[j]])
        }
    }
    rm(p.list)
    rm(p.list2)



    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        df.cage <- as.data.frame(group_by(data.class.list[[i]], by=structural_category) %>%
                                 dplyr::summarise(count=dplyr::n(),
                                                  cage_detected=length(which(within_CAGE_peak)),
                                                  cage_detected_perc=round(cage_detected*100/count) , 
                                                  .groups = 'keep'))
        table.cage <- tableGrob(df.cage, rows = NULL, cols = c("Category","Count","CAGE\nDetected","%"))
        title.cage <- textGrob(sample.names[[i]], gp=gpar(fontface="italic", fontsize=15), vjust=-18)
        p.list[[i]] <- gTree(children=gList(table.cage, title.cage))
    }
    page_title = textGrob("Number of CAGE Detected", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)


    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        df.cage_subc <- as.data.frame(group_by(data.class.list[[i]], by=subcategory) %>%
                                      dplyr::summarise(count=dplyr::n(),
                                                       cage_detected=length(which(within_CAGE_peak)),
                                                       cage_detected_perc=round(cage_detected*100/count) ,
                                                       .groups = 'keep'))
        table.cage_subc <- tableGrob(df.cage_subc, rows = NULL, cols = c("Subcategory","Count","CAGE\nDetected","%"))
        title.cage_subc <- textGrob(sample.names[[i]], gp=gpar(fontface="italic", fontsize=15), vjust=-18)
        p.list[[i]] <- gTree(children=gList(table.cage_subc, title.cage_subc))
    }
    page_title = textGrob("Number of CAGE Detected", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}



s <- textGrob("Redundancy Analysis", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)


p.list = vector("list", length(class.files))
p.list2 = vector("list", length(class.files))
p.list3 = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (nrow(data.ISM.list[[i]]) > 0 || nrow(data.FSM.list[(i)]) > 0) {
        ism_per_transcript = data.ISM.list[[i]] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
        names(ism_per_transcript)[3] <- "ISM_per_tr"
        fsm_per_transcript = data.FSM.list[[i]] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
        names(fsm_per_transcript)[3] <- "FSM_per_tr"

        iso_per_knownTr = merge(x = fsm_per_transcript ,y = ism_per_transcript, by = "associated_transcript", all=T)
        iso_per_knownTr$ISM_per_tr[is.na(iso_per_knownTr$ISM_per_tr)] <- 0
        iso_per_knownTr$FSM_per_tr[is.na(iso_per_knownTr$FSM_per_tr)] <- 0

        iso_per_knownTr$total_iso = apply(iso_per_knownTr, 1, function(X) (as.integer(X[3]) + as.integer(X[5])) )

        #
        iso_per_knownTr$FSM_cat = NA
        iso_per_knownTr$FSM_cat = apply(iso_per_knownTr, 1, function(X){
             if(as.numeric(X["FSM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["FSM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr$FSM_bin = apply(iso_per_knownTr, 1, function(X){
            if(as.numeric(X["FSM_per_tr"]) < 8){
                return(as.character(X["FSM_per_tr"]))
            }else{
                return("8+")
            }})

        iso_per_knownTr$ISM_cat = NA
        iso_per_knownTr$ISM_cat = apply(iso_per_knownTr, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["ISM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr$ISM_bin = apply(iso_per_knownTr, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) < 8){
                return(as.character(X["ISM_per_tr"]))
            }else{
                return("8+")
            }})


        iso_per_knownTr$total_cat = NA
        iso_per_knownTr$total_cat = apply(iso_per_knownTr, 1, function(X){
            if(as.numeric(X["total_iso"]) == 1){
                return("Unique")
            }else if(as.numeric(X["total_iso"]) > 1){
                return("Multiple")
                }else{
                    return("NULL")
                }})

        iso_per_knownTr$total_bin=apply(iso_per_knownTr, 1, function(X){
            if(as.numeric(X["total_iso"]) < 8){
                return(as.character(X["total_iso"]))
            }else{
                return("8+")
            }})
        iso_per_knownTr$FSM_cat=factor(iso_per_knownTr$FSM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr$ISM_cat=factor(iso_per_knownTr$ISM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr$total_cat=factor(iso_per_knownTr$total_cat, levels=c("Unique", "Multiple"))

        max_y = max(table(iso_per_knownTr[which(iso_per_knownTr$FSM_cat != "NULL"), "FSM_cat"]))+10
        p.list[[i]] <- ggplot(iso_per_knownTr[which(iso_per_knownTr$FSM_cat != "NULL"),])+
            geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
            mytheme+
            geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
            scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
            scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
            labs(x="FSM per reference transcript",
                 y="Count of reference transcripts",
                 title=sample.names[[i]],
                 subtitle="Only FSM")

        max_y = max(table(iso_per_knownTr[which(iso_per_knownTr$ISM_cat != "NULL"), "ISM_cat"]))+10
        p.list2[[i]] <- ggplot(iso_per_knownTr[which(iso_per_knownTr$ISM_cat != "NULL"),])+
            geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
            mytheme+
            geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
            scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
            scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
            labs(x="ISM per reference transcript",
                 y="Count of Reference Transcripts",
                 title=sample.names[[i]],
                 subtitle="Only ISM")

        max_y = max(table(iso_per_knownTr[which(iso_per_knownTr$total_cat != "NULL"), "total_cat"]))+10
        p.list3[[i]] <- ggplot(iso_per_knownTr[which(iso_per_knownTr$total_cat != "NULL"),])+
            geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
            mytheme+
            geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
            scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
            scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
            labs(x="FSM+ISM per reference transcript",
                 y="Count of reference transcripts",
                 title=sample.names[[i]],
                 subtitle="FSM+ISM")
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Reference Transcript Redundancy", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    organize_in_grid_with_title(page_title, p.list2)
    organize_in_grid_with_title(page_title, p.list3)
    rm(iso_per_knownTr)
}
rm(p.list)
rm(p.list2)
rm(p.list3)



p.list = vector("list", length(class.files))
p.list2 = vector("list", length(class.files))
p.list3 = vector("list", length(class.files))
#### Now only with cage + isoforms
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$dist_to_cage_peak))) {
        ism_per_transcript_cage = data.ISM.list[[i]][which(data.ISM.list[[i]]$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
        names(ism_per_transcript_cage)[3] <- "ISM_per_tr"
        fsm_per_transcript_cage = data.FSM.list[[i]][which(data.FSM.list[[i]]$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n())
        names(fsm_per_transcript_cage)[3] <- "FSM_per_tr"

        iso_per_knownTr_cage = merge(x = fsm_per_transcript_cage ,y=ism_per_transcript_cage, by = "associated_transcript", all=T)
        iso_per_knownTr_cage$ISM_per_tr[is.na(iso_per_knownTr_cage$ISM_per_tr)] <- 0
        iso_per_knownTr_cage$FSM_per_tr[is.na(iso_per_knownTr_cage$FSM_per_tr)] <- 0

        iso_per_knownTr_cage$total_iso = apply(iso_per_knownTr_cage, 1, function(X) as.integer(X[3]) + as.integer(X[5]) )

        iso_per_knownTr_cage$FSM_cat = NA
        iso_per_knownTr_cage$FSM_cat = apply(iso_per_knownTr_cage, 1, function(X){
            if(as.numeric(X["FSM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["FSM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_cage$FSM_bin = apply(iso_per_knownTr_cage, 1, function(X){
            if(as.numeric(X["FSM_per_tr"]) < 8){
                return(as.character(X["FSM_per_tr"]))
            }else{
                return("8+")
            }})

        iso_per_knownTr_cage$ISM_cat = NA
        iso_per_knownTr_cage$ISM_cat = apply(iso_per_knownTr_cage, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["ISM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_cage$ISM_bin = apply(iso_per_knownTr_cage, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) < 8){
                return(as.character(X["ISM_per_tr"]))
            }else{
                return("8+")
            }})


        iso_per_knownTr_cage$total_cat = NA
        iso_per_knownTr_cage$total_cat = apply(iso_per_knownTr_cage, 1, function(X){
            if(as.numeric(X["total_iso"]) == 1){
                return("Unique")
            }else if(as.numeric(X["total_iso"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_cage$total_bin = apply(iso_per_knownTr_cage, 1, function(X){
            if(as.numeric(X["total_iso"]) < 8){
                return(as.character(X["total_iso"]))
            }else{
                return("8+")
            }})

        iso_per_knownTr_cage$FSM_cat = factor(iso_per_knownTr_cage$FSM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr_cage$ISM_cat = factor(iso_per_knownTr_cage$ISM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr_cage$total_cat = factor(iso_per_knownTr_cage$total_cat, levels=c("Unique", "Multiple"))


        max_y = max(table(iso_per_knownTr_cage[which(iso_per_knownTr_cage$FSM_cat!="NULL"), "FSM_cat"]))+10
        p.list[[i]] <- ggplot(iso_per_knownTr_cage[which(iso_per_knownTr_cage$FSM_cat!="NULL"),])+
                       geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
                       mytheme+
                       geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
                       scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
                       scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
                       labs(x="FSM per reference transcript",
                            y="Count of reference transcripts",
                            title=sample.names[[i]],
                            subtitle="Only FSM with CAGE support")

        max_y = max(table(iso_per_knownTr_cage[which(iso_per_knownTr_cage$ISM_cat!="NULL"), "ISM_cat"]))+10
        p.list2[[i]] <- ggplot(iso_per_knownTr_cage[which(iso_per_knownTr_cage$ISM_cat!="NULL"),])+
                        geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
                        mytheme+
                        geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
                        scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
                        scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
                        labs(x="ISM per reference transcript",
                             y="Count of reference transcripts",
                             title=sample.names[[i]],
                             subtitle="Only ISM with CAGE support")

        max_y = max(table(iso_per_knownTr_cage[which(iso_per_knownTr_cage$total_cat!="NULL"), "total_cat"]))+10
        p.list3[[i]] <- ggplot(iso_per_knownTr_cage[which(iso_per_knownTr_cage$total_cat!="NULL"),])+
                        geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
                        mytheme+
                        geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
                        scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
                        scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
                        labs(x="FSM+ISM per reference transcript",
                             y="Count of reference transcripts",
                             title=sample.names[[i]],
                             subtitle="FSM+ISM with CAGE support")
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Reference Transcript Redundancy", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    organize_in_grid_with_title(page_title, p.list2)
    organize_in_grid_with_title(page_title, p.list3)
    rm(iso_per_knownTr_cage)
}
rm(p.list)
rm(p.list2)
rm(p.list3)


p.list = vector("list", length(class.files))
p.list2 = vector("list", length(class.files))
p.list3 = vector("list", length(class.files))
### Now only with polyA motif = T isoforms
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$polyA_motif))) {
        ism_per_transcript_polya = data.ISM.list[[i]][which(!is.na(data.ISM.list[[i]]$polyA_motif)),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
        names(ism_per_transcript_polya)[3] <- "ISM_per_tr"
        fsm_per_transcript_polya = data.FSM.list[[i]][which(!is.na(data.FSM.list[[i]]$polyA_motif)),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
        names(fsm_per_transcript_polya)[3] <- "FSM_per_tr"

        iso_per_knownTr_polya=merge(x = fsm_per_transcript_polya ,y=ism_per_transcript_polya, by = "associated_transcript", all=T)
        iso_per_knownTr_polya$ISM_per_tr[is.na(iso_per_knownTr_polya$ISM_per_tr)] <- 0
        iso_per_knownTr_polya$FSM_per_tr[is.na(iso_per_knownTr_polya$FSM_per_tr)] <- 0
        iso_per_knownTr_polya$total_iso=apply(iso_per_knownTr_polya, 1, function(X) as.integer(X[3]) + as.integer(X[5]) )

        iso_per_knownTr_polya$FSM_cat = NA
        iso_per_knownTr_polya$FSM_cat = apply(iso_per_knownTr_polya, 1, function(X){
            if(as.numeric(X["FSM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["FSM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_polya$FSM_bin = apply(iso_per_knownTr_polya, 1, function(X){
            if(as.numeric(X["FSM_per_tr"]) < 8){
                return(as.character(X["FSM_per_tr"]))
            }else{
                return("8+")
            }})

        iso_per_knownTr_polya$ISM_cat = NA
        iso_per_knownTr_polya$ISM_cat = apply(iso_per_knownTr_polya, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["ISM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_polya$ISM_bin = apply(iso_per_knownTr_polya, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) < 8){
                return(as.character(X["ISM_per_tr"]))
            }else{
                return("8+")
            }})


        iso_per_knownTr_polya$total_cat = NA
        iso_per_knownTr_polya$total_cat = apply(iso_per_knownTr_polya, 1, function(X){
            if(as.numeric(X["total_iso"]) == 1){
                return("Unique")
            }else if(as.numeric(X["total_iso"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_polya$total_bin = apply(iso_per_knownTr_polya, 1, function(X){
            if(as.numeric(X["total_iso"]) < 8){
                return(as.character(X["total_iso"]))
            }else{
                return("8+")
            }})

        iso_per_knownTr_polya$FSM_cat = factor(iso_per_knownTr_polya$FSM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr_polya$ISM_cat = factor(iso_per_knownTr_polya$ISM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr_polya$total_cat = factor(iso_per_knownTr_polya$total_cat, levels=c("Unique", "Multiple"))

        max_y = max(table(iso_per_knownTr_polya[which(iso_per_knownTr_polya$FSM_cat!="NULL"), "FSM_cat"]))+10
        p.list[[i]] <- ggplot(iso_per_knownTr_polya[which(iso_per_knownTr_polya$FSM_cat!="NULL"),])+
                       geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
                       mytheme+
                       geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
                       scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
                       scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
                       labs(x="FSM per reference transcript",
                            y="Count of reference transcripts",
                            title=sample.names[[i]],
                            subtitle="Only FSM with a polyA motif found")

        max_y = max(table(iso_per_knownTr_polya[which(iso_per_knownTr_polya$ISM_cat!="NULL"), "ISM_cat"]))+10
        p.list2[[i]] <- ggplot(iso_per_knownTr_polya[which(iso_per_knownTr_polya$ISM_cat!="NULL"),])+
                        geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
                        mytheme+
                        geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
                        scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
                        scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
                        labs(x="ISM per reference transcript",
                             y="Count of reference transcripts",
                             title=sample.names[[i]],
                             subtitle="Only ISM with a polyA motif found")

        max_y = max(table(iso_per_knownTr_polya[which(iso_per_knownTr_polya$total_cat!="NULL"), "total_cat"]))+10
        p.list3[[i]] <- ggplot(iso_per_knownTr_polya[which(iso_per_knownTr_polya$total_cat!="NULL"),])+
                        geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
                        mytheme+
                        geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
                        scale_y_continuous(breaks = pretty_breaks(6), limits = c(0,max_y)) +
                        scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
                        labs(x="FSM+ISM per reference transcript",
                             y="Count of reference transcripts",
                             title=sample.names[[i]],
                             subtitle="FSM+ISM with a polyA motif found")
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Reference Transcript Redundancy", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    organize_in_grid_with_title(page_title, p.list2)
    organize_in_grid_with_title(page_title, p.list3)
    rm(iso_per_knownTr_polya)
}
rm(p.list)
rm(p.list2)
rm(p.list3)




p.list = vector("list", length(class.files))
p.list2 = vector("list", length(class.files))
p.list3 = vector("list", length(class.files))
#### Now with just isoforms polyA and Cage +
for (i in seq_along(class.files)) {
    if (!all(is.na(data.class.list[[i]]$polyA_motif)) && !all(is.na(data.class.list[[i]]$dist_to_cage_peak))) {
        ism_per_transcript_cage_polya = data.ISM.list[[i]][which(!is.na(data.ISM.list[[i]]$polyA_motif) & data.ISM.list[[i]]$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
        names(ism_per_transcript_cage_polya)[3] <- "ISM_per_tr"
        fsm_per_transcript_cage_polya = data.FSM.list[[i]][which(!is.na(data.FSM.list[[i]]$polyA_motif) & data.FSM.list[[i]]$within_CAGE_peak),] %>% group_by(associated_transcript, structural_category) %>% dplyr::summarize(dplyr::n(), .groups='keep')
        names(fsm_per_transcript_cage_polya)[3] <- "FSM_per_tr"

        iso_per_knownTr_cage_polya = merge(x = fsm_per_transcript_cage_polya , y=ism_per_transcript_cage_polya, by = "associated_transcript", all=T)
        iso_per_knownTr_cage_polya$ISM_per_tr[is.na(iso_per_knownTr_cage_polya$ISM_per_tr)] <- 0
        iso_per_knownTr_cage_polya$FSM_per_tr[is.na(iso_per_knownTr_cage_polya$FSM_per_tr)] <- 0
        iso_per_knownTr_cage_polya$total_iso=apply(iso_per_knownTr_cage_polya, 1 , function(X) as.integer(X[3]) + as.integer(X[5]) )

        iso_per_knownTr_cage_polya$FSM_cat = NA
        iso_per_knownTr_cage_polya$FSM_cat = apply(iso_per_knownTr_cage_polya, 1, function(X){
            if(as.numeric(X["FSM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["FSM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_cage_polya$FSM_bin = apply(iso_per_knownTr_cage_polya, 1, function(X){
            if(as.numeric(X["FSM_per_tr"]) < 8){
                return(as.character(X["FSM_per_tr"]))
            }else{
                return("8+")
            }})

        iso_per_knownTr_cage_polya$ISM_cat = NA
        iso_per_knownTr_cage_polya$ISM_cat = apply(iso_per_knownTr_cage_polya, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) == 1){
                return("Unique")
            }else if(as.numeric(X["ISM_per_tr"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_cage_polya$ISM_bin = apply(iso_per_knownTr_cage_polya, 1, function(X){
            if(as.numeric(X["ISM_per_tr"]) < 8){
                return(as.character(X["ISM_per_tr"]))
            }else{
                return("8+")
            }})


        iso_per_knownTr_cage_polya$total_cat = NA
        iso_per_knownTr_cage_polya$total_cat = apply(iso_per_knownTr_cage_polya, 1, function(X){
            if(as.numeric(X["total_iso"]) == 1){
                return("Unique")
            }else if(as.numeric(X["total_iso"]) > 1){
                return("Multiple")
            }else{
                return("NULL")
            }})

        iso_per_knownTr_cage_polya$total_bin = apply(iso_per_knownTr_cage_polya, 1, function(X){
            if(as.numeric(X["total_iso"]) < 8){
                return(as.character(X["total_iso"]))
            }else{
                return("8+")
            }})

        iso_per_knownTr_cage_polya$FSM_cat=factor(iso_per_knownTr_cage_polya$FSM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr_cage_polya$ISM_cat=factor(iso_per_knownTr_cage_polya$ISM_cat, levels=c("Unique", "Multiple"))
        iso_per_knownTr_cage_polya$total_cat=factor(iso_per_knownTr_cage_polya$total_cat, levels=c("Unique", "Multiple"))

        max_y=max(table(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$FSM_cat!="NULL"),"FSM_cat"])) + 10
        p.list[[i]] <- ggplot(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$FSM_cat!="NULL"),]) +
            geom_bar(aes(x=FSM_cat, fill=FSM_bin), color = "black", width = 0.5) +
            mytheme+
            geom_text(aes(x = FSM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
            scale_y_continuous(breaks = pretty_breaks(6), limits = c(0, max_y)) +
            scale_fill_brewer("Total FSM \nper reference ID", palette = "Blues") +
            labs(x="FSM per reference transcript",
                 y="Count of reference transcripts",
                 title=sample.names[[i]],
                 subtitle="Only FSM with CAGE support and polyA motif")

        max_y=max(table(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$ISM_cat!="NULL"),"ISM_cat"])) + 10
        p.list2[[i]] <- ggplot(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$ISM_cat!="NULL"),]) +
            geom_bar(aes(x=ISM_cat, fill=ISM_bin), color = "black", width = 0.5) +
            mytheme+
            geom_text(aes(x = ISM_cat, label = stat(count)), stat = "count", vjust = -0.5) +
            scale_y_continuous(breaks = pretty_breaks(6), limits = c(0, max_y)) +
            scale_fill_brewer("Total ISM \nper reference ID", palette = "Oranges") +
            labs(x="ISM per reference transcript",
                 y="Count of reference transcripts",
                 title=sample.names[[i]],
                 subtitle="Only ISM with CAGE support and polyA motif")

        max_y=max(table(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$total_cat!="NULL"),"total_cat"])) + 10
        p.list3[[i]] <- ggplot(iso_per_knownTr_cage_polya[which(iso_per_knownTr_cage_polya$total_cat!="NULL"),]) +
            geom_bar(aes(x=total_cat, fill=total_bin), color = "black", width = 0.5) +
            mytheme+
            geom_text(aes(x = total_cat, label = stat(count)), stat = "count", vjust = -0.5) +
            scale_y_continuous(breaks = pretty_breaks(6), limits = c(0, max_y)) +
            scale_fill_brewer("Total FSM+ISM \nper reference", palette = "Greens") +
            labs(x="FSM+ISM per reference transcript",
                 y="Count of reference transcripts",
                 title=sample.names[[i]],
                 subtitle="FSM+ISM with CAGE support and polyA motif")
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Reference Transcript Redundancy", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
    organize_in_grid_with_title(page_title, p.list2)
    organize_in_grid_with_title(page_title, p.list3)
    rm(iso_per_knownTr_cage)
}
rm(p.list)
rm(p.list2)
rm(p.list3)




s <- textGrob("Intra-Priming Quality Check", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)
# PLOT p30,p31,p32: percA by subcategory
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.FSMISM.list[[i]], aes(y=perc_A_downstream_TTS, x=subcategory, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.1, position = "dodge") +
                   mytheme +
                   xlab("Structural category") +
                   ylab("'A's, %") +
                   labs(title = paste0(sample.names[i], "\n\n"),
                        subtitle="Percent of genomic 'A's in downstream 20 bp\n\n") +
                   theme(legend.position="right", legend.title=element_blank()) +
                   theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                         axis.ticks.x = element_blank()) +
                   theme(strip.background = element_rect(color = "white"),
                         strip.placement = "outside", strip.text = element_text(size = 10)) +
                   scale_fill_manual(values=subcat.palette, drop=T) +
                   facet_grid(~structural_category, scales = "free", switch = "x")
}
page_title = textGrob("Possible Intra-Priming by Structural Category", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.NICNNC.list[[i]], aes(y=perc_A_downstream_TTS, x=subcategory, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.1, position = "dodge") +
                   mytheme +
                   xlab("Structural category") +
                   ylab("'A's, %") +
                   labs(title = paste0(sample.names[i], "\n\n"),
                        subtitle="Percent of genomic 'A's in downstream 20 bp\n\n") +
                   theme(legend.position="right", legend.title=element_blank()) +
                   theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                         axis.ticks.x = element_blank()) +
                   theme(strip.background = element_rect(color = "white"),
                         strip.placement = "outside", strip.text = element_text(size = 10)) +
                   scale_fill_manual(values=subcat.palette, drop=T) +
                   facet_grid(~structural_category, scales = "free", switch = "x")
}
page_title = textGrob("Possible Intra-Priming by Structural Category", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.other.list[[i]], aes(y=perc_A_downstream_TTS, x=subcategory, fill=subcategory)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.1, position = "dodge") +
                   mytheme +
                   xlab("Structural category") +
                   ylab("'A's, %") +
                   labs(title = paste0(sample.names[i], "\n\n"),
                        subtitle="Percent of genomic 'A's in downstream 20 bp\n\n") +
                   theme(legend.position="right", legend.title=element_blank()) +
                   theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
                         axis.ticks.x = element_blank()) +
                   theme(strip.background = element_rect(color = "white"),
                         strip.placement = "outside", strip.text = element_text(size = 10)) +
                   scale_fill_manual(values=subcat.palette, drop=T) +
                   facet_grid(~structural_category, scales = "free", switch = "x")
                   theme(axis.text.x = element_text(angle = 60, hjust = 1))
}
page_title = textGrob("Possible Intra-Priming by Structural Category", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(y=perc_A_downstream_TTS, x=structural_category, fill=exonCat)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) +
                   mytheme +
                   scale_fill_manual(breaks=c("Mono-Exon", "Multi-Exon"),
                                     labels=c("Mono-Exon Isoforms", "Multi-Exon Isoforms"), values=myPalette) +
                   ylab("'A's, %") +
                   theme(legend.position="bottom", legend.title=element_blank()) +
                   theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   labs(title = paste0(sample.names[i], "\n\n"),
                        subtitle = "Percent of genomic 'A's in downstream 20 bp\n\n") +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Mono- vs Multi-Exon Possible Intra-Priming", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)


p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    p.list[[i]] <- ggplot(data=data.class.list[[i]], aes(y=perc_A_downstream_TTS, x=structural_category, fill=coding)) +
                   geom_boxplot(color="black", size=0.3, outlier.size = 0.2) + mytheme +
                   scale_fill_manual(breaks=c("Coding", "Non coding"),
                                     labels=c("Coding Isoforms", "Non-Coding Isoforms"), values=myPalette[3:4]) +
                   ylab("'A's, % ") +
                   theme(legend.position="bottom", legend.title=element_blank() ) +
                   theme(axis.text.x = element_text(angle = 45)) +
                   theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12))+
                   labs(title = paste0(sample.names[i], "\n\n"),
                        subtitle = "Percent of genomic 'A's in downstream 20 bp\n\n") +
                   theme(axis.title.x=element_blank())
}
page_title = textGrob("Coding vs Non-Coding Possible Intra-Priming", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)






p.list = vector("list", length(class.files))
p.list2 = vector("list", length(class.files))
p.list3 = vector("list", length(class.files))
p.list4 = vector("list", length(class.files))
p.list5 = vector("list", length(class.files))
p.list6 = vector("list", length(class.files))
p.list7 = vector("list", length(class.files))
p.list8 = vector("list", length(class.files))
#### Rarefraction plots
if (saturation.curves == 'True') {
    for (i in seq_along(class.files)) {
        if (!all(is.na(data.class.list[[i]]$FL))) {
            FL.counts <- as.matrix(data.class.list[[i]]$FL)
            rownames(FL.counts) <- data.class.list[[i]]$isoform
            colnames(FL.counts) <- "FL"
            myfactors <- data.frame(sample = c("FL"))
            rownames(myfactors) = colnames(FL.counts)
            mybiotype = as.matrix(data.class.list[[i]]$coding)
            rownames(mybiotype) = data.class.list[[i]]$isoform
            mycategory = as.matrix(data.class.list[[i]]$structural_category)
            rownames(mycategory) = data.class.list[[i]]$isoform
            mydata = readData(data = FL.counts, factors = myfactors, biotype = mybiotype, category=mycategory)

            rarefact <- LR.rarefaction(mydata , samples = 1)
            res <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 1, depth.increase = 2, break.by = "category"))
            p.list[[i]] = res[[1]] + labs(title = paste0(sample.names[i], "\n\n"))
            p.list2[[i]] = res[[2]] + labs(title = paste0(sample.names[i], "\n\n"))
            res <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 2, depth.increase = 2, break.by = "category"))
            p.list3[[i]] = res[[1]] + labs(title = paste0(sample.names[i], "\n\n"))
            p.list4[[i]] = res[[2]] + labs(title = paste0(sample.names[i], "\n\n"))
            res <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 3, depth.increase = 2, break.by = "category"))
            p.list5[[i]] = res[[1]] + labs(title = paste0(sample.names[i], "\n\n"))
            p.list6[[i]] = res[[2]] + labs(title = paste0(sample.names[i], "\n\n"))
            res <- suppressWarnings(plot.rarefaction(rarefact, sample = 1, k = 5, depth.increase = 2, break.by = "category"))
            p.list7[[i]] = res[[1]] + labs(title = paste0(sample.names[i], "\n\n"))
            p.list8[[i]] = res[[2]] + labs(title = paste0(sample.names[i], "\n\n"))
        }
    }
    if (!all(sapply(p.list, is.null))) {
        page_title = textGrob("Saturation plot per structural category", gp=gpar(fontface="italic", fontsize=17))
        page_title2 = textGrob("Increments of detected isoforms per structural category", gp=gpar(fontface="italic", fontsize=17))
        organize_in_grid_with_title(page_title, p.list)
        organize_in_grid_with_title(page_title2, p.list2)
        organize_in_grid_with_title(page_title, p.list3)
        organize_in_grid_with_title(page_title2, p.list4)
        organize_in_grid_with_title(page_title, p.list5)
        organize_in_grid_with_title(page_title2, p.list6)
        organize_in_grid_with_title(page_title, p.list7)
        organize_in_grid_with_title(page_title2, p.list8)
        rm(res)
    }
    rm(p.list)
    rm(p.list2)
    rm(p.list3)
    rm(p.list4)
    rm(p.list5)
    rm(p.list6)
    rm(p.list7)
    rm(p.list8)
}



    
s <- textGrob("Features of Bad Quality", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)


x.list = vector("list", length(class.files))
x.min_cov.list = vector("logical", length(class.files))
x.predicted_NMD.list = vector("logical", length(class.files))
t2.RTS.list = vector("list", length(class.files))
t3.NMD.list = vector("list", length(class.files))
t3.Cov.list = vector("list", length(class.files))
t3.Cage.list = vector("list", length(class.files))
t3.PolyA.list = vector("list", length(class.files))

# t3.list = vector("list", length(class.files))
t3.a.Cov.list = vector("list", length(class.files))
t3.RTS.list = vector("list", length(class.files))
t3.annot.list = vector("list", length(class.files))
t3.SJ.list = vector("list", length(class.files))
t3.a.SJ.list = vector("list", length(class.files))
n_t3.SJ.list = vector("integer", length(class.files))
n_t3.RTS.list = vector("integer", length(class.files))

t3.list = vector("list", length(class.files))
t3.data.sets.list = vector("list", length(class.files))

p.list = vector("list", length(class.files))
### Bad quality control attributes
# setup a bunch of intermediate results reused until the end

#p28.RTS
for (i in seq_along(class.files)) {
    if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
        t3.data.sets.list[[i]] <- list()
        t3.list[[i]] <- list()
        # (Fran) ToDo: USE COVERAGE DATA LATER
        # for FSM, ISM, NIC, and NNC, plot the percentage of RTS and non-canonical junction
        x.list[[i]] <- filter(data.class.list[[i]], structural_category %in% c("FSM", "ISM", "NIC", "NNC" ) & exons > 1)
        
        t1.RTS <- group_by(x.list[[i]], structural_category, RTS_stage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
        t2.RTS.list[[i]] <- group_by(x.list[[i]], structural_category) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
        t3.RTS.list[[i]] <- merge(t1.RTS, t2.RTS.list[[i]], by="structural_category")
        rm(t1.RTS)
        t3.RTS.list[[i]] <- t3.RTS.list[[i]][-which(t3.RTS.list[[i]]$structural_category=="ISM"),]
        t3.RTS.list[[i]]$perc <- t3.RTS.list[[i]]$count.x / t3.RTS.list[[i]]$count.y * 100
        t3.RTS.list[[i]] <- subset(t3.RTS.list[[i]], RTS_stage=='TRUE');
        n_t3.RTS.list[[i]] <- dim(t3.RTS.list[[i]])[1];
        if (n_t3.RTS.list[[i]] > 0) {
            t3.RTS.list[[i]]$Var <- "RT switching"
        }

        p.list[[i]] <- ggplot(t3.RTS.list[[i]], aes(x=structural_category, y=perc)) +
                       geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[11], color="black") +
                       geom_text(label=paste(round(t3.RTS.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
                       scale_fill_manual(values = myPalette[9:11]) +
                       scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                       ylab("Isoforms, %") +
                       xlab("") +
                       mytheme +
                       theme(legend.position="bottom", axis.title.x = element_blank()) +
                       ggtitle(paste0(sample.names[[i]], "\n\n"))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("RT-switching", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


#p28.SJ
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
        # Liz: this is a placeholder for dealing with all_canonical being NA instead of "Non-canonical"
        x.list[[i]][is.na(x.list[[i]]$all_canonical), "all_canonical"] <- "Non-canonical"
        t1.SJ <- group_by(x.list[[i]], structural_category, all_canonical) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
        t3.SJ.list[[i]] <- merge(t1.SJ, t2.RTS.list[[i]], by="structural_category")
        rm(t1.SJ)
        t3.SJ.list[[i]]$perc <- t3.SJ.list[[i]]$count.x / t3.SJ.list[[i]]$count.y * 100
        t3.a.SJ.list[[i]] <- subset(t3.SJ.list[[i]], all_canonical=='Canonical');
        t3.SJ.list[[i]] <- subset(t3.SJ.list[[i]], all_canonical=='Non-canonical');
        n_t3.SJ.list[[i]] <- dim(t3.SJ.list[[i]])[1];
        if (n_t3.SJ.list[[i]] > 0) {
            t3.SJ.list[[i]]$Var <- "Non-canonical"
            t3.a.SJ.list[[i]]$Var <- 'Canonical'
        }

        p.list[[i]] <- ggplot(t3.SJ.list[[i]], aes(x=structural_category, y=perc)) +
                       geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[9] ,color="black") +
                       geom_text(label=paste(round(t3.SJ.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
                       scale_fill_manual(values = myPalette[9:11]) +
                       scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                       ylab("Isoforms, %") +
                       xlab("") +
                       mytheme +
                       theme(legend.position="bottom", axis.title.x = element_blank()) +
                       ggtitle(paste0(sample.names[[i]], "\n\n"))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Non-Canonical Junctions", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


# some more variable setup
for (i in seq_along(class.files)) {
    if (!all(is.na(x.list[[i]]$within_CAGE_peak))) {
        x.list[[i]][which(!x.list[[i]]$within_CAGE_peak),"Coverage_Cage"] <- "No Coverage CAGE"
        x.list[[i]][which(x.list[[i]]$within_CAGE_peak),"Coverage_Cage"] <- "Has Coverage CAGE"
        t1.Cage <- group_by(x.list[[i]], structural_category, Coverage_Cage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
        t3.Cage.list[[i]] <- merge(t1.Cage, t2.RTS.list[[i]], by="structural_category")
        rm(t1.Cage)
        t3.Cage.list[[i]]$perc <- t3.Cage.list[[i]]$count.x / t3.Cage.list[[i]]$count.y * 100
        t3.Cage.list[[i]] <- subset(t3.Cage.list[[i]], Coverage_Cage=='Has Coverage CAGE');
        t3.Cage.list[[i]]$Var <- t3.Cage.list[[i]]$Coverage_Cage
        t3.data.sets.list[[i]][[length(t3.data.sets.list[[i]]) + 1]] <- !all(is.na(data.class.list[[i]]$dist_to_cage_peak))
        t3.list[[i]][[length(t3.list[[i]]) + 1]] <- t3.Cage.list[[i]]
    }

    if (!all(is.na(data.class.list[[i]]$polyA_motif))) {
        x.list[[i]][which(is.na(x.list[[i]]$polyA_motif)),"Coverage_PolyA"] <- "No Coverage PolyA"
        x.list[[i]][which(!is.na(x.list[[i]]$polyA_motif)),"Coverage_PolyA"] <- "Has Coverage PolyA"
        t1.PolyA <- group_by(x.list[[i]], structural_category, Coverage_PolyA) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
        t3.PolyA.list[[i]] <- merge(t1.PolyA, t2.RTS.list[[i]], by="structural_category")
        rm(t1.PolyA)
        t3.PolyA.list[[i]]$perc <- t3.PolyA.list[[i]]$count.x / t3.PolyA.list[[i]]$count.y * 100
        t3.PolyA.list[[i]] <- subset(t3.PolyA.list[[i]], Coverage_PolyA=='Has Coverage PolyA');
        t3.PolyA.list[[i]]$Var <- t3.PolyA.list[[i]]$Coverage_PolyA
        t3.data.sets.list[[i]][[length(t3.data.sets.list[[i]]) + 1]] = !all(is.na(data.class.list[[i]]$polyA_motif))
        t3.list[[i]][[length(t3.list[[i]]) + 1]] = t3.PolyA.list[[i]]

    }
    
    x.list[[i]][which(x.list[[i]]$diff_to_gene_TSS<=50),"Annotation"] <- "Annotated"
    x.list[[i]][which(x.list[[i]]$diff_to_gene_TSS>50),"Annotation"] <- "Not annotated"
    t1.annot <- group_by(x.list[[i]], structural_category, Annotation) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')

    x.min_cov.list[[i]] = all(is.na(x.list[[i]]$min_cov))
    x.predicted_NMD.list[[i]] = all(is.na(x.list[[i]]$predicted_NMD))

    t3.annot.list[[i]] <- merge(t1.annot, t2.RTS.list[[i]], by="structural_category")
    rm(t1.annot)
    t3.annot.list[[i]]$perc <- t3.annot.list[[i]]$count.x / t3.annot.list[[i]]$count.y * 100
    t3.annot.list[[i]] <- subset(t3.annot.list[[i]], Annotation=='Annotated');
    t3.annot.list[[i]]$Var = t3.annot.list[[i]]$Annotation
}




# p28.Cov
if (sum(n_t3.SJ.list) > 0 & sum(n_t3.RTS.list) > 0 & !all(is.na(data.class.list[[1]]$min_cov))) {
    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
            if (!all(is.na(x.list[[i]]$min_cov))){
                x.list[[i]][which(x.list[[i]]$min_cov==0),"Coverage_SJ"]="Not Coverage SJ"
                x.list[[i]][which(x.list[[i]]$min_cov>0),"Coverage_SJ"]="Coverage SJ"
                t1.Cov <- group_by(x.list[[i]], structural_category, Coverage_SJ) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
                t3.Cov.list[[i]] <- merge(t1.Cov, t2.RTS.list[[i]], by="structural_category")
                rm(t1.Cov)
                t3.Cov.list[[i]]$perc <- t3.Cov.list[[i]]$count.x / t3.Cov.list[[i]]$count.y * 100
                t3.a.Cov.list[[i]] <- subset(t3.Cov.list[[i]], Coverage_SJ=='Coverage SJ');
                t3.Cov.list[[i]] <- subset(t3.Cov.list[[i]], Coverage_SJ=='Not Coverage SJ');
                t3.Cov.list[[i]]$Var = t3.Cov.list[[i]]$Coverage_SJ
                t3.a.Cov.list[[i]]$Var = t3.a.Cov.list[[i]]$Coverage_SJ
                t3.data.sets.list[[i]][[length(t3.data.sets.list[[i]]) + 1]] = !all(is.na(x.list[[i]]$min_cov))
                t3.list[[i]][[length(t3.list[[i]]) + 1]] = t3.a.Cov.list[[i]]
            }

            if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0 & !x.min_cov.list[[i]] & x.predicted_NMD.list[[i]]) {
                p.list[[i]] <- ggplot(t3.Cov.list[[i]], aes(x=structural_category, y=perc)) +
                               geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
                               geom_text(label=paste(round(t3.Cov.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
                               scale_fill_manual(values = myPalette[9:11]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Isoforms, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n"))
            } else if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0) {
                p.list[[i]] <- ggplot(t3.Cov.list[[i]], aes(x=structural_category, y=perc)) +
                               geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
                               geom_text(label=paste(round(t3.Cov.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
                               scale_fill_manual(values = myPalette[9:11]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Isoforms, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n"))

            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        page_title = textGrob("Splice Junctions Without Short Read Coverage", gp=gpar(fontface="italic", fontsize=17))
        organize_in_grid_with_title(page_title, p.list)
    }
    rm(p.list)
}



# p28.NMD
if (sum(n_t3.SJ.list) > 0 & sum(n_t3.RTS.list) > 0 & !all(is.na(data.class.list[[1]]$predicted_NMD))) {
    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
            if (!all(is.na(x.list[[i]]$predicted_NMD))){
                x.list[[i]][which(x.list[[i]]$predicted_NMD=="TRUE"),"predicted_NMD"] = "Predicted NMD"
                x.list[[i]][which(x.list[[i]]$predicted_NMD=="FALSE"),"predicted_NMD"] = "Not NMD predicted"
                t1.NMD <- group_by(x.list[[i]], structural_category, predicted_NMD) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
                t3.NMD.list[[i]] <- merge(t1.NMD, t2.RTS.list[[i]], by="structural_category")
                rm(t1.NMD)
                t3.NMD.list[[i]]$perc <- t3.NMD.list[[i]]$count.x / t3.NMD.list[[i]]$count.y * 100
                t3.NMD.list[[i]] <- subset(t3.NMD.list[[i]], predicted_NMD=='Predicted NMD');
                t3.NMD.list[[i]]$Var = t3.NMD.list[[i]]$predicted_NMD
            }

            if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0 & x.min_cov.list[[i]] & !x.predicted_NMD.list[[i]]) {
                p.list[[i]] <- ggplot(t3.NMD.list[[i]], aes(x=structural_category, y=perc)) +
                               geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[5], color="black") +
                               geom_text(label=paste(round(t3.NMD.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
                               scale_fill_manual(values = myPalette[9:11]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Isoforms, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n"))

            } else if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0) {
                p.list[[i]] <- ggplot(t3.NMD.list[[i]], aes(x=structural_category, y=perc)) +
                               geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[5], color="black") +
                               geom_text(label=paste(round(t3.NMD.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
                               scale_fill_manual(values = myPalette[9:11]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Isoforms, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n")) 
            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        page_title = textGrob("Nonsense-Mediated Decay by Structural Category", gp=gpar(fontface="italic", fontsize=17))
        organize_in_grid_with_title(page_title, p.list)
    }
    rm(p.list)
}


# p28
if (sum(n_t3.SJ.list) > 0 & sum(n_t3.RTS.list) > 0) {
    p.list = vector("list", length(class.files))
    page_title = textGrob("Quality Control Attributes Across Structural Categories", gp=gpar(fontface="italic", fontsize=17))
    for (i in seq_along(class.files)) {
        if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
           if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0 & !x.min_cov.list[[i]] & x.predicted_NMD.list[[i]]) {
                t3 <- rbind(t3.RTS.list[[i]][,c(1,5,6)], t3.SJ.list[[i]][,c(1,5,6)], t3.Cov.list[[i]][,c(1,5,6)])
                p.list[[i]] <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
                               geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
                               scale_fill_manual(values = myPalette[9:11]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Transcripts, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n")) +
                               theme(legend.title = element_blank())
                page_title = textGrob("Summary Features of Bad Quality", gp=gpar(fontface="italic", fontsize=17))
            } else if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0 & x.min_cov.list[[i]] & x.predicted_NMD.list[[i]]) {
                t3 = rbind(t3.RTS.list[[i]][,c(1,5,6)], t3.SJ.list[[i]][,c(1,5,6)])
                p.list[[i]] <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
                               geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
                               scale_fill_manual(values = myPalette[c(9,11)]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Transcripts, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n")) +
                               theme(legend.title = element_blank())
            } else if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0 & x.min_cov.list[[i]] & !x.predicted_NMD.list[[i]]) {
                t3 = rbind(t3.RTS.list[[i]][,c(1,5,6)], t3.SJ.list[[i]][,c(1,5,6)], t3.NMD.list[[i]][,c(1,5,6)])
                p.list[[i]] <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
                               geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
                               scale_fill_manual(values = myPalette[c(9,5,11)]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Transcripts, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n")) +
                               theme(legend.title = element_blank())
            } else if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0) {
                t3 = rbind(t3.RTS.list[[i]][,c(1,5,6)], t3.SJ.list[[i]][,c(1,5,6)], t3.Cov.list[[i]][,c(1,5,6)], t3.NMD.list[[i]][,c(1,5,6)])
                p.list[[i]] <- ggplot(data=t3, aes(x=structural_category, y=perc, fill= Var)) +
                               geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
                               scale_fill_manual(values = myPalette[c(9,10,5,11)]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Transcripts, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n")) +
                               theme(legend.title = element_blank()) 
            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        organize_in_grid_with_title(page_title, p.list)
    }
    rm(p.list)
}




s <- textGrob("Features of Good Quality", gp=gpar(fontface="italic", fontsize=17), vjust = 0)
grid.arrange(s)

# p28.a.Cage
if (!all(is.na(data.class.list[[1]]$dist_to_cage_peak))) {
    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
            p.list[[i]] <- ggplot(t3.Cage.list[[i]], aes(x=structural_category, y=perc)) +
                           geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[2] ,color="black") +
                           geom_text(label=paste(round(t3.Cage.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
                           scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                           ylab("Isoforms, %") +
                           xlab("") +
                           mytheme +
                           theme(legend.position="bottom", axis.title.x = element_blank()) +
                           ggtitle(paste0(sample.names[[i]], "\n\n")) +
                           theme(legend.title = element_blank())
        }
        if (!all(sapply(p.list, is.null))) {
            page_title = textGrob("CAGE Support", gp=gpar(fontface="italic", fontsize=17))
            organize_in_grid_with_title(page_title, p.list)
        }
        rm(p.list)
    }
}


# p28.a.annot
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
        p.list[[i]] <- ggplot(t3.annot.list[[i]], aes(x=structural_category, y=perc)) +
                       geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[6] ,color="black") +
                       geom_text(label=paste(round(t3.annot.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
                       scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                       ylab("Isoforms, %") +
                       xlab("") +
                       mytheme +
                       theme(legend.position="bottom", axis.title.x = element_blank()) +
                       ggtitle(paste0(sample.names[[i]], "\n\n")) 
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("Annotation Support", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)


# p28.a.PolyA
if (!all(is.na(data.class.list[[1]]$polyA_motif))) {
    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
            p.list[[i]] <- ggplot(t3.PolyA.list[[i]], aes(x=structural_category, y=perc)) +
                           geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[3] ,color="black") +
                           geom_text(label=paste(round(t3.PolyA.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
                           scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                           ylab("Isoforms, %") +
                           xlab("") +
                           mytheme +
                           theme(legend.position="bottom", axis.title.x = element_blank()) +
                           ggtitle(paste0(sample.names[[i]], "\n\n")) +
                           guides(fill = guide_legend(title = "QC Attributes") )
        }
    }
    if (!all(sapply(p.list, is.null))) {
        page_title = textGrob("PolyA Support", gp=gpar(fontface="italic", fontsize=17))
        organize_in_grid_with_title(page_title, p.list)
    }
    rm(p.list)
}


# p28.a.SJ
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
        p.list[[i]] <- ggplot(t3.a.SJ.list[[i]], aes(x=structural_category, y=perc)) +
                       geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[7] ,color="black") +
                       geom_text(label=paste(round(t3.a.SJ.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) + 
                       scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                       ylab("Isoforms, %") +
                       xlab("") +
                       mytheme +
                       theme(legend.position="bottom", axis.title.x = element_blank()) +
                       ggtitle(paste0(sample.names[[i]], "\n\n"))

}
page_title = textGrob("All Canonical Junctions", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)


# p28.a.Cov
if (!all(is.na(data.class.list[[i]]$min_cov))) {
    p.list = vector("list", length(class.files))
    for (i in seq_along(class.files)) {
        if (n.data.junction.list[[i]] > 0) { ### was nrow(data.junction)
            if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0 & !x.min_cov.list[[i]] & x.predicted_NMD.list[[i]]) {
                p.list[[i]] <- ggplot(t3.a.Cov.list[[i]], aes(x=structural_category, y=perc)) +
                               geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
                               geom_text(label=paste(round(t3.a.Cov.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
                               scale_fill_manual(values = myPalette[9:11]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Isoforms, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n")) 
            } else if (n_t3.SJ.list[[i]] > 0 & n_t3.RTS.list[[i]] > 0) {
                p.list[[i]] <- ggplot(t3.a.Cov.list[[i]], aes(x=structural_category, y=perc)) +
                               geom_col(position='dodge', width = 0.7,  size=0.3, fill=myPalette[10], color="black") +
                               geom_text(label=paste(round(t3.a.Cov.list[[i]]$perc, 1),"%",sep=''), position = position_dodge(0.9),vjust = -0.8) +
                               scale_fill_manual(values = myPalette[9:11]) +
                               scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                               ylab("Isoforms, %") +
                               xlab("") +
                               mytheme +
                               theme(legend.position="bottom", axis.title.x = element_blank()) +
                               ggtitle(paste0(sample.names[[i]], "\n\n")) 
            }
        }
    }
    if (!all(sapply(p.list, is.null))) {
        page_title = textGrob("Splice Junctions With Short Read Coverage", gp=gpar(fontface="italic", fontsize=17))
        organize_in_grid_with_title(page_title, p.list)
    }
    rm(p.list)
}


# p28.a
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    t3.aa <-  rbind(t3.annot.list[[i]][,c("structural_category", "perc", "Var")], t3.a.SJ.list[[i]][,c(1,5,6)])

    for(j in 1:length(t3.list[[i]])) {
        c = data.frame(t3.list[[i]][[j]])
        if (t3.data.sets.list[[i]][[j]]) {
            t.temp = t3.aa
            t3.aa = rbind(t.temp, c[,c(1,5,6)])
        }
    }

    p.list[[i]] <- ggplot(data=t3.aa, aes(x=structural_category, y=perc, fill= Var)) +
                   geom_bar(position = position_dodge(), stat="identity", width = 0.7,  size=0.3, color="black") +
                   guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
                   scale_fill_manual(values = c(myPalette)) +
                   scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                   ylab("Transcripts, %") +
                   xlab("") +
                   mytheme +
                   theme(legend.position="bottom", axis.title.x = element_blank()) +
                   ggtitle(paste0(sample.names[[i]], "\n\n")) +
                   theme(axis.text.y = element_text(size=10),
                         axis.text.x  = element_text(size=10))+
                   theme(legend.title = element_blank())
}
page_title = textGrob("Good Quality Control Attributes Across Structural Categories", gp=gpar(fontface="italic", fontsize=17))
organize_in_grid_with_title(page_title, p.list)
rm(p.list)





#TSS ratio
p.list = vector("list", length(class.files))
for (i in seq_along(class.files)) {
    data.ratio = rbind(data.refmatch.list[[i]][,c(6,47)], data.ISM.list[[i]][,c(6,47)])
    if (!all(is.na(data.ratio$ratio_TSS))) {
        require(scales)
        p.list[[i]] = ggplot(data.ratio, aes(x=ratio_TSS, fill=structural_category)) + 
                      geom_density(alpha=0.6)+
                      labs(x="TSS ratio, log2", y="Density", title=paste0(sample.names[[i]], "\n\n")) +
                      scale_fill_manual(values = myPalette, breaks=c("FSM", "ISM"),
                                        labels=c("FSM reference match", "ISM"), drop=F)+
                      geom_vline(xintercept=1, linetype="dashed", color = "red")+
                      scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
                      mytheme +
                      theme(legend.position="bottom", legend.title = element_blank())
        p.list[[i]] <- p.list[[i]] + scale_x_continuous(trans='log2', breaks = trans_breaks("log2", function(x) 2^x),
                                                        labels = trans_format("log2", math_format(2^.x)))
    }
}
if (!all(sapply(p.list, is.null))) {
    page_title = textGrob("TSS Ratio\nFSM Reference Match vs ISM", gp=gpar(fontface="italic", fontsize=17))
    organize_in_grid_with_title(page_title, p.list)
}
rm(p.list)

dev.off()

    ##### % of FSM or ISM associated to the same transcript ( histogram )







    # PLOT pn1.2: Splice Junction relative coverage (if coverage and expression provided)
    ##### NEEDS TRANSCRIPT_COORD VALUES IN JUNCTIONS FILE (???)

    #if (nrow(data.junction) > 0){
    #  if (!all(is.na(data.junction$total_coverage)) & !all(is.na(data.class$iso_exp))){

    #   data.junction$isoExp = data.class[data.junction$isoform, "iso_exp"]

    #    total = aggregate(cbind(total_coverage,isoExp,transcript_coord) ~ junctionLabel, data = data.junction,
    #                      FUN = function(x) c(mn = sum(x), n = min(x) ) )

    #    total$relCov = total$total_coverage[,"n"] / total$isoExp[,"mn"]
    #    total$minTSS = total$transcript_coord[,"n"]

    #    uniqJunc = unique(data.junction[,c("junctionLabel", "canonical_known", "total_coverage")])
    #    uniqJunc$notCov = uniqJunc$total_coverage == 0

    #    uniqueJunc_nonCov = as.data.frame(table(uniqJunc[uniqJunc$totalCoverage==0,"canonical_known"])/table(uniqJunc$canonical_known)*100)

    #    uniqJunc2 = merge(total, uniqJunc, by=1)
    #    uniqJunc2$TSSrange =cut(uniqJunc2$minTSS, breaks = c(0,40,80,120,160,200,10000000), labels = c("0-40", "41-80", "81-120", "121-160", "161-200",">200"))



            # calculate total expression associated to each unique junction
    #    sumExpPerJunc = tapply(data.junction$isoExp, data.junction$junctionLabel, sum)

    #    data.junction$sumIsoExp = sumExpPerJunc[data.junction$junctionLabel]

    #    data.junction$relCov = data.junction$total_coverage / data.junction$sumIsoExp

    #    max_dist = max(data.junction$transcript_coord) +1

    #    data.junction$TSSrange = cut(data.junction$transcript_coord, breaks = c(0,20,40,60,80,100,120,140,160,180,200,max_dist), labels = c("0-20", "21-40","41-80","61-80", "81-100","101-120", "121-140","141-160", "161-180", "181-200", ">200"))

    #    pn1.2 <-ggplot(data=data.junction[data.junction$relCov<1,], aes(y=relCov,x=TSSrange,fill=canonical_known)) +
    #      geom_boxplot(outlier.size = 0.2, size=0.3) +
    #      scale_fill_manual(values = myPalette[c(1,7,3,2)], drop=FALSE) +
    #      ylab("Relative coverage") +
    #      xlab("# TSS distance range") +
    #      mytheme_bw +
    #      theme(legend.position="bottom", legend.title=element_blank())  +
    #      ggtitle( "Junctions Relative Coverage\n\n\n") +
    #      theme(axis.text.x = element_text(angle = 45,margin=margin(15,0,0,0), size=12))


    #  }else{    uniqJunc = unique(data.junction[,c("junctionLabel", "canonical_known")])
    #  }
    #}





    #
    #df.polyA_freq.list[[i]] = df.polyA_freq
    #df.polyA_subcat.list[[i]] = df.polyA_subcat
    #df.cage.list[[i]] = df.cage
    #df.cage_subc.list[[i]] = df.cage_subc


###** Output plots

# if (report.format == 'both') {
#     invisible(generatePDFreport())
#     rmarkdown::render(input = paste(utilities.path, "/report_qc/SQANTI3_report.Rmd", sep = "/"), 
#                       intermediates_dir = output_directory, 
#                       output_dir =  output_directory, 
#                       output_file=html.report.file)
# } else if (args[6] == 'pdf' & args[6] != 'html'){
#    invisible(generatePDFreport())
####generatePDFreport()
# } else {
#     rmarkdown::render(input = paste(utilities.path, "/report_qc/SQANTI3_report.Rmd", sep = "/"), 
#                       intermediates_dir = output_directory,
#                       output_dir =  output_directory, 
#                       output_file=html.report.file)
# }

