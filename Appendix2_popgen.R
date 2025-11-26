library(diveRsity)
library(adegenet)
library(dartR)
library(LEA)
library(ggplot2)
library(ggrepel)

snp2gen(infile="infiles_genepop/elep_bypop_snps.txt", prefix_length = 5) # This is just a text file of SNPs. They are formatted with one individual per column, and one snp per row (with SNPs in "nucleotide format"). NOTE! This gives the file a weird name, and there isnʻt a built-in way to change it. So you have to manually change the file name for the next step to work. 
	
snps <- read.genepop(file="infiles_genepop/elep_bypop.gen") #reads in genepop file, and convers to genind format needed for adegenet. 

######################
## dapc in adegenet ##
######################

grp <- find.clusters(snps, max.n.clust = 10, criterion = "diffNgroup", perc.pca = 75)
dapc1 <- dapc(snps, grp$grp)

pops <- read.table("pops.txt", sep="\t", header=T)

x <- dapc1$ind.coord[,1]
y <- dapc1$ind.coord[,2]
df <- data.frame(x = dapc1$ind.coord[,1], y = dapc1$ind.coord[,2], ind=pops$ind_id, group= pops$group)

ggplot(df, aes(x=x, y=y, color=group)) +
    geom_point()+
    geom_label_repel(aes(label=ind), max.overlaps = 30)+ 
    theme_bw()+
    guides(color=guide_legend(title="Islands"))+ 
    geom_vline(xintercept=0) + geom_hline(yintercept=0)

#compoplot(dapc1, legend=F, show.lab=T, lab=pops$ind_id) #This makes a "structure-like" bar plot with the population assignments according to the above PCA. We are doing this independently in LEA below, so not plotting the dapc barchart.

##########
## Fst ###
##########

genlight_snps <- gi2gl(snps) #For the Reich Fst calcualtion to work, we need to get the genid file (needed for adegenet) into genlight format, this script does that. 

source("reich_fst.R")
fst <- reich.fst(genlight_snps, bootstrap=100, plot = T)

#write.csv(file="reich_fst_pairwise_bypop.csv", fst[[1]])

###################
## fixed alleles ##
###################

fixed <- gl.fixed.diff(genlight_snps)
write.csv(file="fixed_diffs_bypop.csv", fixed$fd)
    
#various other stats 
gl.basic.stats(genlight_snps, digits=2)
gl.report.heterozygosity(genlight_snps)
private <- gl.report.pa(genlight_snps)
private_df <- as_tibble(private) %>% 
    select("pop1", "pop2", "fixed", "totalpriv")

#################
## snmf in LEA ##
#################

source("genid2structure.R")
genind2structure(snps, file="bypop_structure.txt", pops=F) # We ultimately need the SNPs in "geno" format for using LEA, but they donʻt have a built-in way to do that for genlight or genid files. So first we convert to structure format using this code. 
struct2geno("bypop_structure.txt", extra.row=1, extra.column =1, ploidy=2, FORMAT=2) # Now read the structure-formatted file and convert to geno format.

input.file = "snmf/bypop_structure.txt.geno"
input.file <- read.geno(input.file)

project = NULL
project = snmf(input.file, K = 1:10, entropy = T, repetitions=10, project="new")

plot(project, col = "blue", pch = 19, cex = 1.2)

best = which.min(cross.entropy(project, K = 2)) #This picks the best model for a given K based on cross-entropy score. Need to change the "K=" to whatever you want to plot, and need to give more colors depending on number of Ks. 

my.colors <- c("tomato", "lightblue",
               "olivedrab", "gold")
               
barchart(project, K = 2, run = best,
         border = NA, space = 0,
         col = my.colors,
         xlab = "Individuals",
         ylab = "Ancestry proportions",
         main = "Ancestry matrix") -> bp
         
axis(1, at = 1:length(bp$order),
     labels = bp$order, las=1,
     cex.axis = .4) # This labels the bars by number, corresponding to the order of the individuals in the input file (i.e., individual 15 could be in column 1 on the plot)
