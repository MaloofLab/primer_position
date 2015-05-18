---
  title: "gene structure primer position"
  output: html_document
---
  
###########
## How to draw gene structure graph and mapping primer position/direction?
##      ver0.1 (try and error)
################
```{r eval=FALSE} # run this section once
source("http://bioconductor.org/biocLite.R")
biocLite(ggbio)
# making Arabidopsis thaliana TxDb package
library("org.At.tair.db")
class(org.At.tair.db)
# columns(org.At.tair.db) # http://www.bioconductor.org/help/workflows/annotation/annotation/
# keytypes(org.At.tair.db)
# ids<-head(keys(org.At.tair.db, keytype="ENTREZID"))
# select(org.At.tair.db, keys=ids, columns="SYMBOL", keytype="ENTREZID")

 TxDb.At<-makeTxDbFromGFF("/Volumes/Data8/NGS_related/Arabidopsis_analysis/reference/TAIR10_GFF3_genes.gff",format="gff3",organism="Arabidopsis thaliana")
 makeTxDbPackage(txdb=TxDb.At,
                 version="1.0.0",
                maintainer="Kazunari Nozue <knozue@ucavis.edu>",
                author="Kazunari Nozue",
                destDir="/Volumes/Data6/bioconductor_R")
system("tar -zcvf TxDb.Athaliana.tar.gz TxDb.Athaliana") #in /Volumes/Data6/bioconductor_R
install.packages("/Volumes/Data6/bioconductor_R/TxDb.Athaliana.tar.gz", repos = NULL, type = "source")
```


library("TxDb.Athaliana", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")

library(ggbio)

# an example
transcripts(TxDb.Athaliana)
transcripts.Athaliana.selected2<-transcripts(TxDb.Athaliana,vals<-list(gene_id="AT1G01040"))
gr.Athaliana <- crunch(TxDb.Athaliana,which=transcripts.Athaliana.selected2)
## change column to 'model'
colnames(values(gr.Athaliana))[4] <- "model"
grl.Athaliana <- split(gr.Athaliana, gr.Athaliana$tx_name)
autoplot(grl.Athaliana)
autoplot(grl.Athaliana,aes(type=model)) # good!

# multiple genes
transcripts.Athaliana.selected2<-transcripts(TxDb.Athaliana,vals<-list(gene_id=c("AT1G01040","AT1G01050","AT1G01060")))
gr.Athaliana <- crunch(TxDb.Athaliana,which=transcripts.Athaliana.selected2)
## change column to 'model'
colnames(values(gr.Athaliana))[4] <- "model"
grl.Athaliana <- split(gr.Athaliana, gr.Athaliana$tx_name)
autoplot(grl.Athaliana,aes(type=model)) # good!
ggplot() + geom_alignment(grl.Athaliana,type=model)
#autoplot(grl.Athaliana) + geom_alignment(type="model") # ??
# facet by gene_id
#p <-autoplot(grl.Athaliana,aes(type=model),facets=.~strand)  # does no work

## use tracks
# subset by AGI name
transcripts.Athaliana.selected.AT1G01040<-transcripts(TxDb.Athaliana,vals<-list(gene_id="AT1G01040"))
transcripts.Athaliana.selected.AT1G01050<-transcripts(TxDb.Athaliana,vals<-list(gene_id="AT1G01050"))
transcripts.Athaliana.selected.AT1G01060<-transcripts(TxDb.Athaliana,vals<-list(gene_id="AT1G01060"))
# crunch
gr.Athaliana.AT1G01040 <- crunch(TxDb.Athaliana,which=transcripts.Athaliana.selected.AT1G01040)
gr.Athaliana.AT1G01050 <- crunch(TxDb.Athaliana,which=transcripts.Athaliana.selected.AT1G01050)
gr.Athaliana.AT1G01060 <- crunch(TxDb.Athaliana,which=transcripts.Athaliana.selected.AT1G01060)
## change column to 'model'
colnames(values(gr.Athaliana.AT1G01040))[4] <- "model"
grl.Athaliana.AT1G01040 <- split(gr.Athaliana.AT1G01040, gr.Athaliana.AT1G01040$tx_name)

colnames(values(gr.Athaliana.AT1G01050))[4] <- "model"
grl.Athaliana.AT1G01050 <- split(gr.Athaliana.AT1G01050, gr.Athaliana.AT1G01050$tx_name)

colnames(values(gr.Athaliana.AT1G01060))[4] <- "model"
# add primer
#gr.Athaliana.AT1G01060<-c(gr.Athaliana.AT1G01060,primers) # primier became "model" so this is not good idea
grl.Athaliana.AT1G01060 <- split(gr.Athaliana.AT1G01060, gr.Athaliana.AT1G01060$tx_name)

p1<-autoplot(grl.Athaliana.AT1G01040,aes(type=model))
p2<-autoplot(grl.Athaliana.AT1G01050,aes(type=model))
p3<-autoplot(grl.Athaliana.AT1G01060,aes(type=model))

# add primer info for AT1G01060 (fake primers)
primer1 <- GRanges("Chr1", IRanges(34000, 34000 + 30),strand="+",tx_id="",tx_name="oJM1000",gene_id="AT1G01060",model="primer")
primer2 <- GRanges("Chr1", IRanges(35000 -30, 35000),strand="-",tx_id="",tx_name="oJM1001",gene_id="AT1G01060",model="primer")
primers<-c(primer1,primer2)
primers

primers.df<-as.data.frame(ranges(primers))
rownames(primers.df)<-c("oJM1000","oJM1001")

p3.primer<-ggplot() + geom_text(data=primers.df,aes(x=primers.df[,"start"],y=rep(1.5,2),label=rownames(primers.df),size=0.5)) + scale_y_continuous(limits=c(0.5,1.7))
p3.primer<- p3.primer + geom_arrowrect(primers,aes(fill=strand,y=100),arrow.head=1,rect.height=0.4) 
q<-tracks(p3,p3.primer,heights = c(14, 1))
q<-q + theme(axis.text.y=element_blank())
q<-q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
             panel.background = element_blank(), axis.line.y = element_line(colour = "white"),
             axis.ticks.y=element_line(color="white"))
q

#tracks('AT1G01040'=p1,'AT1G01050'=p2,'AT1G01060'=p3)
#alignPlots('AT1G01040'=p1,'AT1G01050'=p2,'AT1G01060'=p3) # does not work
# better to plot in separate window (eg use grid for each gene)

# cf http://www.tengfei.name/ggbio/docs/man/autoplot-method.html


### making general funciton ######

primer.info<-data.frame(primername=c("oJM1000","oJM1001"),
                        AGI=c("AT1G01060","AT1G01060"),
                        Chr=c("Chr1","Chr1"),
                        start=c(34000,34970),
                        end=c(34030,35000),
                        strand=c("+","-"),
                        pairedprimer=c("oJM1001","oJM1000"),
                        gPCRsize=c("1000","1000"),
                        RTPCRsize=c("800","800")                        
                            )
primers<-primer.info[1:2,] #you can select primer of interest here
AGI<-unique(primer.info[1:2,"AGI"]) # this should be only one AGI

primer.gene.plot<-function(primers,AGI) {
  # making primer GRange object
  # add primer info for AT1G01060 (fake primers)
  primer1 <- GRanges(primers[1,"Chr"], IRanges(primers[1,"start"], primers[1,"end"]),strand=primers[1,"strand"],tx_id="",tx_name=primers[1,"primername"],gene_id=primers[1,"AGI"],model="primer")
  primer2 <- GRanges(primers[2,"Chr"], IRanges(primers[2,"start"], primers[2,"end"]),strand=primers[2,"strand"],tx_id="",tx_name=primers[2,"primername"],gene_id=primers[2,"AGI"],model="primer") 
  primers.temp<-c(primer1,primer2)
  print(primers.temp) 
  primers.df<-as.data.frame(ranges(primers.temp))
  rownames(primers.df)<-primers[1:2,"primername"]
  print(primers.df)
  print(rownames(primers.df))
  # extracting gene structure information
  transcripts.Athaliana.selected.temp<-transcripts(TxDb.Athaliana,vals<-list(gene_id=AGI))
  gr.Athaliana.temp <- crunch(TxDb.Athaliana,which=transcripts.Athaliana.selected.temp)
  ## change column to 'model'
  colnames(values(gr.Athaliana.temp))[4] <- "model"
  grl.Athaliana.temp <- split(gr.Athaliana.temp, gr.Athaliana.temp$tx_name)
  # plot gene structure
  p.temp<-autoplot(grl.Athaliana.temp,aes(type=model))
  # plot primer
  p.primer<-ggplot() + geom_text(data=primers.df,aes(x=primers.df[,"start"],y=rep(1.5,2),label=primers[1:2,"primername"],size=0.5)) + scale_y_continuous(limits=c(0.5,1.7))
  p.primer<- p.primer + geom_arrowrect(primers.temp,aes(fill=strand,y=100),arrow.head=1,rect.height=0.4)  
  # merge plots
  q<-tracks(p.temp,p.primer,heights = c(14, 1))
  q<-q + theme(axis.text.y=element_blank())
  q<-q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
               panel.background = element_blank(), axis.line.y = element_line(colour = "white"),
               axis.ticks.y=element_line(color="white"))
  return(q)
  
} # this function retunrs ggplot object

p<-primer.gene.plot(primers, AGI)
p



########################## trials below ####
library(rtracklayer)
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/gffMod.R") # Imports gffFeat function. It can also be called by old name: 'gffMod'. 
#system("wget ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff") # Or manually download this file to your current working directory. 
system("curl -O ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff") # Or manually download this file to your current working directory.  (050515)


gff <- import.gff("/Volumes/Data8/NGS_related/Arabidopsis_analysis/reference/TAIR10_GFF3_genes.gff", asRangedData=FALSE)
gffmod <- getFeat(x=gff, format="gff", range_types=c("intergenic", "gene_red")) # Argument 'range_types' supports the following values: intergenic, gene_red, intron, gene and exon. The input format ("gff" or "gtf") can be specified under the 'format' argument. Note: the computation of the intron ranges will take some time.
gffmod # Returns GRanges object.

### subsetting GRobject by ID (093012)
 gff[grep("ID=AT1G01010",elementMetadata(gff)$group),]
# GRanges with 3 ranges and 4 elementMetadata values:
      # seqnames       ranges strand |     type   source    phase
         # <Rle>    <IRanges>  <Rle> | <factor> <factor> <factor>
  # [1]     Chr1 [3631, 5899]      + |     gene   TAIR10     <NA>
  # [2]     Chr1 [3631, 5899]      + |     mRNA   TAIR10     <NA>
  # [3]     Chr1 [3760, 5630]      + |  protein   TAIR10     <NA>
                                                                 # group
                                                           # <character>
  # [1]             ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
  # [2]         ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
  # [3] ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
  # ---
  # seqlengths:
   # Chr1 Chr2 Chr3 Chr4 Chr5 ChrC ChrM
     # NA   NA   NA   NA   NA   NA   NA
##################
# ggbio package (050515)
source("http://bioconductor.org/biocLite.R")
biocLite("arabidopsis.db0")
biocLite("Homo.sapiens")
library(Homo.sapiens)
class(Homo.sapiens)

data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)


library(ggbio)
range.test<-range(gff[grep("ID=AT1G01010",elementMetadata(gff)$group),])
p<-autoplot(gff[grep("ID=AT1G01010",elementMetadata(gff)$group),],which=range.test, label.color = "black", color = "brown",fill = "white")
p
p<-autoplot(gff,which=range.test)
p # range.test did not work 
# example of GRL
library(biovizBase)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr.txdb <- crunch(txdb, which = wh)
## change column to 'model'
colnames(values(gr.txdb))[4] <- "model"
grl <- split(gr.txdb, gr.txdb$tx_id)
## fake some randome names
names(grl) <- sample(LETTERS, size = length(grl), replace = TRUE)
grl
autoplot(grl,aes(type=model))
autoplot(grl)
ggplot() +geom_alignment(grl,type="model")




##

#
range.test2<-range(gff[grep("AT1G01010",elementMetadata(gff)$group),])
p<-autoplot(gff[grep("AT1G01010",elementMetadata(gff)$group),],which=range.test2)
p # range.test2 did not work 

colnames(values(gff))
colnames(values(gff))<-c("source","model","score","phase","group")
range.test2<-range(gff[grep("AT1G01030",elementMetadata(gff)$group),])

gff.s<-gff[c(grep("AT1G01030",elementMetadata(gff)$group),grep("AT1G01040",elementMetadata(gff)$group)),]
p<-autoplot(gff.s,aes(type=model))
p # range.test2 did not work 
ggplot() + geom_alignment(gff.s,type="model")

# 
library("org.At.tair.db")
class(org.At.tair.db)
# [1] "OrgDb"
# attr(,"package")
# [1] "AnnotationDbi"
columns(org.At.tair.db) # http://www.bioconductor.org/help/workflows/annotation/annotation/
keytypes(org.At.tair.db)
ids<-head(keys(org.At.tair.db, keytype="ENTREZID"))
select(org.At.tair.db, keys=ids, columns="SYMBOL", keytype="ENTREZID")
# Making an OrganismDb package
library(OrganismDbi)
gd <- list(join1 = c(GO.db="GOID", org.Hs.eg.db="GO"),
           join2 = c(org.Hs.eg.db="ENTREZID",
           TxDb.Hsapiens.UCSC.hg19.knownGene="GENEID"))


# http://www.bioconductor.org/help/workflows/annotation/annotation/#sample-workflow-TxDb
TxDb.At<-makeTxDbFromGFF("/Volumes/Data8/NGS_related/Arabidopsis_analysis/reference/TAIR10_GFF3_genes.gff",format="gff3",organism="Arabidopsis thaliana")
makeTxDbPackage(txdb=TxDb.At,
                version="1.0.0",
                maintainer="Kazunari Nozue <knozue@ucavis.edu>",
                author="Kazunari Nozue",
                destDir="/Volumes/Data6/bioconductor_R")
system("tar -zcvf TxDb.Athaliana.tar.gz TxDb.Athaliana") #in /Volumes/Data6/bioconductor_R
install.packages("/Volumes/Data6/bioconductor_R/TxDb.Athaliana.tar.gz", repos = NULL, type = "source")
library("TxDb.Athaliana", lib.loc="/Library/Frameworks/R.framework/Versions/3.2/Resources/library")
gd <- list(join1 = c(GO.db="GOID", org.At.tair.db="GO"),
           join2 = c(org.At.tair.db="ENTREZID",TxDb.Athaliana="GENEID"))

makeOrganismPackage(pkgname = "Arabidopsis.thaliana",
                    graphData = gd,
                    organism = "Arabidopsis thaliana",
                    version = "1.0.0",
                    maintainer = "Kazunari Nozue <knozue@ucdavis.edu>",
                    author = "Kazunari Nozue",
                    destDir = "/Volumes/Data6/bioconductor_R",
                    license = "Artistic-2.0")


system("tar -zcvf Arabidopsis.thaliana.tar.gz Arabidopsis.thaliana")
install.packages("/Volumes/Data6/bioconductor_R/Arabidopsis.thaliana.tar.gz", repos = NULL, type = "source")
library(Arabidopsis.thaliana)
keytypes(Arabidopsis.thaliana)
p<-autoplot(Arabidopsis.thaliana)
p
q<-autoplot(Homo.sapiens,which=GRanges("chr2", IRanges(1e6, 1e6+100000)))
q # OK
p<-autoplot(Arabidopsis.thaliana,which=GRanges("Chr2", IRanges(1e6, 1e6+100000)))
p # error
p2<-autoplot(TxDb.Athaliana,which=GRanges("Chr2", IRanges(1e6, 1e6+1000)))
p2 # works! 

# how to see CDS, exons? and strand info
p2<-autoplot(TxDb.Athaliana,which=GRanges("Chr2", IRanges(1e6, 1e6+5000)),label.color="red",aes(color=TXSTRAND))
p2<-autoplot(TxDb.Athaliana,which=GRanges("Chr2", IRanges(1e6, 1e6+5000)),label.color="red")
p2



# how to select region of interest?
select(TxDb.Athaliana,keys="AT1G01040.1",columns="TXNAME",keytype="TXNAME")
select(TxDb.Athaliana,keys="AT1G01040.1",columns=columns(TxDb.Athaliana),keytype="TXNAME")
#p2<-autoplot(TxDb.Athaliana,which=GRange(select(TxDb.Athaliana,keys="AT1G01040.1",columns=columns(TxDb.Athaliana)),keytype="TXNAME")
    ,label.color="red") # does not work



TxDb.Athaliana.selected<-select(TxDb.Athaliana,keys="AT1G01040.1",columns=columns(TxDb.Athaliana),keytype="TXNAME")
class(TxDb.Athaliana.selected) # this is data.frame
p3<-autoplot(TxDb.Athaliana.selected) # does not work because TxDb.Athaliana.selected is a data.frame.
#transcripts(TxDb.Athaliana,select(TxDb.Athaliana,keys="AT1G01040.1",columns=columns(TxDb.Athaliana),keytype="TXNAME"))
TxDb.Athaliana.selected

# 
transcripts(TxDb.Athaliana)
transcripts.Athaliana.selected2<-transcripts(TxDb.Athaliana,vals<-list(gene_id="AT1G01040"))


TxDb.Athaliana.selected2<-genes(TxDb.Athaliana,vals<-list(gene_id="AT1G01040"))
TxDb.Athaliana.selected2<-exons(TxDb.Athaliana,vals<-list(gene_id="AT1G01040"))
exons(TxDb.Athaliana)
GRL.Athaliana.exon<-exonsBy(TxDb.Athaliana,by="gene")

#TxDb.Athaliana.selected2<-disjointExons(TxDb.Athaliana)
autoplot(GRL.Athaliana.exon[4:5])

# another method
gr.Athaliana <- crunch(TxDb.Athaliana,which=transcripts.Athaliana.selected2)
## change column to 'model'
colnames(values(gr.Athaliana))[4] <- "model"
grl.Athaliana <- split(gr.Athaliana, gr.Athaliana$tx_name)
autoplot(grl.Athaliana)
autoplot(grl.Athaliana,aes(type=model)) # good!


# how to add additional plot for primer position?
## option1: by simply use ggplot function?
## option 2: by "2.6 building your tracks"
primer1 <- GRanges("chr2", IRanges(1e6+300, 1e6+330),strand="-")
primer2 <- GRanges("chr2", IRanges(1e6+120, 1e6+150),strand="+")
c(primer1,primer2)

autoplot(primer1) #error
p2<-autoplot(TxDb.Athaliana,which=GRanges("Chr2", IRanges(1e6, 1e6+1000)),label.color="red",fill="gray")
p2
p2 + geom_arrowrect(primer1,aes(fill=strand),arrow.head=0.5,rect.height=0.01) + geom_arrowrect(primer2,arrow.head=0.5,rect.height=0.01,aes(fill=strand))
p2 + geom_arrowrect(c(primer1,primer2),aes(fill=strand,y=100),arrow.head=0.5,rect.height=0.01)
# trial 2
p<-autoplot(gff[grep("AT1G01040",elementMetadata(gff)$group),])
p
#

### combine
primer1 <- GRanges("chr2", IRanges(1e6+300, 1e6+330),strand="-")
primer2 <- GRanges("chr2", IRanges(1e6+120, 1e6+150),strand="+")






# example in geom_arrowrect()
set.seed(1)
N <- 100
require(GenomicRanges)
## ======================================================================
##  simmulated GRanges
## ======================================================================
gr <- GRanges(seqnames = 
                sample(c("chr1", "chr2", "chr3"),
                       size = N, replace = TRUE),
              IRanges(
                start = sample(1:300, size = N, replace = TRUE),
                width = sample(70:75, size = N,replace = TRUE)),
              strand = sample(c("+", "-", "*"), size = N, 
                              replace = TRUE),
              value = rnorm(N, 10, 3), score = rnorm(N, 100, 30),
              sample = sample(c("Normal", "Tumor"), 
                              size = N, replace = TRUE),
              pair = sample(letters, size = N, 
                            replace = TRUE))


## ======================================================================
##  default
## ======================================================================
ggplot(gr) + geom_arrowrect()
## or
ggplot() + geom_arrowrect(gr)
## ======================================================================
##  facetting and aesthetics
## ======================================================================
ggplot(gr) + geom_arrowrect(facets = sample ~ seqnames, aes(color = strand, fill = strand))


## ======================================================================
##  stat:identity
## ======================================================================
ggplot(gr) + geom_arrowrect(stat = "identity", aes(y = value))


## ======================================================================
##  stat:stepping
## ======================================================================
ggplot(gr) + geom_arrowrect(stat = "stepping", aes(y = value, group = pair))
## ======================================================================
##  group.selfish controls when 
## ======================================================================
ggplot(gr) + geom_arrowrect(gr, stat = "stepping", aes(y = value, group = pair), group.selfish = FALSE)
ggplot(gr) + geom_arrowrect(gr, stat = "stepping", aes(y = value, group = pair))#, group.selfish = FALSE)


