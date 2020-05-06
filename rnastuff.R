# ------------------------------
#       basics  - IDE
# ------------------------------

# what you see around you (IDE):
# console (command line, basically raw R language)
# environment (saved variables)
# help/viewer/files/plots
# code (where we are typing now)

# --------------------------------
#     basics - math and variables
# ------------------------------
# math:
4 + 3

# we obviously grew up with formulas
# a + b = c

# which in this case we can write as:
a <- 4
b <- 3
a + b

# what we've done here is assigned values to variabels.
# you can also assign it with '=', but that is bad practice.

# --------------------------------
#     basics - Types of variables
# ----------------------------------

vector <- c(1,3,7,'A','B')
matrix <- matrix(data = NA,nrow = 3,ncol=3)
df <- data.frame(matrix, row.names = NULL,stringsAsFactors = F)

lists<- list(vector, matrix, df)

# -------------------------------
#   Basics - functions
# ------------------------------

calcu <- function(a, b) {
  c <- a + b
  return(c)
}

calcu(a = 3, b = 4)


# -----------------------
# Basics - Packages
# -----------------------

# Hundreds of packages of functions are created with various forms and goals.
# for our goal we need 3 packages installed.

install.packages('edgeR')
install.packages('tweeDEseqCountData')
instal.packages('statmod')

install.packages('TCGAutils')
install.packages('data.table')

library(edgeR)
library(tweeDEseqCountData)
library(statmod)
library(data.table)
library(TCGAutils)

# edgeR is the RNA seq analysis package.
# tweeDeseqCountData and statmode are two packages to support edgeR

# data.table is a modified, more efficient and reliabel take on data.frames.
# it is especially good on your RAM memory management.

# ---------------------------------------------------
#   RNA seq differential gene expression analysis
# ---------------------------------------------------

# ---- Getting data ----
# TCGA-MESO:
# https://portal.gdc.cancer.gov/projects/TCGA-MESO

# TCGA-PAAD: 
# https://portal.gdc.cancer.gov/projects/TCGA-PAAD

# download:
createRNAdat <- function(filepath, ext="*.counts.gz") {
  
  library(TCGAutils)

  RNAlist <- list.files(path = filepath, recursive = T,
                        full.names = F,
                        pattern = ext)
  # get patient names:
  names <- lapply(strsplit(unlist(RNAlist),split = '/'),
                  function(x){ save <- x[2] })
  patname <- UUIDtoBarcode(unlist(names),from_type = 'file_id')
  setwd(filepath)
  rawRNA <- fread(file = RNAlist[1])
  rows <- rawRNA[,grep(pattern = "ENSG*",x = V1)]
  ENSMBL_RNA_namelist <- rawRNA[rows, c(V1)]
  
  # remove 2 digit version info
  nENSMBL_RNA_name <-
    sapply(ENSMBL_RNA_namelist, function(x) {
      sub("\\..*", "", x)
    }, USE.NAMES = F, simplify = T)
  
  for (i in RNAlist) {
    #import and merge all the reads from patients.
    if (!exists("rawRNA")) {
      rawRNA <- fread(file = i, col.names =  c("ENSMBLname", i))
    } else
      if (exists("rawRNA")) {
        temp_RNA <- fread(file = i,
                           # drop = c((tail(rows,n = 1)+1) : nrow(rawRNA)),
                           select = 2,
                           col.names = i
                          )
        rawRNA <- cbind(rawRNA, temp_RNA)
        rm(temp_RNA)
      }
  }

  
  RNA_dat <- rawRNA[rows,-2]
  RNA_dat[, V1 := nENSMBL_RNA_name]
  colnames(RNA_dat) <- c("nENSMBL_RNA_name", patname[[2]])
  # colnames(RNA_dat) <- sub(x = colnames(RNA_dat), pattern = ".tsv",
  #                             replacement = "")
  return(RNA_dat)
}


setwd("~/Documents/TCGA-PAAD/Data/RNA_htseq_counts/")

fwrite(RNADAT,file = "~/Documents/TCGA-PAAD/Data/htseq_counts_comb.csv", sep = ",")
# ---- cleaning data -----


RNADAT <- createRNAdat(filepath = "~/Documents/TCGA-PAAD/Data/RNA_htseq_counts/")
RNADAT[1:10, 1:10]


# ---- creating groups ----

# a little bit on data.table

DT[i , j , by]

# take DT, subset rows by 'i', then calculate 'j' grouped by 'by'
# https://www.ensembl.org/index.html
ASS1 <- "ENSG00000130707"

assdt <- RNADAT[nENSMBL_RNA_name == ASS1,-1]

plot(unlist(assdt, use.names = F))
summary(unlist(assdt, use.names = F))


# how many transcripts does each patient have? normalize to that.
sumreads <- RNADAT[,colSums(.SD), .SDcols=2:length(RNADAT)]
plot(sumreads)
plot(unlist(assdt, use.names = F))

new_assdt <- unlist(assdt / sumreads)
plot(new_assdt)
summary(new_assdt)

# This is our border.
median_assdt <- summary(new_assdt)[3] 
# Now create patients into two groups: high and low.

# our first test! 
high <- names(new_assdt[new_assdt > median_assdt])
low <- names(new_assdt[new_assdt <= median_assdt])


# ---- performing differential gene expression analysis -----

RNAmat <- as.matrix(RNADAT[,-c(1)])

# create a DGEList object which stores read counts and associated information such as gene annotation.
h <- DGEList(counts=RNAmat, genes=RNADAT[[1]]) 
isexpr <- rowSums(cpm(h)>1) >= 20  # keep reads with at least 1cpm in at least 20 samples.

hasannot <- rowSums(is.na(h$genes))==0  # keep only genes with annotations
h <- h[isexpr & hasannot, ,keep.lib.sizes=F]
dim(h)

png(filename= "q1q4 HLA-comb library size in million reads")
barplot(h$samples$lib.size*1e-6, ylab="Library size (millions)")
h <- calcNormFactors(h)


grp <- rownames(h$sample)
grp <- as.data.frame(grp)
grp[,"val"] <- NA
grp$val <- ifelse(grp$grp %in% high, "HIGH", "LOW")
table(grp$val)
Grouping <- grp$val


design <- model.matrix(~Grouping)
h <- estimateDisp(h, design, robust=TRUE)
tiff(filename= "Estimate of HLA combined q1q4 negative binomial distribution", width = 720, height = 720)
plotBCV(h)

fit <- glmQLFit(h, design, robust=TRUE)
tiff(filename = "QL dispersions HLA Combined q1q4", width = 720, height = 720)
plotQLDisp(fit)

#top differentially expressed genes:
qlf <- glmQLFTest(fit)
topTags(qlf,n=50)
summary(decideTests(qlf))
write.csv(topTags(qlf,n=50), file = "differential gene expression q1q4 HLA comb reads.csv")



# ---- Getting data ----
# TCGA-MESO:
# https://portal.gdc.cancer.gov/projects/TCGA-MESO

# TCGA-PAAD: 
# https://portal.gdc.cancer.gov/projects/TCGA-PAAD
