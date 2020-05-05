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
