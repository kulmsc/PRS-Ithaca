my.subset.bigSNP <- function(x,
                          ind.row = rows_along(G),
                          ind.col = cols_along(G),
                          ...) {

  G <- x$genotypes
  # Support for negative indices
  ind.row <- rows_along(G)[ind.row]
  ind.col <- cols_along(G)[ind.col]

  #check_args()
  getNewFile <- function(x, type) {

    root <- sub("\\.bk$", "", x$genotypes$backingfile)
    EXTS <- c("bk", "rds")

    number <- 1
    repeat {
      files <- sprintf("%s_%s%d.%s", root, type, number, EXTS)
      if (all(!file.exists(files))) break
      number <- number + 1
    }

    sprintf("%s_%s%d", root, type, number)
  }

  replaceSNP <- function(BM, BM2, rowInd, colInd) {
    invisible(.Call(`_bigsnpr_replaceSNP`, BM, BM2, rowInd, colInd))
  }


  # Create new FBM and fill it
  G2 <- FBM.code256(
    nrow = length(ind.row),
    ncol = length(ind.col),
    code = G$code256,
    init = NULL,
    backingfile = getNewFile(x, "sub"),
    create_bk = TRUE
  )
  replaceSNP(G2, G, rowInd = ind.row, colInd = ind.col)

  # http://stackoverflow.com/q/19565621/6103040
  newfam <- x$fam[ind.row, , drop = FALSE]
  rownames(newfam) <- rows_along(newfam)
  newmap <- x$map[ind.col, , drop = FALSE]
  rownames(newmap) <- rows_along(newmap)

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = G2,
                             fam = newfam,
                             map = newmap),
                        class = "bigSNP")

  # save it and return the path of the saved object
  #rds <- sub("\\.bk$", ".rds", G2$backingfile)
  rds <- "./subset.bigsnp.rds"
  saveRDS(snp.list, rds)
  rds
}

