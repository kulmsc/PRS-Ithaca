splitvalidate.mine <- function(ls.pipeline, test.bfile=NULL, 
                                       keep=NULL, remove=NULL, 
                                       pheno=NULL, covar=NULL, 
                                       trace=1, split=NULL, 
                                       rematch=!is.null(test.bfile)) {
  stopifnot(class(ls.pipeline) == "lassosum.pipeline")
  
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  
  rematch <- rematch # Forces an evaluation at this point
  if(is.null(test.bfile)) {
    test.bfile <- ls.pipeline$test.bfile
    keep.through.pheno <- !is.null(pheno) && 
      ((is.data.frame(pheno)) || 
         (is.character(pheno) && length(pheno) == 1))
    if(is.null(keep) && is.null(remove) && !keep.through.pheno)
      keep <- ls.pipeline$keep.test
  }
  
  ### Pheno & covar ### 
  parsed.test <- parseselect(test.bfile, keep=keep, remove=remove, export=TRUE)
  phcovar <- parse.pheno.covar(pheno=pheno, covar=covar, parsed=parsed.test,trace=trace)
  parsed.test <- phcovar$parsed
  pheno <- phcovar$pheno
  covar <- phcovar$covar
  print(str(parsed.test))
  # recal <- !identical(ls.pipeline$test.bfile, test.bfile) || 
  #   !identical(parsed.test$keep, ls.pipeline$keep.test)
  
  ### Split ###
  if(is.null(split)) {
    split <- sample(1:parsed.test$n %% 2 + 1)
  } else {
    stopifnot(length(split) == parsed.test$n)
    stopifnot(all(sort(unique(split)) == 1:2))
  }

  ### Split-validation ###
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  best.s <- best.lambda <- best.validation.result <- numeric(0)
  best.pgs <- pheno * NA
  best.beta <- numeric(0)
  validation.table <- data.frame()
  PGS <- list()
  for(s in 1:2) {
    if(is.null(parsed.test$keep)) {
      keep <- split == s
      pheno2 <- pheno[keep]
      covar2 <- if(!is.null(covar)) covar[keep,] else NULL
    } else {
      keep <- parsed.test$keep
      keeps <- split == s
      keep[keep] <- keeps
      pheno2 <- pheno[keeps]
      covar2 <- if(!is.null(covar)) covar[keeps,] else NULL
    }
    if(trace) cat(paste0("Split ", s, ":\n")) 
    v <- validate(ls.pipeline, keep=keep, pheno=pheno2, covar=covar2, 
                  test.bfile=test.bfile, trace=trace, rematch=rematch)
    best.s <- c(best.s, v$best.s)
    best.lambda <- c(best.lambda, v$best.lambda)
    PGS[[s]] <- v$pgs
    best.beta <- cbind(best.beta, v$best.beta)
    validation.table <- rbind(validation.table, v$validation.table)
    # best.validation.result <- c(best.validation.result, v$best.validation.result)
  }
  S <- v$s; L <- v$lambda
  for(s in 1:2) {
    best.pgs[split == 3-s] <- PGS[[3-s]][[which(S == best.s[s])]][,L == best.lambda[s]]
  }


  print(head(parsed.test$fam))
  print(is.null(parsed.test$fam))

  #### Results table ####
  if(is.null(phcovar$table)) {
    #results.table <- (if(is.null(parsed.test$fam)) read.table2(parsed.test$famfile) else
    #  parsed.test$fam)[,1:2]
    results.table <- read.table2(parsed.test$famfile)[,1:2]
    if(!is.null(parsed.test$keep)) results.table <- results.table[parsed.test$keep,]
    colnames(results.table) <- c("FID", "IID")
    results.table$pheno <- pheno
    results.table$best.pgs <- best.pgs
  } else {
    results.table <- phcovar$table
    results.table$best.pgs <- best.pgs[results.table$order]
    results.table$split <- split[results.table$order]
  }
  
    
  results <- c(results, list(split=split,
                             best.s=best.s, 
                             best.lambda=best.lambda,
                             best.pgs=best.pgs, 
                             best.beta=best.beta, 
                             validation.table=validation.table, 
                             validation.type=v$validation.type, 
                             pheno=pheno, 
                             best.validation.result=best.validation.result, 
                             results.table=results.table))
  class(results) <- "validate.lassosum"
  return(results)
  
}



parse.pheno.covar <- function(pheno, covar, parsed, trace=0) {
  #' @keywords internal

  fam <- parsed[['fam']]
  keep <- parsed$keep
  # keep <- NULL
  pheno.df <- NULL
  
  update.keep <- function(old, new) {
    if(all(new)) {
      return(old)
    } else {
      if(is.null(old)) return(new) else {
        if(is.null(new)) return(old) else 
          return(old & new)
      }
    }
  }
  #### pheno ####
  if(!is.null(pheno) && is.character(pheno) && length(pheno) == 1) {
    if(file.exists(pheno)) pheno <- read.table2(pheno, header=T) else 
      stop(paste("Cannot find", pheno))
  }
  if(is.data.frame(pheno)) {
    if(ncol(pheno) != 3) {
      stop(paste("A pheno data.frame must have 3 columns exactly",
                 "with the first 2 with headers 'FID' and 'IID'"))
    }
    colnames <- colnames(pheno) 
    if(!all(colnames[1:2] == c("FID", "IID"))) {
      stop(paste("The first two columns of the pheno", 
                 "data.frame must have headers 'FID' and 'IID'"))
    }
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    pheno.df <- pheno
    colnames(pheno.df)[3] <- "pheno"
    rownames(pheno) <- paste(pheno$FID, pheno$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(pheno))
    Pheno <- pheno[,3]
    names(Pheno) <- rownames(pheno)
  } else {
    if(!is.null(pheno)) {
      stopifnot(length(pheno) == parsed$n)
    } else {
      fam <- read.table2(parsed$famfile)
      if(is.null(parsed$keep)) pheno <- fam$V6 else 
        pheno <- fam$V6[parsed$keep]
    }
  }
  
  #### covar ####
  user.covar <- FALSE
  if(!is.null(covar) && is.character(covar) && length(covar) == 1) {
    if(file.exists(covar)) covar <- read.table2(covar, header=T) else 
      stop(paste("Cannot find", covar))
  }
  if(is.data.frame(covar) & all(colnames(covar)[1:2] == c("FID", "IID"))) {
    user.covar <- TRUE
    colnames <- colnames(covar) 
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    rownames(covar) <- paste(covar$FID, covar$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(covar))
    Covar <- covar[,-(1:2), drop=FALSE]
  } else {
    if(!is.null(covar)) {
      if(is.vector(covar)) covar <- matrix(covar, ncol=1)
      if(is.matrix(covar)) covar <- as.data.frame(covar)
      Covar <- covar
    } 
  }
  
  #### updates ####
  parsed$keep <- update.keep(parsed$keep, keep)
  if(!is.null(parsed$keep)) parsed$n <- sum(parsed$keep)
  if(is.data.frame(pheno)) {
    if(!is.null(parsed$keep)) {
      names <- rownames(fam)[parsed$keep]
    } else {
      names <- rownames(fam)
    }
    pheno <- Pheno[names] 
    if(trace) {
      message(length(pheno), " out of ", length(Pheno), " samples kept in pheno.")
      # message(paste("Note that the order of best.pgs is the order given in the .fam file", 
      #               " rather than the pheno data.frame. Use v$best.pgs[v$order] to get", 
      #               " the pgs in the order of the phenotype."))
    }
    Order <- 1:length(pheno)
    names(Order) <- names
    pheno.df$order <- Order[names(Pheno)]
  } 

  if(user.covar) {
    if(!is.null(parsed$keep)) covar <- Covar[rownames(fam)[parsed$keep],,drop=F] else 
      covar <- Covar[rownames(fam),,drop=F]
    if(trace) message(nrow(covar), " out of ", nrow(Covar), " samples kept in covar.")
  } 
  
  if(length(pheno) == 0) {
    stop("No phenotype left. Perhaps the FID/IID do not match?")
  } else if(length(pheno) != parsed$n) {
    stop("The length of pheno does not match the number of samples.")
  }
  if(!is.null(covar) && nrow(covar) != parsed$n) {
    stop(paste("The dimension of the covar matrix does not match the number of samples used.", 
               "If your covariate is a data.frame with FID and IID, make sure they have headers."))
  }
  # if(sd(pheno, na.rm = TRUE) == 0) stop("There's no variation in phenotype")
  parsed$fam <- fam

  return(list(pheno=pheno, covar=covar, parsed=parsed, table=pheno.df))
  
}


read.table2 <- function(file, header=F, data.table=F, check.names=TRUE, ...) {
  return(data.table::fread(file, header=header, data.table=data.table, 
                           check.names=check.names, ...))
}
