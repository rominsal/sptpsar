fit_terms <- function(sptsarfit,variables){
  # Function to get fitted values of each variable
  # Inputs:
    # sptsarfit: Object of class psar fitted
    # variables: Vector of character including variables names
  X <- sptsarfit$X
  Z <- sptsarfit$Z
  mt <- sptsarfit$terms
  terms_mt <- attr(mt,"term.labels")
  bfixed <- sptsarfit$bfixed
  brandom <- sptsarfit$brandom
  psanova <- sptsarfit$psanova
  cov_b <-sptsarfit$vcov_b
  row_cov_fixed <- c(grepl("fixed",rownames(cov_b)))
  col_cov_fixed <- c(grepl("fixed",colnames(cov_b)))
  cov_b_fixed <- cov_b[row_cov_fixed,col_cov_fixed]
  row_cov_random <- c(grepl("random",rownames(cov_b)))
  col_cov_random <- c(grepl("random",colnames(cov_b)))
  cov_b_random <- cov_b[row_cov_random,col_cov_random]
  fitted_terms_fixed <- fitted_terms_random <- fitted_terms <- NULL
  se_fitted_terms <- se_fitted_terms_fixed <- se_fitted_terms_random <- NULL
  # var_fitted_terms <- list()
  # var_fitted_terms_fixed <- list()
  # var_fitted_terms_random <- list()
  for (i in 1:length(variables)) {
    var_name <- variables[i]
    if (grepl("spttrend",var_name)) {
       time <- sptsarfit$time
       if (psanova){# psanova=TRUE
         if (!is.null(time)){ # 3d
           eff_spttrend_psanova <- c("Intercept","f1_main","f2_main","ft_main",
                                     "f12_int", "f1t_int", "f2t_int", "f12t_int")
         } else { # 2d
           eff_spttrend_psanova <- c("Intercept","f1_main","f2_main","f12_int")
         }
         for (j in 1:length(eff_spttrend_psanova)) {
           eff_spttrend_psanova_j <- eff_spttrend_psanova[j]
           Xj <- as.matrix(X[,grepl(eff_spttrend_psanova_j,colnames(X))])
           bfixed_j <- bfixed[grepl(eff_spttrend_psanova_j,names(bfixed))]
           Zj <- as.matrix(Z[,grepl(eff_spttrend_psanova_j,colnames(Z))])
           brandom_j <- brandom[grepl(eff_spttrend_psanova_j,names(brandom))]
           term_fixed_j <- Xj %*% bfixed_j
           term_random_j <- Zj %*% brandom_j
           term_j <- term_fixed_j + term_random_j
           names(term_fixed_j) <- names(term_random_j) <- eff_spttrend_psanova_j
           names(term_j) <- eff_spttrend_psanova_j
           fitted_terms_fixed <- Matrix::cBind(term_fixed_j, fitted_terms_fixed)
           fitted_terms_random <- Matrix::cBind(term_random_j, fitted_terms_random)
           fitted_terms <- Matrix::cBind(term_j, fitted_terms)
           colnames(fitted_terms_fixed)[1] <- eff_spttrend_psanova_j
           colnames(fitted_terms_random)[1] <- eff_spttrend_psanova_j
           colnames(fitted_terms)[1] <- eff_spttrend_psanova_j
           row_cov_j <- c(grepl(eff_spttrend_psanova_j,rownames(cov_b)))
           col_cov_j <- c(grepl(eff_spttrend_psanova_j,colnames(cov_b)))
           cov_b_j <- cov_b[row_cov_j,col_cov_j]
           se_term_j <- Matrix::rowSums( (Matrix::cBind(Xj,Zj) %*% cov_b_j)
                                        * Matrix::cBind(Xj,Zj) )^0.5
           se_fitted_terms <- Matrix::cBind(se_term_j,se_fitted_terms)
           colnames(se_fitted_terms)[1] <- eff_spttrend_psanova_j
           #var_fitted_terms[[j]] <- Matrix::cBind(Xj,Zj) %*%
           #                          (cov_b_j %*% Matrix::t(Matrix::cBind(Xj,Zj)))
           #names(var_fitted_terms)[j] <- eff_spttrend_psanova_j
           row_cov_j_fixed <- c(grepl(eff_spttrend_psanova_j,rownames(cov_b_fixed)))
           col_cov_j_fixed <- c(grepl(eff_spttrend_psanova_j,colnames(cov_b_fixed)))
           cov_b_j_fixed <- cov_b_fixed[row_cov_j_fixed,col_cov_j_fixed]
           se_term_j_fixed <- Matrix::rowSums( (Xj %*% cov_b_j_fixed)
                                               * Xj )^0.5
           se_fitted_terms_fixed <- Matrix::cBind(se_term_j_fixed,
                                          se_fitted_terms_fixed)
           colnames(se_fitted_terms_fixed)[1] <- eff_spttrend_psanova_j
           #var_fitted_terms_fixed[[j]] <- Xj %*% cov_b_j_fixed %*% t(Xj)
           #names(var_fitted_terms_fixed)[j] <-
           #                   paste("fixed_",eff_spttrend_psanova_j,sep="")
           row_cov_j_random <- c(grepl(eff_spttrend_psanova_j,rownames(cov_b_random)))
           col_cov_j_random <- c(grepl(eff_spttrend_psanova_j,colnames(cov_b_random)))
           cov_b_j_random <- cov_b_random[row_cov_j_random,col_cov_j_random]
           se_term_j_random <- Matrix::rowSums( (Zj %*% cov_b_j_random)
                                               * Zj )^0.5
           se_fitted_terms_random <- Matrix::cBind(se_term_j_random,
                                          se_fitted_terms_random)
           colnames(se_fitted_terms_random)[1] <- eff_spttrend_psanova_j
           #var_fitted_terms_random[[j]] <- Zj %*% cov_b_j_random %*% t(Zj)
           #names(var_fitted_terms_random)[j] <-
           #  paste("random_",eff_spttrend_psanova_j,sep="")
         }
         match_fixed <- unique(grep(paste(eff_spttrend_psanova,
                                               collapse="|"),
                                         names(bfixed), value=TRUE))
         bfixed_spt <- bfixed[match_fixed]
         match_random <- unique(grep(paste(eff_spttrend_psanova,
                                          collapse="|"),
                                    names(brandom), value=TRUE))
         brandom_spt <- brandom[match_random]
         match_X <- unique(grep(paste(eff_spttrend_psanova,
                                          collapse="|"),
                                    colnames(X), value=TRUE))
         X_spt <- X[,match_X]
         match_Z <- unique(grep(paste(eff_spttrend_psanova,
                                      collapse="|"),
                                colnames(Z), value=TRUE))
         Z_spt <- Z[,match_Z]
         fitted_term_fixed_spt <- X_spt %*% bfixed_spt
         fitted_term_random_spt <- Z_spt %*% brandom_spt
         fitted_term_spt <- fitted_term_fixed_spt + fitted_term_random_spt

         fitted_terms_fixed <- Matrix::cBind(fitted_term_fixed_spt,
                                     fitted_terms_fixed)
         fitted_terms_random <- Matrix::cBind(fitted_term_random_spt,
                                      fitted_terms_random)
         fitted_terms <- Matrix::cBind(fitted_term_spt, fitted_terms)
         colnames(fitted_terms_fixed)[1] <- "spttrend"
         colnames(fitted_terms_random)[1] <- "spttrend"
         colnames(fitted_terms)[1] <- "spttrend"
         match_cov_bspt <- unique(grep(paste(eff_spttrend_psanova,
                                     collapse="|"),
                               colnames(cov_b), value=TRUE))
         cov_b_spt <- cov_b[match_cov_bspt,match_cov_bspt]
         se_term_spt <- Matrix::rowSums( (Matrix::cBind(X_spt,Z_spt) %*% cov_b_spt)
                                        * Matrix::cBind(X_spt,Z_spt) )^0.5
         se_fitted_terms <- Matrix::cBind(se_term_spt,
                                  se_fitted_terms)
         colnames(se_fitted_terms)[1] <- "spttrend"
         # var_fitted_terms$spttrend <- Matrix::cBind(X_spt,Z_spt) %*%
         #  (cov_b_spt %*% Matrix::t(Matrix::cBind(X_spt,Z_spt)))
         match_cov_bspt_fixed <- unique(grep(paste(eff_spttrend_psanova,
                                             collapse="|"),
                                       colnames(cov_b_fixed), value=TRUE))
         cov_b_spt_fixed <- cov_b[match_cov_bspt_fixed,match_cov_bspt_fixed]
         se_term_spt_fixed <- Matrix::rowSums( (X_spt %*% cov_b_spt_fixed)
                                              * X_spt )^0.5
         se_fitted_terms_fixed <- Matrix::cBind(se_term_spt_fixed,
                                        se_fitted_terms_fixed)
         colnames(se_fitted_terms_fixed)[1] <- "spttrend"
         #var_fitted_terms_fixed$spttrend <- X_spt %*% cov_b_spt_fixed %*% t(X_spt)
         match_cov_bspt_random <- unique(grep(paste(eff_spttrend_psanova,
                                                   collapse="|"),
                                             colnames(cov_b_random), value=TRUE))
         cov_b_spt_random <- cov_b[match_cov_bspt_random,match_cov_bspt_random]
         se_term_spt_random <- Matrix::rowSums( (Z_spt %*% cov_b_spt_random)
                                               * Z_spt )^0.5
         se_fitted_terms_random <- Matrix::cBind(se_term_spt_random,
                                         se_fitted_terms_random)
         colnames(se_fitted_terms_random)[1] <- "spttrend"
         #var_fitted_terms_random$spttrend <- Z_spt %*% cov_b_spt_random %*% t(Z_spt)
      } else { # psanova=FALSE
        Xi <- as.matrix(X[,grepl("spt",colnames(X))])
        bfixed_i <- bfixed[grepl("spt",names(bfixed))]
        Zi <- as.matrix(Z[,grepl("spt",colnames(Z))])
        brandom_i <- brandom[grepl("spt",names(brandom))]
        term_fixed_i <- Xi %*% bfixed_i
        term_random_i <- Zi %*% brandom_i
        term_i <- term_fixed_i + term_random_i
        fitted_terms_fixed <- Matrix::cBind(term_fixed_i,fitted_terms_fixed)
        fitted_terms_random <- Matrix::cBind(term_random_i, fitted_terms_random)
        fitted_terms <- Matrix::cBind(term_i,fitted_terms)
        colnames(fitted_terms_fixed)[1] <- "spttrend"
        colnames(fitted_terms_random)[1] <- "spttrend"
        colnames(fitted_terms)[1] <- "spttrend"
        row_cov_i <- c(grepl("spt",rownames(cov_b)))
        col_cov_i <- c(grepl("spt",colnames(cov_b)))
        cov_b_i <- cov_b[row_cov_i,col_cov_i]
        se_term_i <- Matrix::rowSums( (Matrix::cBind(Xi,Zi) %*% cov_b_i)
                                     * Matrix::cBind(Xi,Zi) )^0.5
        se_fitted_terms <- Matrix::cBind(se_term_i,
                                 se_fitted_terms)
        colnames(se_fitted_terms)[1] <- "spttrend"
        #var_fitted_terms[[i]] <- Matrix::cBind(Xi,Zi) %*%
        #                         (cov_b_i %*% t(Matrix::cBind(Xi,Zi)))
        #names(var_fitted_terms)[i] <- "spttrend"
        row_cov_i_fixed <- c(grepl("spt",rownames(cov_b_fixed)))
        col_cov_i_fixed <- c(grepl("spt",colnames(cov_b_fixed)))
        cov_b_i_fixed <- cov_b_fixed[row_cov_i_fixed,col_cov_i_fixed]
        se_term_i_fixed <- Matrix::rowSums( (Xi %*% cov_b_i_fixed)
                                           * Xi )^0.5
        se_fitted_terms_fixed <- Matrix::cBind(se_term_i_fixed,
                                       se_fitted_terms_fixed)
        colnames(se_fitted_terms_fixed)[1] <- "spttrend"
        # var_fitted_terms_fixed[[i]] <- Xi %*% cov_b_i_fixed %*% t(Xi)
        # names(var_fitted_terms_fixed)[i] <-
        #   paste("fixed_","spttrend",sep="")
        row_cov_i_random <- c(grepl("spt",rownames(cov_b_random)))
        col_cov_i_random <- c(grepl("spt",colnames(cov_b_random)))
        cov_b_i_random <- cov_b_random[row_cov_i_random,col_cov_i_random]
        se_term_i_random <- Matrix::rowSums( (Zi %*% cov_b_i_random)
                                            * Zi )^0.5
        se_fitted_terms_random <- Matrix::cBind(se_term_i_random,
                                       se_fitted_terms_random)
        colnames(se_fitted_terms_random)[1] <- "spttrend"
        # var_fitted_terms_random[[i]] <- Zi %*% cov_b_i_random %*% t(Zi)
        # names(var_fitted_terms_random)[i] <-
        #   paste("random_","spttrend",sep="")
      }
    } else { # No spttrend
      Xi <- as.matrix(X[,grepl(var_name,colnames(X))])
      bfixed_i <- bfixed[grepl(var_name,names(bfixed))]
      Zi <- as.matrix(Z[,grepl(var_name,colnames(Z))])
      brandom_i <- brandom[grepl(var_name,names(brandom))]
      term_fixed_i <- Xi %*% bfixed_i
      term_random_i <- Zi %*% brandom_i
      term_i <- term_fixed_i + term_random_i
      fitted_terms_fixed <- Matrix::cBind(term_fixed_i,fitted_terms_fixed)
      fitted_terms_random <- Matrix::cBind(term_random_i, fitted_terms_random)
      fitted_terms <- Matrix::cBind(term_i,fitted_terms)
      colnames(fitted_terms_fixed)[1] <- var_name
      colnames(fitted_terms_random)[1] <- var_name
      colnames(fitted_terms)[1] <- var_name
      row_cov_i <- c(grepl(var_name,rownames(cov_b)))
      col_cov_i <- c(grepl(var_name,colnames(cov_b)))
      cov_b_i <- cov_b[row_cov_i,col_cov_i]
      se_term_i <- Matrix::rowSums( (Matrix::cBind(Xi,Zi) %*% cov_b_i)
                                    * Matrix::cBind(Xi,Zi) )^0.5
      se_fitted_terms <- Matrix::cBind(se_term_i,se_fitted_terms)
      colnames(se_fitted_terms)[1] <- var_name
      # var_fitted_terms[[i]] <- Matrix::cBind(Xi,Zi) %*%
      #                          (cov_b_i %*% t(Matrix::cBind(Xi,Zi)))
      # names(var_fitted_terms)[i] <- var_name
      row_cov_i_fixed <- c(grepl(var_name,rownames(cov_b_fixed)))
      col_cov_i_fixed <- c(grepl(var_name,colnames(cov_b_fixed)))
      cov_b_i_fixed <- cov_b_fixed[c(row_cov_i_fixed),c(col_cov_i_fixed)]
      se_term_i_fixed <- Matrix::rowSums( (Xi %*% cov_b_i_fixed)
                                          * Xi )^0.5
      se_fitted_terms_fixed <- Matrix::cBind(se_term_i_fixed,
                                     se_fitted_terms_fixed)
      colnames(se_fitted_terms_fixed)[1] <- var_name
      # var_fitted_terms_fixed[[i]] <- Xi %*% cov_b_i_fixed %*% t(Xi)
      # names(var_fitted_terms_fixed)[i] <-
      #   paste("fixed_",var_name,sep="")
      row_cov_i_random <- c(grepl(var_name,rownames(cov_b_random)))
      col_cov_i_random <- c(grepl(var_name,colnames(cov_b_random)))
      cov_b_i_random <- cov_b_random[row_cov_i_random,col_cov_i_random]
      se_term_i_random <- Matrix::rowSums( (Zi %*% cov_b_i_random)
                                           * Zi )^0.5
      se_fitted_terms_random <- Matrix::cBind(se_term_i_random,
                                      se_fitted_terms_random)
      colnames(se_fitted_terms_random)[1] <- var_name
      # var_fitted_terms_random[[i]] <- Zi %*% cov_b_i_random %*% t(Zi)
      # names(var_fitted_terms_random)[i] <-
      #   paste("random_",var_name,sep="")
    }
  } # end for (i in 1:length(variables))

  # se_fitted_terms <- matrix(0,nrow=nrow(fitted_terms),ncol=ncol(fitted_terms))
  # colnames(se_fitted_terms) <- colnames(fitted_terms)
  #
  # for (i in 1:ncol(se_fitted_terms)){
  #   se_fitted_terms[,i] <- sqrt(diag(as.matrix(var_fitted_terms[[i]])))
  # }
  # se_fitted_terms_fixed <- matrix(0,nrow=nrow(fitted_terms_fixed),
  #                                 ncol=ncol(fitted_terms_fixed))
  # colnames(se_fitted_terms_fixed) <- colnames(fitted_terms_fixed)
  # for (i in 1:ncol(se_fitted_terms_fixed)){
  #   se_fitted_terms_fixed[,i] <- sqrt(diag(as.matrix(var_fitted_terms_fixed[[i]])))
  # }
  # se_fitted_terms_random <- matrix(0,nrow=nrow(fitted_terms_random),
  #                                 ncol=ncol(fitted_terms_random))
  # colnames(se_fitted_terms_random) <- colnames(fitted_terms_random)
  # for (i in 1:ncol(se_fitted_terms_random)){
  #   se_fitted_terms_random[,i] <- sqrt(diag(as.matrix(var_fitted_terms_random[[i]])))
  # }
  res <- list(fitted_terms = as.matrix(fitted_terms),
              se_fitted_terms = as.matrix(se_fitted_terms),
              fitted_terms_fixed = as.matrix(fitted_terms_fixed),
              se_fitted_terms_fixed = as.matrix(se_fitted_terms_fixed),
              fitted_terms_random = as.matrix(fitted_terms_random),
              se_fitted_terms_random = as.matrix(se_fitted_terms_random) )
  # if (list_cov){
  #   res$lcov_fitted <- var_fitted_terms
  #   res$lcov_fitted_fixed <- var_fitted_terms_fixed
  #   res$lcov_fitted_random <- var_fitted_terms_random
  # }
  res
}
