
pspl <- function(x, xl=min(x)-0.01, xr=max(x)+0.01,
                nknots=10, bdeg=3, pord=2, decom=2){
  dx <- (xr - xl)/nknots
  knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by=dx)
  B <- splines::spline.des(knots, x, bdeg+1, 0*x)$design
  a <- list(nknots=nknots,knots=knots,bdeg=bdeg,pord=pord,decom=decom)
  attributes(B) <- c(attributes(B),a)
  class(B) <- c("bs","basis","matrix")
  B
}

#####################################################################################

pspt <- function(sp1,sp2,time=NULL,scale=TRUE,ntime=NULL,
                xl_sp1=min(sp1)-0.01, xr_sp1=max(sp1)+0.01,
                xl_sp2=min(sp2)-0.01, xr_sp2=max(sp2)+0.01,
                xl_time=min(time)-0.01, xr_time=max(time)+0.01,
                nknots=c(10,10,5), bdeg=c(3,3,3),pord=c(2,2,2),
                decom=1,psanova=FALSE,
                nest_sp1=1,nest_sp2=1,nest_time=1,
                f1_main=TRUE,f2_main=TRUE,ft_main=TRUE,
                f12_int=TRUE,f1t_int=TRUE,f2t_int=TRUE,
                f12t_int=TRUE){

  if (length(sp1) != length(sp2)) stop("variables must have same length")
  nsp <- length(sp1)
  if (!is.null(time)){
    if (length(sp1) != length(time)) stop("variables must have same length")
    if (is.null(ntime)) stop("ntime is needed as argument")
  }
  if (scale){
    sp1 <- scale(sp1); sp2 <- scale(sp2)
    xl_sp1=min(sp1)-0.01; xr_sp1=max(sp1)+0.01
    xl_sp2=min(sp2)-0.01; xr_sp2=max(sp2)+0.01
    if (!is.null(time)) {
      time <- scale(time)
      xl_time=min(time)-0.01; xr_time=max(time)+0.01
    }
  }
  B <- NULL
  if (psanova){
    decom <- 2
    if (length(nest_sp1) != length(nest_sp2)){
      stop("nest_sp1 and nest_sp2 must have same length")
    }
    if (!is.null(time)){
      if (length(nest_sp1) != length(nest_time)) {
         stop("nest_sp and nest_time must have same length")
       }
    }
    if (length(nest_sp1)==1) nest_sp1 <- rep(nest_sp1,3)
    if (length(nest_sp2)==1) nest_sp2 <- rep(nest_sp2,3)
    if (!is.null(time)){ # Efectos Interacción hasta orden 3
      if (length(nest_time)==1) nest_time <- rep(nest_time,3)
      order_effects <- 1:3
    } else { # Efectos Interacción hasta orden 2
    order_effects <- 1:2
    nest_time <- rep(1,length(nest_sp1))
    }
    for (i in 1:length(order_effects)) {
       Bsp1 <- pspl(sp1,xl=xl_sp1,xr=xr_sp1,
                   nknots=nknots[1]/nest_sp1[i],
                   bdeg=bdeg[1],pord=pord[1],decom=decom)
       colnames(Bsp1) <- paste("Bsp1",1:ncol(Bsp1),sep=".")
       Bsp2 <- pspl(sp2,xl=xl_sp2,xr=xr_sp2,
                    nknots=nknots[2]/nest_sp2[i],
                    bdeg=bdeg[2],pord=pord[2],decom=decom)
       colnames(Bsp2) <- paste("Bsp2",1:ncol(Bsp2),sep=".")
       Bi <- cbind(Bsp1,Bsp2)
       if(!is.null(time)) {
         Btime <- pspl(time,xl=xl_time,xr=xr_time,
                      nknots=nknots[3]/nest_time[i],
                      bdeg=bdeg[3],pord=pord[3],decom=decom)
         colnames(Btime) <- paste("Btime",1:ncol(Btime),sep=".")
         Bi <- cbind(Bi,Btime)
       } else {
         ft_main <- f1t_int <- f2t_int <- f12t_int <- FALSE
         nest_time <- NULL; bdeg <- bdeg[1:2]; pord <- pord[1:2]
       }
       if(i==1) colnames(Bi) <- paste(colnames(Bi),"main",sep=".")
       if(i==2) colnames(Bi) <- paste(colnames(Bi),"int2ord",sep=".")
       if(i==3) colnames(Bi) <- paste(colnames(Bi),"int3ord",sep=".")
       B <- cbind(B,Bi)
    }
  } else { # PSANOVA=FALSE
    f1_main <- f2_main <- ft_main <- f12_int <- f1t_int <- FALSE
    f2t_int <- f12t_int <- FALSE
    Bsp1 <- pspl(sp1,xl=xl_sp1,xr=xr_sp1,nknots=nknots[1],
                 bdeg=bdeg[1],pord=pord[1],decom=decom)
    colnames(Bsp1) <- paste("Bsp1",1:ncol(Bsp1),sep=".")
    Bsp2 <- pspl(sp2,xl=xl_sp2,xr=xr_sp2, nknots=nknots[2],
                 bdeg=bdeg[2],pord=pord[2],decom=decom)
    colnames(Bsp2) <- paste("Bsp2",1:ncol(Bsp2),sep=".")
    B <- cbind(Bsp1,Bsp2)
    if (!is.null(time)){
       Btime <- pspl(time,xl=xl_time,xr=xr_time, nknots=nknots[3],
                 bdeg=bdeg[3],pord=pord[3],decom=decom)
       colnames(Btime) <- paste("Btime",1:ncol(Btime),sep=".")
       B <- cbind(B,Btime)
    } else { nest_time <- NULL; bdeg <- bdeg[1:2]; pord <- pord[1:2] }
  }
  a <- list(sp1=sp1,sp2=sp2,time=time,ntime=ntime,
            nknots=nknots, bdeg=bdeg,pord=pord,decom=decom,
            psanova=psanova,nest_sp1=nest_sp1,nest_sp2=nest_sp2,
            nest_time=nest_time,
            f1_main=f1_main,f2_main=f2_main,ft_main=ft_main,
            f12_int=f12_int,f1t_int=f1t_int,f2t_int=f2t_int,
            f12t_int=f12t_int)
  attributes(B) <- c(attributes(B),a)
  class(B) <- c("bs","basis","matrix")
  B
}

#####################################################################################

B_XZ <- function (B, x=NULL, pord=2, decom=1){
  m <- ncol(B)
  n <- nrow(B)
  D <- diff(diag(m), differences=pord)
  P.svd <- svd(crossprod(D))
  U <- (P.svd$u)[,1:(m-pord)] # eigenvectors
  d <- (P.svd$d)[1:(m-pord)]  # eigenvalues
  Z <- B %*% U
  if(decom == 1) {
    X <- B %*% ((P.svd$u)[,-(1:(m-pord))])
  } else if (decom == 2){
    X <- NULL
    for(i in 0:(pord-1)){ X <- cbind(X,x^i) }
  }
  # else if(decom == 3) {
  #   Xf <- NULL
  #   for(i in 0:(pord-1)){
  #     Xf <- cbind(Xf,knots[-c((1:pord),
  #         (length(knots)- pord + 1):length(knots))]^i)
  #   }
  #   X <- B %*% Xf
  # }
  list(X = X, Z = Z, d = d, B = B, m = m, D = D, U = P.svd$u)
}

#####################################################################################

Bspt <- function(sp1,sp2,time,nfull,ntime,psanova,
                 Bi,bdegspt){
  Bi_col <- colnames(Bi)
  Bsp1_col <- Bi_col[grepl("sp1",Bi_col)]
  Bsp2_col <- Bi_col[grepl("sp2",Bi_col)]
  Bsp1 <- Bi[,c(Bsp1_col)]
  Bsp2 <- Bi[,c(Bsp2_col)]
  rm(Bsp1_col,Bsp2_col)
  if (!is.null(time)){
    Btime_col <- Bi_col[grepl("time",Bi_col)]
    Btime <- Bi[,c(Btime_col)]
    time <- time[1:ntime]
    Btime <- Btime[1:ntime,]
    rm(Btime_col)
    if ((nfull %% ntime) != 0) stop("ntime is not a divisor of the sample size")
    seq_sp <- seq(from=1,to=nfull,by=ntime)
    sp1 <- sp1[seq_sp]
    sp2 <- sp2[seq_sp]
    nsp <- length(sp1)
    Bsp1 <- Bsp1[seq_sp,]
    Bsp2 <- Bsp2[seq_sp,]
  }  else Btime <- NULL
  rm(Bi,Bi_col)
  if (!psanova) { # PSANOVA=FALSE
    Bsptlist <- list(Bsp1=Bsp1, Bsp2=Bsp2, Btime=Btime)
    #spt_names <- c("sp1","sp2","time")
    rm(Bsp1,Bsp2,Btime)
  } else {# PSANOVA=TRUE
    Bsp1_col <- colnames(Bsp1)
    Bsp1_col_main <- Bsp1_col[grepl("main", Bsp1_col)]
    Bsp1_col_int2ord <- Bsp1_col[grepl("int2ord", Bsp1_col)]
    Bsp1_main <- Bsp1[,c(Bsp1_col_main)]
    Bsp1_int2ord <- Bsp1[,c(Bsp1_col_int2ord)]
    Bsp2_col <- colnames(Bsp2)
    Bsp2_col_main <- Bsp2_col[grepl("main", Bsp2_col)]
    Bsp2_col_int2ord <- Bsp2_col[grepl("int2ord", Bsp2_col)]
    Bsp2_main <- Bsp2[,c(Bsp2_col_main)]
    Bsp2_int2ord <- Bsp2[,c(Bsp2_col_int2ord)]
    if (!is.null(time)){
      Bsp1_col_int3ord <- Bsp1_col[grepl("int3ord", Bsp1_col)]
      Bsp1_int3ord <- Bsp1[,c(Bsp1_col_int3ord)]
      Bsp2_col_int3ord <- Bsp2_col[grepl("int3ord", Bsp2_col)]
      Bsp2_int3ord <- Bsp2[,c(Bsp2_col_int3ord)]
      Btime_col <- colnames(Btime)
      Btime_col_main <- Btime_col[grepl("main", Btime_col)]
      Btime_col_int2ord <- Btime_col[grepl("int2ord", Btime_col)]
      Btime_main <- Btime[,c(Btime_col_main)]
      Btime_int2ord <- Btime[,c(Btime_col_int2ord)]
      Btime_col_int3ord <- Btime_col[grepl("int3ord", Btime_col)]
      Btime_int3ord <- Btime[,c(Btime_col_int3ord)]
    } else {
      Bsp1_int3ord <- Bsp2_int3ord <- NULL
      Bsp1_col_int3ord <- Bsp2_col_int3ord <- NULL
      Btime_main <- Btime_int2ord <- Btime_int3ord <- NULL
      Btime_col_main <- Btime_col_int2ord <- Btime_col_int3ord <- NULL
      ft_main <- f1t_int <- f2t_int <- f12t_int <- FALSE
    }
    Bsptlist <- list(Bsp1_main=Bsp1_main,
                     Bsp2_main=Bsp2_main,
                     Btime_main=Btime_main,
                     Bsp1_int2ord=Bsp1_int2ord,
                     Bsp2_int2ord=Bsp2_int2ord,
                     Btime_int2ord=Btime_int2ord,
                     Bsp1_int3ord=Bsp1_int3ord,
                     Bsp2_int3ord=Bsp2_int3ord,
                     Btime_int3ord=Btime_int3ord)
  }
  res <- list(sp1=sp1,sp2=sp2,time=time,
              nsp=length(sp1),ntime=ntime,
              Bsptlist=Bsptlist)
}

#####################################################################################

B_XZ_spt <- function(sp1,sp2,time,pordspt,psanova,decomspt,
                     f1_main,f2_main,ft_main,
                     f12_int,f1t_int,f2t_int,f12t_int,Bsptlist){

  Xsptlist <- Zsptlist <- dsptlist <- csptlist <- list()
  names_Xj <- names_Zj <- names_dj <- names_cj <-NULL
  cont <- 1
  for (j in 1:length(names(Bsptlist)))
  {
    Bj <- Bsptlist[[j]]
    if (is.null(Bj))  next
    name_Bj <- names(Bsptlist)[j]
    spt_term <- sub("B","",name_Bj)
    is_sp1_term <- grepl("sp1",spt_term)
    is_sp2_term <- grepl("sp2",spt_term)
    is_time_term <- grepl("time",spt_term)
    if (is_sp1_term){
      x_term <- sp1
      pord_term <- pordspt[1]
    } else if (is_sp2_term) {
      x_term <- sp2
      pord_term <- pordspt[2]
    } else if (is_time_term){
      x_term <- time
      pord_term <- pordspt[3]
    }
    #if (!is.null(Bj)){
      #print(name_Bj)
      BtoXZ <- B_XZ(Bj,x=x_term,pord=pord_term,decom=decomspt)
      Xj <- BtoXZ$X
      names_Xj <- c(names_Xj,paste("X",spt_term,sep=""))
      Xsptlist[[cont]] <- Xj
      Zj <- BtoXZ$Z
      names_Zj <- c(names_Zj,paste("Z",spt_term,sep=""))
      Zsptlist[[cont]] <- Zj
      dj <- BtoXZ$d
      names_dj <- c(names_dj,paste(spt_term,sep=""))
      dsptlist[[cont]] <- dj
      names_cj <- c(names_cj,paste(spt_term,sep=""))
      csptlist[[cont]] <- ncol(Bj)
      cont <- cont + 1
    #} else cont <- NULL
  } # end for (j in 1:length(names(Bsptlist)))
  rm(Bj,name_Bj,x_term,pord_term,spt_term,cont)
  rm(is_sp1_term,is_sp2_term,is_time_term)

  names(Xsptlist) <- names_Xj
  names(Zsptlist) <- names_Zj
  names(dsptlist) <- names_dj
  names(csptlist) <- names_cj
  rm(Xj,Zj,dj,names_Xj,names_Zj,names_dj,names_cj)

  # REVISAR ESTA PARTE....
  # g1u <- rep(1,pord[2])%x%d1
  # g2u <- d2%x%rep(1,pord[1])
  # g1b <- rep(1,c2-pord[2])%x%d1
  # g2b <- d2%x%rep(1,c1-pord[1])

  if (!psanova){  # PS-ANOVA=FALSE
    X1 <- Xsptlist$Xsp1; X2 <- Xsptlist$Xsp2
    Z1 <- Zsptlist$Zsp1; Z2 <- Zsptlist$Zsp2
    if (!is.null(time)) { #SPATIO-TEMPORAL TREND. NO PS-ANOVA
      Xt <- Xsptlist$Xtime
      Zt <- Zsptlist$Ztime
      Xspt <- kronecker(Rten2(X1,X2),Xt)
      colnames(Xspt) <- paste("Xspt",1:ncol(Xspt),sep=".")
      Zspt <- cbind(kronecker(Rten2(Z1,X2), Xt),
                    kronecker(Rten2(X1,Z2), Xt),
                    kronecker(Rten2(X1,X2), Zt),
                    kronecker(Rten2(Z1,Z2), Xt),
                    kronecker(Rten2(Z1,X2), Zt),
                    kronecker(Rten2(X1,Z2), Zt),
                    kronecker(Rten2(Z1,Z2), Zt))
      colnames(Zspt) <- paste("Zspt",1:ncol(Zspt),sep=".")
      rm(Xt,Zt)
    } else { #SPATIAL TREND. NO PS-ANOVA
      Xspt <- Rten2(X2, X1)
      colnames(Xspt) <- paste("Xspt",1:ncol(Xspt),sep=".")
      Zspt <- cbind(Rten2(X2, Z1), Rten2(Z2, X1),Rten2(Z2, Z1))
      colnames(Zspt) <- paste("Zspt",1:ncol(Zspt),sep=".")
    }
    rm(X1,Z1)
  } else { # PS-ANOVA=TRUE
    X1 <- Xsptlist$Xsp1_main; X2 <- Xsptlist$Xsp2_main
    Z1 <- Zsptlist$Zsp1_main; Z2 <- Zsptlist$Zsp2_main
    Z1.2 <- Zsptlist$Zsp1_int2ord; Z2.2 <- Zsptlist$Zsp2_int2ord
    one1 <- X1[,1,drop=FALSE]
    one2 <- X2[,1,drop=FALSE]
    x1 <- X1[,-1,drop=FALSE]
    x2 <- X2[,-1,drop=FALSE]
    if (!is.null(time)) { #SPATIO-TEMPORAL TREND. PS-ANOVA
      Xt <- Xsptlist$Xtime_main;
      Zt <- Zsptlist$Ztime_main; Zt.2 <- Zsptlist$Ztime_int2ord
      Z1.3 <- Zsptlist$Zsp1_int3ord; Z2.3 <- Zsptlist$Zsp2_int3ord
      Zt.3 <- Zsptlist$Ztime_int3ord
      onet <- Xt[,1,drop=FALSE]
      xt <- Xt[,-1,drop=FALSE]
      Xones <- kronecker(Rten2(one1,one2),onet)
      X_f1_main <-  kronecker(Rten2(x1,one2),onet)
      X_f2_main <-  kronecker(Rten2(one1,x2),onet)
      X_ft_main <-  kronecker(Rten2(one1,one2),xt)
      X_f12_int <- kronecker(Rten2(x1,x2),onet)
      X_f1t_int <- kronecker(Rten2(x1,one2),xt)
      X_f2t_int <- kronecker(Rten2(one1,x2),xt)
      X_f12t_int <- kronecker(Rten2(x1,x2),xt)
      colnames(Xones) <- c("Intercept")
      colnames(X_f1_main) <- paste("X_f1_main",1:ncol(X_f1_main),sep=".")
      colnames(X_f2_main) <- paste("X_f2_main",1:ncol(X_f2_main),sep=".")
      colnames(X_ft_main) <- paste("X_ft_main",1:ncol(X_ft_main),sep=".")
      colnames(X_f12_int) <- paste("X_f12_int",1:ncol(X_f12_int),sep=".")
      colnames(X_f1t_int) <- paste("X_f1t_int",1:ncol(X_f1t_int),sep=".")
      colnames(X_f2t_int) <- paste("X_f2t_int",1:ncol(X_f2t_int),sep=".")
      colnames(X_f12t_int) <- paste("X_f12t_int",1:ncol(X_f12t_int),sep=".")
      Xspt <- Xones
      if(f1_main) Xspt <- cbind(Xspt,X_f1_main)
      if(f2_main) Xspt <- cbind(Xspt,X_f2_main)
      if(ft_main) Xspt <- cbind(Xspt,X_ft_main)
      if(f12_int) Xspt <- cbind(Xspt,X_f12_int)
      if(f1t_int) Xspt <- cbind(Xspt,X_f1t_int)
      if(f2t_int) Xspt <- cbind(Xspt,X_f2t_int)
      if(f12t_int) Xspt <- cbind(Xspt,X_f12t_int)
      Z_f1_main <- kronecker(Rten2(Z1, one2), onet) # g1u
      Z_f2_main <- kronecker(Rten2(one1, Z2), onet) # g2u
      Z_ft_main <- kronecker(Rten2(one1, one2), Zt) # g3u
      colnames(Z_f1_main) <- paste("Z_f1_main",1:ncol(Z_f1_main),sep=".")
      colnames(Z_f2_main) <- paste("Z_f2_main",1:ncol(Z_f2_main),sep=".")
      colnames(Z_ft_main) <- paste("Z_ft_main",1:ncol(Z_ft_main),sep=".")

      # g12u | g21u | g12b+g21b
      Z_f12_int <- cbind(kronecker(Rten2(Z1.2, x2), onet),
                     kronecker(Rten2(x1, Z2.2), onet),
                     kronecker(Rten2(Z1.2, Z2.2), onet))
      colnames(Z_f12_int) <- paste("Z_f12_int",1:ncol(Z_f12_int),sep=".")
      # g13u | g31u | g13b+g31b
      Z_f1t_int <- cbind(kronecker(Rten2(Z1.2, one2), xt),
                     kronecker(Rten2(x1, one2), Zt.2),
                     kronecker(Rten2(Z1.2, one2), Zt.2))
      colnames(Z_f1t_int) <- paste("Z_f1t_int",1:ncol(Z_f1t_int),sep=".")
      # g23u | g32u | g23b+g32b
      Z_f2t_int <- cbind(kronecker(Rten2(one1, Z2.2), xt),
                     kronecker(Rten2(one1, x2), Zt.2),
                     kronecker(Rten2(one1, Z2.2), Zt.2))
      colnames(Z_f2t_int) <- paste("Z_f2t_int",1:ncol(Z_f2t_int),sep=".")
      # g123u | g213u | g321u | g123b+g213b | g132b+g312b |
      # g231b+g321b | g1t+g2t+g3t
      Z_f12t_int <- cbind(kronecker(Rten2(Z1.3, x2), xt),
                      kronecker(Rten2(x1, Z2.3), xt),
                      kronecker(Rten2(x1, x2), Zt.3),
                      kronecker(Rten2(Z1.3, Z2.3), xt),
                      kronecker(Rten2(Z1.3, x2), Zt.3),
                      kronecker(Rten2(x1, Z2.3), Zt.3),
                      kronecker(Rten2(Z1.3, Z2.3), Zt.3))
      colnames(Z_f12t_int) <- paste("Z_f12t_int",1:ncol(Z_f12t_int),sep=".")
      Zspt <- NULL
      if(f1_main) Zspt <- cbind(Zspt,Z_f1_main)
      if(f2_main) Zspt <- cbind(Zspt,Z_f2_main)
      if(ft_main) Zspt <- cbind(Zspt,Z_ft_main)
      if(f12_int) Zspt <- cbind(Zspt,Z_f12_int)
      if(f1t_int) Zspt <- cbind(Zspt,Z_f1t_int)
      if(f2t_int) Zspt <- cbind(Zspt,Z_f2t_int)
      if(f12t_int) Zspt <- cbind(Zspt,Z_f12t_int)
      rm(X_f1_main,X_f2_main,X_ft_main,X_f12_int,X_f1t_int,X_f2t_int,X_f12t_int)
      rm(Z_f1_main,Z_f2_main,Z_ft_main,Z_f12_int,Z_f1t_int,Z_f2t_int,Z_f12t_int)
      rm(Xt,Zt,Zt.2,Z1.3,Z2.3,Zt.3)
    } else { # SPATIAL TREND. PS-ANOVA
      Xones <- Rten2(one1,one2)
      X_f1_main <- Rten2(x1,one2)
      X_f2_main <- Rten2(one1,x2)
      X_f12_int <- Rten2(x1,x2)
      colnames(Xones) <- c("Intercept")
      colnames(X_f1_main) <- paste("X_f1_main",1:ncol(X_f1_main),sep=".")
      colnames(X_f2_main) <- paste("X_f2_main",1:ncol(X_f2_main),sep=".")
      colnames(X_f12_int) <- paste("X_f12_int",1:ncol(X_f12_int),sep=".")
      Xspt <- Xones
      if(f1_main) Xspt <- cbind(Xspt,X_f1_main)
      if(f2_main) Xspt <- cbind(Xspt,X_f2_main)
      if(f12_int) Xspt <- cbind(Xspt,X_f12_int)
      Z_f1_main <- Rten2(Z1, one2)
      Z_f2_main <- Rten2(one1, Z2)
      colnames(Z_f1_main) <- paste("Z_f1_main",1:ncol(Z_f1_main),sep=".")
      colnames(Z_f2_main) <- paste("Z_f2_main",1:ncol(Z_f2_main),sep=".")
      Z_f12_int <- cbind(Rten2(Z1.2, x2),Rten2(x1, Z2.2),Rten2(Z1.2, Z2.2))
      colnames(Z_f12_int) <- paste("Z_f12_int",1:ncol(Z_f12_int),sep=".")
      Zspt <- NULL
      if(f1_main) Zspt <- cbind(Zspt,Z_f1_main)
      if(f2_main) Zspt <- cbind(Zspt,Z_f2_main)
      if(f12_int) Zspt <- cbind(Zspt,Z_f12_int)
      rm(X_f1_main,X_f2_main,X_f12_int)
      rm(Z_f1_main,Z_f2_main,Z_f12_int)
    }
    rm(X1,X2,Z1,Z2,Z1.2,Z2.2)
  }
  res <- list(Xsptlist=Xsptlist,
              Zsptlist=Zsptlist,
              dsptlist=dsptlist,
              csptlist=csptlist,
              Xspt=Xspt,
              Zspt=Zspt)
}

