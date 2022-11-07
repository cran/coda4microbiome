#' explore_zeros
#'
#' Provides the proportion of zeros for a pair of variables (taxa) in table x and the proportion of samples
#' with zero in both variables. A bar plot with this information is also provided.
#' Results can be stratified by a categorical variable.
#'
#'
#' @param x abundance table (rows are samples, columns are variables (taxa))
#' @param id1 column number in x for the first taxa
#' @param id2 column number in x for the second taxa
#' @param strata stratification variable (default = NULL)
#'
#' @return a list with the frequency table and the associated bar plot
#' @export
#' @import ggpubr
#' @import ggplot2
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(HIV, package = "coda4microbiome")
#'
#' explore_zeros(x_HIV,5,6)
#'
#' explore_zeros(x_HIV,5,6, strata=y_HIV)
#'
#'
#-------------------------------------------------------------------------------
explore_zeros<-function(x,id1,id2,strata=NULL){
  #options(warn=-1)
  # library(ggpubr)
  m <- nrow(x)
  n <- ncol(x)

  if (sum(x==0) == 0){
    print("there are no zeros")
    return()
  } else {
    if(is.null(strata)){ #no strata in the input
      x1 <- x[,id1]
      x1<-ifelse((x1 ==0),0,1)
      x2 <- x[,id2]
      x2<-ifelse((x2 ==0),0,1)
      x1x2 <- ifelse(((x1 ==0)&(x2==0)),0,1)
      tx1 <- table(x1)/m
      if (tx1[[1]]==1) tx1<-cbind(0,tx1)
      tx2 <- table(x2)/m
      if (tx2[[1]]==1) tx2<-cbind(0,tx2)
      tx1x2 <- table(x1x2)/m
      if (tx1x2[[1]]==1) tx1x2<-cbind(0,tx1x2)

      tfreq<-data.frame(variables=c("x1","x2","x1x2"),freq=c(0,0,0))
      tfreq[1,2] <- round(tx1[[1]],2)
      tfreq[2,2] <- round(tx2[[1]],2)
      tfreq[3,2] <- round(tx1x2[[1]],2)
      outP <- ggpubr::ggbarplot(tfreq, x="variables", y="freq",
                        color = "steelblue",
                        fill ="steelblue",   # fill the bars
                        label = TRUE,             # show labels (counts for each category)
                        lab.pos = "in",
                        lab.col = "white")
      print(outP)
    } else { #strata is given in the input
      y<- strata
      x1 <- x[,id1]
      x1<-ifelse((x1 ==0),0,1)
      x2 <- x[,id2]
      x2<-ifelse((x2 ==0),0,1)
      x1x2 <- ifelse(((x1 ==0)&(x2==0)),0,1)
      tx1y <- table(x1,y)
      tx1y<-prop.table(table(x1,y),2)
      if (dim(tx1y)[1]==1){
          if (rownames(tx1y)=="1") tx1y<-rbind(c(0,0),tx1y)
      }
      tx2y <- table(x2,y)
      tx2y<-prop.table(table(x2,y),2)
      if (dim(tx2y)[1]==1){
          if (rownames(tx2y)=="1") tx2y<-rbind(c(0,0),tx2y)
      }
      tx1x2y <- table(x1x2,y)
      tx1x2y<-prop.table(table(x1x2,y),2)
      if (dim(tx1x2y)[1]==1){
           if (rownames(tx1x2y)=="1") tx1x2y<-rbind(c(0,0),tx1x2y)
      }

      tfreq <- rbind(rbind(tx1y[1,],tx2y[1,]),tx1x2y[1,]);
      nom1 <- paste("x",id1,sep="")
      nom2 <- paste("x",id2,sep="")
      nom12 <- paste(paste("x",id1,sep=""),paste("x",id2,sep=""),sep="")
      rownames(tfreq)<-c(nom1,nom2,nom12)
      col0 <- rep(rownames(tfreq),2)
      col1 <-c(rep(colnames(tfreq)[1],3),rep(colnames(tfreq)[2],3))
      col2 <-c(tfreq[,1],tfreq[,2]) #relative frequencies
      freq <- data.frame(col0,col1,col2)
      freq[,3]<-round(freq[,3],2)
      colnames(freq) <- c("variables","y","freq")
      outP <- ggpubr::ggbarplot(freq, x = "variables", y = "freq",
                        fill = "y", color = "y",
                        label = TRUE,
                        lab.pos = "out",
                        lab.col = "black",
                        palette = "hue",
                        position = ggplot2::position_dodge(0.9))

    }
  }
  print(outP)
  outTab <- tfreq
  outVal <- list("freq_table" = outTab, "freq_plot"=outP )


  return(outVal)
}
#-------------------------------------------------------------------------------

#' impute_zeros
#'
#' Simple imputation: When the abundance table contains zeros, a positive value is added
#' to all the values in the table. It adds 1 when the minimum of table is larger than 1 (i.e. tables of counts) or
#' it adds half of the minimum value of the table, otherwise.
#'
#'
#' @param x abundance table (rows are samples, columns are variables (taxa))
#'
#' @return x abundance table with zeros substituted by imputed positive values
#'
#' @export
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(HIV, package = "coda4microbiome")
#'
#' x<-impute_zeros(x_HIV)
#
#'
#-------------------------------------------------------------------------------
impute_zeros<-function(x){
  if (min(as.numeric(unlist(x)))< 0) {
    stop("Negative abundance values (check your data)")
  } else {
    if (sum(x==0)>0){
      #xmin = min(x[x > 0]);
      xmin = min(as.numeric(unlist(x))[as.numeric(unlist(x))>0])
      if (xmin >= 1){
        x = x + 1;
      }else{
        x = x+xmin/2;
      }
    }
  }
  return(x)
}
#-------------------------------------------------------------------------------

#' logratios_matrix
#'
#' Computes a large matrix with all the log-ratios between pairs of taxa (columns) in the abundance table
#'
#' @param x abundance table  (rows are samples, columns are variables (taxa))
#'
#' @return list with matrix of log-ratios, matrix indicating the pairs of variables involved in each log-ratio,
#' and a matrix indicating the names of the variables involved in each log-ratio.
#'
#' @export
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(HIV, package = "coda4microbiome")
#'
#' lrHIV<-logratios_matrix(x_HIV[,(1:4)])
#'
#' # matrix of log-ratios (first 6 rows and 6 columns):
#'
#' lrHIV[[1]][1:6,1:6]
#'
#' # variables involved in each log-ratio
#'
#' head(lrHIV[[2]])
#'
#' # names of the variables involved in each log-ratio
#'
#' head(lrHIV[[3]])
#'

#-------------------------------------------------------------------------------
logratios_matrix<-function(x){

  x<-impute_zeros(x)

  if(is.null(colnames(x))) colnames(x)<-(1:ncol(x))

  k<-ncol(x)
  m<-nrow(x)

  logx <- log(x)
  lrcolnames<-NULL
  lrX <- matrix(0,m,k*(k-1)/2)
  idlrX <- matrix(0,k*(k-1)/2,2)
  nameslrX <- matrix(0,k*(k-1)/2,2)
  colnamesx <- colnames(x)
  lloc <-0
  for(i in (1:(k-1))){
    for(j in ((i+1):k)) {
      lloc=lloc+1
      idlrX[lloc,]<-c(i,j)
      nameslrX[lloc,] <- c(colnamesx[i],colnamesx[j])
      lrX[,lloc] <- logx[,i]-logx[,j]
      lrcolnames<-c(lrcolnames,paste(paste("lr",i,sep=""),j,sep="."))
    }
  }
  colnames(lrX)<-lrcolnames
  results <- list(
    "matrix of log-ratios" = lrX,
    "pairs of variables in the logratio" = idlrX,
    "names of the variables in the logratio" = nameslrX
  )

  return(results)
}
#-------------------------------------------------------------------------------

#' explore_logratios
#'
#' Explores the association of each log-ratio with the outcome.
#' Summarizes the importance of each variable (taxa) as the aggregation of
#' the association measures of those log-ratios involving the variable. The output includes a plot
#' of the association of the log-ratio with the outcome where the variables (taxa) are ranked by importance
#'
#' @param x abundance table (rows are samples, columns are variables (taxa))
#' @param y outcome (binary or continuous)
#' @param decreasing order of importance (default = TRUE)
#' @param measure association measures "AUC", "Pearson","Spearman", "glm" (default = "AUC")
#' @param covar data frame with covariates (default = NULL)
#' @param shownames logical, if TRUE, shows the names of the variables in the rows of the plot (default = FALSE)
#' @param maxrow maximum number of rows to display in the plot (default = 15)
#' @param maxcol maximum number of columns to display in the plot (default = 15)
#' @param showtitle  logical, if TRUE, shows the title of the plot (default = TRUE)
#' @param mar mar numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot (default mar=c(0,0,1,0))
#'
#' @return list with "max log-ratio","names max log-ratio", "order of importance",
#' "name of most important variables", "association log-ratio with y" and
#' "top log-ratios plot"
#' @export
#'
#' @importFrom pROC auc
#' @importFrom pROC roc
#' @import glmnet
#' @import stats
#' @importFrom grDevices colorRampPalette
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(HIV, package = "coda4microbiome")
#'
#' explore_logratios(x_HIV,y_HIV)
#'
#-------------------------------------------------------------------------------
explore_logratios<-function(x,y,decreasing=TRUE, measure="AUC", covar=NULL, shownames=FALSE,maxrow=15, maxcol=15, showtitle=TRUE, mar=c(0,0,1,0)){
  # library(corrplot)
  # library(pROC)
  # library(glmnet)

  title<-NULL
  GnBu8 <- c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#08589e")
  Blues4 <- c("#eff3ff", "#bdd7e7", "#6baed6", "#2171b5")

  col2<-grDevices::colorRampPalette(GnBu8, space = "Lab")

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  k<-ncol(x)
  if (maxrow>k) maxrow<-k
  if (maxcol>k) maxcol<-k

  x<-impute_zeros(x)

  lrmatrix<-logratios_matrix(x)
  lrX<-lrmatrix[[1]]
  idlrX<-lrmatrix[[2]]

  logratio_cor<-matrix(rep(0,k*k),nrow=k)

  if (measure=="Pearson"){
    if (showtitle==TRUE){
      title<-"Pearson correlation between y and log(xi/xj)"
    }
    s=0
    for(i in (1:(k-1))){
      for(j in ((i+1):k)){
        s<-s+1
        lr_ij<-lrX[,s]
        logratio_cor[i,j]<-cor(y,lr_ij)
        logratio_cor[j,i]<- (-cor(y,lr_ij))
      }
    }
  }

  if (measure=="Spearman"){
    if (showtitle==TRUE){
      title<-"Spearman correlation between y and log(xi/xj)"
    }

    s=0
    for(i in (1:(k-1))){
      for(j in ((i+1):k)){
        s<-s+1
        lr_ij<-lrX[,s]
        logratio_cor[i,j]<-cor(y,lr_ij, method = "spearman")
        logratio_cor[j,i]<- (-logratio_cor[i,j])
      }
    }
  }
  if (measure=="AUC"){
    if (showtitle==TRUE){
      title<-"AUC: log(xi/xj) discrimination accuracy of y"
    }

    s=0
    for(i in (1:(k-1))){
      for(j in ((i+1):k)){
        s<-s+1
        lr_ij<-lrX[,s]
        if (requireNamespace("pROC", quietly = TRUE)) {
          logratio_cor[i,j]<-pROC::auc(pROC::roc(y, lr_ij,quiet = TRUE))[[1]]
        } else {
          stop("pROC package not loaded")
        }

        logratio_cor[j,i]<- logratio_cor[i,j]
      }
    }
  }

  if (measure=="glm"){

    if (y.binary==TRUE){
      if (showtitle==TRUE){
        title<-"AUC logistic regression y~log(xi/xj)"
      }
      s=0
      for(i in (1:(k-1))){
        for(j in ((i+1):k)){
          s<-s+1
          lr_ij<-lrX[,s]
          if (is.null(covar)){
            model1<-glm(y~lr_ij, family = "binomial")
          } else {
            df<-data.frame(y,lr_ij, covar)
            model1<-glm(y~., family = "binomial", data=df)
          }
          if (requireNamespace("pROC", quietly = TRUE)) {
            logratio_cor[i,j]<-pROC::auc(pROC::roc(y, as.numeric(predict(model1)),quiet = TRUE))[[1]]
          } else {
            stop("pROC package not loaded")
          }

          logratio_cor[j,i]<- logratio_cor[i,j]

        }
      }
    } else {
      if (showtitle==TRUE){
        title<-"Spearman cor linear regression y~log(xi/xj)"
      }

      s=0
      for(i in (1:(k-1))){
        for(j in ((i+1):k)){
          s<-s+1
          lr_ij<-lrX[,s]
          if (is.null(covar)){
            model1<-lm(y~lr_ij)
          } else {
            df<-data.frame(y,lr_ij, covar)
            model1<-lm(y~., data=df)
          }
          beta1<-model1$coefficients[[2]]
          z<-sign(beta1)*predict(model1)
          logratio_cor[i,j]<- cor(y, z,method="spearman")
          logratio_cor[j,i]<- - logratio_cor[i,j]
        }
      }
    }
  }
  o<-(1:k)

  if (decreasing ==T){
    o<-order(colSums(abs(logratio_cor)), decreasing = T)
  }

  M<-logratio_cor[o,o]
  colnames(M)<-o
  rownames(M)<-colnames(M)
  if (shownames==TRUE) {
    rownames(M)<-colnames(x)[o]
  }


  if (y.binary==TRUE){
    if (requireNamespace("corrplot", quietly = TRUE)) {
      (top_lr_plot<- corrplot::corrplot(M[(1:maxrow),(1:maxcol)],tl.pos = "lt",title=title, mar=mar, cl.lim=c(min(M),max(M)), col = col2(200),is.corr = FALSE))
    } else {
      stop("corrplot package not loaded")
    }
  } else {
    if (requireNamespace("corrplot", quietly = TRUE)) {
      (top_lr_plot<- corrplot::corrplot(M[(1:maxrow),(1:maxcol)],tl.pos = "lt",title=title, mar=mar, cl.lim=c(min(M),max(M)),is.corr = FALSE))
    } else {
      stop("corrplot package not loaded")
    }


  }


  results <- list(
    "max log-ratio" = colnames(M)[which(M == max(abs(M)), arr.ind = TRUE)[(2:1)]],
    "names max log-ratio" = colnames(x)[as.numeric(colnames(M)[which(M == max(abs(M)), arr.ind = TRUE)[(2:1)]])],
    "order of importance" = o,
    "name of most important variables" = colnames(x)[o[1:maxrow]],
    "association log-ratio with y"=M,
    "top log-ratios plot"=top_lr_plot)

  return(results)

}
#-------------------------------------------------------------------------------


#' coda_glmnet
#'
#' Microbial signatures in cross-sectional studies.
#' The algorithm performs variable selection through penalized regression on the set of all pairwise log-ratios.
#' The result is expressed as the (weighted) balance between two groups of taxa.
#' It allows the use of non-compositional covariates.
#'
#' @param x abundance table (rows are samples, columns are variables (taxa))
#' @param y outcome (binary or continuous)
#' @param covar data frame with covariates (default = NULL)
#' @param lambda penalization parameter (default = "lambda.1se")
#' @param nvar number of variables to use in the glmnet.fit function (default = NULL)
#' @param alpha elastic net parameter (default = 0.9)
#' @param nfolds number of folds
#' @param showPlots if TRUE, shows the plots (default = TRUE)
#'
#' @return if y is binary: list with "taxa.num","taxa.name","log-contrast coefficients","predictions","apparent AUC","mean cv-AUC","sd cv-AUC","predictions plot","signature plot"
#' if not:list with "taxa.num","taxa.name","log-contrast coefficients","predictions","apparent Rsq","mean cv-MSE","sd cv-MSE","predictions plot","signature plot"
#' @export
#'
#' @importFrom pROC auc
#' @importFrom pROC roc
#' @import glmnet
#' @import stats
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(Crohn, package = "coda4microbiome")
#'
#' set.seed(123)
#'
#' coda_glmnet(x_Crohn[,(1:10)],y_Crohn,showPlots= FALSE)
#'
#'
#-------------------------------------------------------------------------------
coda_glmnet<-function(x,y, covar=NULL, lambda="lambda.1se",nvar=NULL,alpha=0.9, nfolds=10,showPlots= TRUE){

  # library(glmnet)
  # library(pROC)
  # library(ggpubr)

  x<-impute_zeros(x)

  kselect<-ncol(x)

  taxaselect<-(1:ncol(x))

  lrmatrix<-logratios_matrix(x)

  lrX<-lrmatrix[[1]]
  idlrX<-lrmatrix[[2]]
  nameslrX<-lrmatrix[[3]]


  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)


  if (y.binary==TRUE){
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrX,y, family = "binomial" , alpha=alpha, type.measure = "auc", nfolds = nfolds, keep=TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-glm(y~., family = "binomial", data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrX,y, family = "binomial" , offset=x0,alpha=alpha, type.measure = "auc", nfolds = nfolds, keep=TRUE)
    }
  } else {
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrX,y , alpha=alpha, type.measure = "deviance", nfolds = nfolds,keep=TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-lm(y~., data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrX,y , offset=x0, alpha=alpha, type.measure = "deviance", nfolds = nfolds,keep=TRUE)
    }
  }

  if (showPlots==TRUE){
    plot(lassocv)
  }

  if (!is.null(nvar)){
    rowlasso<-max(which(lassocv$glmnet.fit$df<=nvar))  # rowlasso= row in glmnet.fit corresponding to the specified number of variables (nvar)
    lambda<-lassocv$glmnet.fit$lambda[rowlasso]
  }

  lambdavalue<-lambda
  if (is.character(lambda)){
    if (lambda=="lambda.1se") lambdavalue <-lassocv$lambda.1se
    if (lambda=="lambda.min") lambdavalue <-lassocv$lambda.min
  }
  idrow<-max(which(lassocv$glmnet.fit$lambda>=lambdavalue))  # idrow= row in glmnet.fit object corresponding to the specified lambda

  coeflr<-as.vector(coef(lassocv, s = lambda))[-1]
  lrselect<-which(coeflr!=0)


  coeflogcontrast<-rep(0,ncol(x))
  for (i in (1:length(coeflr))){
    coeflogcontrast[idlrX[i,1]]<-coeflogcontrast[idlrX[i,1]]+coeflr[i]
    coeflogcontrast[idlrX[i,2]]<-coeflogcontrast[idlrX[i,2]]-coeflr[i]
  }

  coeflogcontrast<-2*coeflogcontrast/sum(abs(coeflogcontrast))
  varlogcontrast<-which(abs(coeflogcontrast)>0)
  coeflogcontrast<-coeflogcontrast[varlogcontrast]

  (names.select<-colnames(x)[varlogcontrast])

  (positive<-ifelse(coeflogcontrast>0,1,0))

  positive<-factor(positive, levels = c(0,1), labels = c("negative","positive"))

  logcontrast=as.matrix(log(x)[,varlogcontrast])%*%coeflogcontrast
  # logcontrast<-logcontrast-mean(logcontrast)

  if (is.null(covar)){
    predictions<-logcontrast
  } else {
    if (y.binary==TRUE){
      df1<-data.frame(y,logcontrast, covar)
      m1<-glm(y~., family = "binomial", data=df1)
      predictions<-predict(m1)

    } else {
      df1<-data.frame(y,logcontrast, covar)
      m1<-lm(y~., data=df1)
      predictions<-predict(m1)
    }

    # predictions<-predictions-mean(predictions)

  }

  if (y.binary==TRUE){
    AUC_signature<-pROC::auc(pROC::roc(y, as.numeric(predictions),quiet = TRUE))[[1]]
    if (length(varlogcontrast)==0) AUC_signature<- 0.5
    mcvAUC<-lassocv$cvm[idrow]
    sdcvAUC<-lassocv$cvsd[idrow]

  } else {
    mcvMSE<-lassocv$cvm[idrow]
    sdcvMSE<-lassocv$cvsd[idrow]
    Rsq<-as.numeric(cor(predictions,y)^2)
    if (length(varlogcontrast)==0) Rsq <- 0
  }

  plot1 <- plot_prediction(predictions,y,showPlots=showPlots)


  plot2<-plot_signature(names.select,coeflogcontrast,showPlots=showPlots)

  if (y.binary==TRUE){
    results <- list(
      "taxa.num" = varlogcontrast,
      "taxa.name" = names.select,
      "log-contrast coefficients" = coeflogcontrast,
      "predictions"=predictions,
      "apparent AUC"= AUC_signature,
      "mean cv-AUC"= mcvAUC,
      "sd cv-AUC"= sdcvAUC,
      "predictions plot"=plot1,
      "signature plot"=plot2)
  } else {
    results <- list(
      "taxa.num" = varlogcontrast,
      "taxa.name" = names.select,
      "log-contrast coefficients" = coeflogcontrast,
      "predictions"=predictions,
      "apparent Rsq" = Rsq,
      "mean cv-MSE"= mcvMSE,
      "sd cv-MSE"= sdcvMSE,
      "predictions plot"=plot1,
      "signature plot"=plot2)
  }
  return(results)
}

#' plot_prediction
#'
#' Plot of the predictions of a fitted model:
#' Multiple box-plot and density plots for binary outcomes and Regression plot for continuous outcome
#'
#' @param prediction the fitted values of predictions for the model
#' @param y outcome (binary or continuous)
#' @param strata stratification variable (default = NULL)
#' @param showPlots if TRUE, shows the plots (default = TRUE)
#'
#' @return prediction plot
#' @export
#'
#' @import ggpubr
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' # prediction plot for the log-ratio between columns 3 and 32 on HIV status
#'
#' data(HIV, package = "coda4microbiome")
#'
#' x<-impute_zeros(x_HIV)
#'
#' lr<-log(x[,3])-log(x[,32])
#'
#' plot_prediction(lr, y_HIV)
#'
#'
#-------------------------------------------------------------------------------
plot_prediction <- function(prediction, y, strata=NULL, showPlots=TRUE){

  # library(ggpubr)

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  if (is.null(strata)){
    data <- data.frame(prediction,y)
    nameX <- as.character(substitute(prediction))
    nameY <- as.character(substitute(y))
    colnames(data)=c(nameX, nameY)
  } else {
    data <- data.frame(prediction,y,strata)
    nameX <- as.character(substitute(prediction))
    nameY <- as.character(substitute(y))
    nameZ <- as.character(substitute(strata))
    colnames(data)=c(nameX, nameY,nameZ)
  }

  # colores <- c("gold1","orchid")
  colores <- c("chocolate1","slateblue2")
  if (y.binary == TRUE){

    pdens <- ggpubr::ggdensity(data, x=nameX,
                       add = "mean", rug = FALSE,
                       color = nameY, fill = nameY,
                       palette = c(colores[1], colores[2])
                       # palette ="PuOr"
    )
    if (is.null(strata)){
      pbbox <- ggpubr::ggboxplot(data, x = nameY, y = nameX,
                         color = nameY,
                         palette = c(colores[1], colores[2])
                         # palette ="PuOr"
      )+rotate()
    } else {
      pbbox <- ggpubr::ggboxplot(data, x = nameZ, y = nameX,
                         color = nameY,
                         palette = c(colores[1], colores[2])
                         # palette ="PuOr"
                         , shape=nameZ
      )+rotate()
    }

    L <- ggpubr::ggarrange( pbbox, pdens,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)

  }else{
    if (is.null(strata)){
      pscat<- ggpubr::ggscatter(data, x = nameX, y = nameY,
                        add = "reg.line",                                      # regression line
                        color = "deepskyblue3",
                        conf.int = TRUE,                                       # confidence interval
                        add.params = list(color = "coral1", fill = "bisque")   # color
      )+
        stat_cor(method = "pearson")
      # stat_cor(method = "pearson", label.x = .1, label.y = 0.8)      # correlation coefficient
    } else {
      pscat<- ggpubr::ggscatter(data, x = nameX, y = nameY,
                        add = "reg.line",                                      # regression line
                        conf.int = TRUE,                                       # confidence interval
                        color = nameZ,
                        facet.by = nameZ,
                        add.params = list(color = "coral1", fill = "bisque")   # color
      )+
        stat_cor(method = "pearson")
      # stat_cor(method = "pearson", label.x = .1, label.y = 0.8)      # correlation coefficient
    }
    L <- pscat
  }

  if (showPlots==TRUE){
    print(L)
  }

  return(L)

}

#' plot_signature
#'
#' Graphical representation of the variables selected and their coefficients
#'
#' @param vars variables selected
#' @param coeff associated coefficients
#' @param showPlots if TRUE, shows the plots (default = TRUE)
#' @param varnames if TRUE, shows the names of the variables
#'
#' @return bar plot
#'
#' @export
#'
#' @import ggpubr
#' @import ggplot2
#'
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' plot_signature(c(2,10, 3, 15, 4), c(0.8, -0.1, 0.2, -0.6, -0.3))
#'
#'
#-------------------------------------------------------------------------------
plot_signature <- function(vars, coeff, showPlots=TRUE, varnames=NULL){

  # library(ggpubr)


  positive<-ifelse(coeff>0,1,0)

  positive<-factor(positive, levels = c(0,1), labels = c("negative","positive"))


  df<-data.frame(vars,coeff=round(coeff,digits = 2), positive)

  if (!is.null(varnames)){
    df$vars<-varnames
  }


  L<-ggpubr::ggbarplot(df, x = "vars", y = "coeff", color="positive",fill="positive",
               sort.val="asc",  orientation = "horiz",position = ggplot2::position_dodge(),
               label = TRUE, lab.vjust=0.2, lab.hjust=0.5,ylab=FALSE)

  if (showPlots==TRUE){

    print(L)
  }

  return(L)

}

#' coda_glmnet0
#'
#' Internal function for the permutational test
#'
#' @param x .
#' @param lrX .
#' @param idlrX .
#' @param nameslrX .
#' @param y .
#' @param covar .
#' @param lambda .
#' @param alpha .
#'
#' @return .
#'
#' @importFrom pROC auc
#' @importFrom pROC roc
#' @import glmnet
#' @import stats
#'
#' @author M. Calle - T. Susin
#'
#'
#-------------------------------------------------------------------------------
coda_glmnet0<-function(x,lrX,idlrX,nameslrX,y, covar=NULL, lambda="lambda.1se",alpha=0.9){
  #suppressWarnings()

  # library(glmnet)
  # library(pROC)
  # library(ggpubr)


  if (sum(x==0)>0){
    x<-impute_zeros(x)
  }


  kselect<-ncol(x)

  idlrXsub<-idlrX
  lrXsub<-lrX

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  if (y.binary==TRUE){
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y, family = "binomial" , alpha=alpha, type.measure = "auc", keep=TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-glm(y~., family = "binomial", data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y, family = "binomial" , offset=x0,alpha=alpha, type.measure = "auc", keep=TRUE)
    }
  } else {
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y , alpha=alpha, type.measure = "deviance", keep=TRUE)
    } else {
      df0<-data.frame(y,covar)
      model0<-lm(y~., data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y , offset=x0, alpha=alpha, type.measure = "deviance", keep=TRUE)
    }
  }

  lambdavalue<-lambda
  if (is.character(lambda)){
    if (lambda=="lambda.1se") lambdavalue <-lassocv$lambda.1se
    if (lambda=="lambda.min") lambdavalue <-lassocv$lambda.min
  }
  idrow<-max(which(lassocv$glmnet.fit$lambda>=lambdavalue))  # idrow= row in glmnet.fit object corresponding to the specified lambda

  coeflr<-as.vector(coef(lassocv, s = lambda))[-1]
  lrselect<-which(coeflr!=0)

  idlrXsub[lrselect,]

  coeflogcontrast<-rep(0,ncol(x))
  for (i in (1:length(coeflr))){
    coeflogcontrast[idlrXsub[i,1]]<-coeflogcontrast[idlrXsub[i,1]]+coeflr[i]
    coeflogcontrast[idlrXsub[i,2]]<-coeflogcontrast[idlrXsub[i,2]]-coeflr[i]
  }

  coeflogcontrast<-2*coeflogcontrast/sum(abs(coeflogcontrast))
  varlogcontrast<-which(abs(coeflogcontrast)>0)
  coeflogcontrast<-coeflogcontrast[varlogcontrast]

  (names.select<-colnames(x)[varlogcontrast])

  (positive<-ifelse(coeflogcontrast>0,1,0))

  positive<-factor(positive, levels = c(0,1), labels = c("negative","positive"))

  #df<-data.frame(taxa.name=names.select, taxa.num=varlogcontrast, coefficient=round(coeflogcontrast,digits = 2), positive)

  logcontrast=as.matrix(log(x)[,varlogcontrast])%*%coeflogcontrast
  # logcontrast<-logcontrast-mean(logcontrast)

  if (is.null(covar)){
    predictions<-logcontrast
  } else {
    if (y.binary==TRUE){
      df1<-data.frame(y,logcontrast, covar)
      m1<-glm(y~., family = "binomial", data=df1)
      predictions<-predict(m1)

    } else {
      df1<-data.frame(y,logcontrast, covar)
      m1<-lm(y~., data=df1)
      predictions<-predict(m1)
    }

    # predictions<-predictions-mean(predictions)

  }

  if (y.binary==TRUE){
    AUC_signature<-pROC::auc(pROC::roc(y, as.numeric(predictions),quiet = TRUE))[[1]]
    if (length(varlogcontrast)==0) AUC_signature<- 0.5
    mcvAUC<-lassocv$cvm[idrow]
    sdcvAUC<-lassocv$cvsd[idrow]

  } else {
    mcvMSE<-lassocv$cvm[idrow]
    sdcvMSE<-lassocv$cvsd[idrow]
    Rsq<-as.numeric(cor(predictions,y)^2)
    if (length(varlogcontrast)==0) Rsq <- 0
  }

  if (y.binary==TRUE){
    results <- list(
      "taxa.num" = varlogcontrast,
      "taxa.name" = names.select,
      "log-contrast coefficients" = coeflogcontrast,
      "predictions"=predictions,
      "apparent AUC"= AUC_signature,
      "mean cv-AUC"= mcvAUC,
      "sd cv-AUC"= sdcvAUC)
  } else {
    results <- list(
      "taxa.num" = varlogcontrast,
      "taxa.name" = names.select,
      "log-contrast coefficients" = coeflogcontrast,
      "predictions"=predictions,
      "apparent Rsq" = Rsq,
      "mean cv-MSE"= mcvMSE,
      "sd cv-MSE"= sdcvMSE)
  }
  return(results)
}



#' coda_glmnet_null
#'
#' Performs a permutational test for the coda_glmnet() algorithm:
#' It provides the distribution of results under the null hypothesis by
#' implementing the coda_glmnet() on different rearrangements of the response variable.
#'
#' @param x abundance table (rows are samples, columns are variables (taxa))
#' @param y outcome (binary or continuous)
#' @param niter number of iterations (default = 100)
#' @param covar data frame with covariates (default = NULL)
#' @param lambda penalization parameter (default = "lambda.1se")
#' @param alpha elastic net parameter (default = 0.9)
#' @param sig significance level for the confidence interval (default = 0.05)
#'
#' @return a list with "accuracy" and "confidence interval"
#' @export
#'
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(Crohn, package = "coda4microbiome")
#'
#' coda_glmnet_null(x=x_Crohn[,(1:10)], y=y_Crohn, niter=2,covar=NULL,lambda="lambda.1se",
#'                                                 alpha=0.9,sig=0.05)
#'
#'
#-------------------------------------------------------------------------------
coda_glmnet_null<-function(x,y,niter=100,covar=NULL,lambda="lambda.1se", alpha=0.9,sig=0.05){

  alpha0<-alpha
  lambda0<-lambda
  covar0<-covar


  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)
  y1<-y
  lrmatrix<-logratios_matrix(x)

  lrX<-lrmatrix[[1]]
  idlrX<-lrmatrix[[2]]
  nameslrX<-lrmatrix[[3]]

  accuracy<-rep(0,niter)
  for(i in (1:niter)){
    y1<-sample(y1)
    lr<-coda_glmnet0(x=x,lrX=lrX,idlrX=idlrX,nameslrX=nameslrX,y=y1,lambda=lambda0,covar=covar0, alpha=alpha0)
    if (y.binary==TRUE){
      res<-lr$`mean cv-AUC`
    } else {
      res<-lr$`mean cv-MSE`
    }
    accuracy[i]<-res
    print(paste("iter",i))
  }
  results <- list(
    "accuracy"=accuracy,
    "confidence interval"=quantile(accuracy, c((sig/2),(1-(sig/2))))
  )
  return(results)
}


# Shannon effective number of species
#----------------------------------------------------

#' shannon
#'
#' Shannon information
#'
#' @param x abundance composition (vector)
#'
#' @return shannon information
#' @export
#'
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(HIV, package = "coda4microbiome")
#'
#' shannon(x_HIV[1,])
#'
#'
shannon<-function(x){
  p<-x/sum(x)
  plog<-as.numeric(ifelse(p>0,p*log(p),0))
  s<- -sum(plog)
  s
}

#' shannon_effnum
#'
#' Shannon effective number of variables in a composition
#'
#' @param x abundance composition (vector)
#'
#' @return shannon information
#' @export
#'
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(HIV, package = "coda4microbiome")
#'
#' shannon_effnum(x_HIV[1,])
#'
#-------------------------------------------------------
shannon_effnum<-function(x){
  p<-x/sum(x)
  plog<-as.numeric(ifelse(p>0,p*log(p),0))
  s<- -sum(plog)
  exp(s)
}

#' shannon_sim
#'
#' Shannon similarity between two compositions
#'
#' @param x abundance composition (vector)
#'
#' @param y abundance composition (vector)
#'
#' @return shannon similarity (value between 0 and 1)
#' @export
#'
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(HIV, package = "coda4microbiome")
#'
#' shannon_sim(x_HIV[1,],x_HIV[2,])
#'
#-------------------------------------------------------

shannon_sim<-function(x,y){
  wx<-sum(x)/(sum(x)+sum(y))
  wy<-sum(y)/(sum(x)+sum(y))
  px<-x/sum(x)
  py<-y/sum(y)
  ha<-(wx*shannon(x)+wy*shannon(y))
  hg<-shannon(wx*px+wy*py)
  sim<-2*(exp(ha)/exp(hg)-1/2)
  sim
}

