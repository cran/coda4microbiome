#' integralFun
#'
#' Integral of the curve trajectory of subject_id in the interval a,b
#'
#' @param x abundance matrix or data frame in long format (several rows per individual)
#' @param y outcome (binary); data type: numeric, character or factor vector
#' @param id subjects-ids
#' @param a interval initial time
#' @param b interval final time
#'
#' @return matrix with integrals for each individual (rows) and each taxa (columns)
#'
#'
#' @author M. Calle - T. Susin
#'
#'
#-------------------------------------------------------------------------------
integralFun <- function(x, y, id, a, b){
  df<-data.frame(x,y,id)
  inputa<-mean(df$y[abs(df$x-a)<0.1*(b-a)])
  inputb<-mean(df$y[abs(df$x-b)<0.1*(b-a)])
  subjects=unique(df$id);
  n=length(subjects); #num of subjects
  integral=NULL;
  j=0;
  for (k in subjects){
    j=j+1;
    curve <- subset(df, id==k)
    xref<-curve$x
    ycurve<-curve$y
    if (min(xref)>a){
      xref<-c(a,xref)
      ycurve<-c(inputa,ycurve)
    }
    if (max(xref)<b){
      xref<-c(xref,b)
      ycurve<-c(ycurve,inputb)
    }
    integral[j]=integrate(approxfun(xref,ycurve),a,b,stop.on.error = FALSE)$value
  }

  return(integral)
}

#' explore_lr_longitudinal
#'
#' Explores the association of summary (integral) of each log-ratio trajectory with the outcome.
#' Summarizes the importance of each variable (taxa) as the aggregation of
#' the association measures of those log-ratios involving the variable. The output includes a plot
#' of the association of the log-ratio with the outcome where the variables (taxa) are ranked by importance
#'
#' @param x abundance matrix or data frame in long format (several rows per individual)
#' @param y outcome (binary); data type: numeric, character or factor vector
#' @param x_time observation times
#' @param subject_id subject id
#' @param ini_time initial time to be analyzed
#' @param end_time end time to be analyzed
#' @param showPlots if TRUE, shows the plot (default = FALSE)
#' @param decreasing order of importance (default = TRUE)
#' @param covar data frame with covariates (default = NULL)
#' @param shownames if TRUE, shows the names of the variables in the rows of the plot (default = FALSE)
#' @param maxrow maximum number of rows to display in the plot (default = 15)
#' @param maxcol maximum number of columns to display in the plot (default = 15)
#' @param showtitle  logical, if TRUE, shows the title of the plot (default = TRUE)
#' @param mar mar numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot (default mar=c(0,0,1,0))
#'
#' @return list with "max log-ratio","names max log-ratio","order of importance","name of most important variables","association log-ratio with y","top log-ratios plot"
#'
#' @export
#'
#' @importFrom corrplot corrplot
#' @importFrom pROC auc
#' @importFrom pROC roc
#' @import stats
#' @importFrom grDevices colorRampPalette
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' set.seed(123) # to reproduce the results
#'
#' data(ecam_filtered, package = "coda4microbiome")   # load the data
#'
#' x=x_ecam # microbiome abundance
#' x_time = metadata$day_of_life    # observation times
#' subject_id = metadata$studyid   # subject id
#' y= metadata$diet           # diet ("bd"= breast diet, "fd"=formula diet)
#' ini_time = 0
#' end_time = 90
#'
#' ecam_logratios<-explore_lr_longitudinal(x,y,x_time,subject_id,ini_time,end_time)
#'
#-------------------------------------------------------------------------------
explore_lr_longitudinal<-function(x,y, x_time, subject_id, ini_time, end_time,
                                  showPlots=FALSE, decreasing=TRUE, covar=NULL,
                                  shownames=FALSE,maxrow=15, maxcol=15,
                                  showtitle=TRUE, mar=c(0,0,1,0)){

  # library(corrplot)
  # library(pROC)

  title<-NULL
  X1<-impute_zeros(x)


  # Data without zeros......
  logX1 = log(X1);
  subject_id<-as.numeric(as.factor(subject_id))
  nsubjects=length(unique(subject_id)); #num of subjects

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  if (y.binary==T) y=factor(y)

  indexUser=seq_along(subject_id)[!duplicated(subject_id)];
  y_unique<-as.numeric(y[indexUser])-1

  if (!is.null(covar)){
    covar=covar[indexUser,]
  }

  # Compute all the column integrals
  intLogX <-NULL
  for (ki in (1:(ncol(logX1)))){
    print(paste('ind=', ki))
    yy=as.numeric(logX1[,ki]);
    integrals=integralFun(x_time, yy, subject_id, a=ini_time, b=end_time)
    intLogX<-cbind(intLogX, matrix(integrals));
  }

  GnBu8 <- c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#08589e")
  Blues4 <- c("#eff3ff", "#bdd7e7", "#6baed6", "#2171b5")

  col2<-grDevices::colorRampPalette(GnBu8, space = "Lab")



  k<-ncol(x)
  if (maxrow>k) maxrow<-k

  if (maxcol>k) maxcol<-k


  m<-length(y_unique)
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
      lrX[,lloc] <- intLogX[,i]-intLogX[,j]
      lrcolnames<-c(lrcolnames,paste(paste("lr",i,sep=""),j,sep="."))
    }
  }
  colnames(lrX)<-lrcolnames


  logratio_cor<-matrix(rep(0,k*k),nrow=k)




  if (y.binary==TRUE){
    y_unique<-factor(y_unique)
    if (showtitle==TRUE){
          title<-"AUC logistic regression y~Integral(log(xi/xj))"
    }
    s=0
    for(i in (1:(k-1))){
      for(j in ((i+1):k)){
        s<-s+1
        lr_ij<-lrX[,s]
        if (is.null(covar)){
          model1<-glm(y_unique~lr_ij, family = "binomial")
        } else {
          df<-data.frame(y_unique,lr_ij, covar)
          model1<-glm(y_unique~., family = "binomial", data=df)
        }
        logratio_cor[i,j]<-pROC::auc(pROC::roc(y_unique, predict(model1),quiet = TRUE))[[1]]
        logratio_cor[j,i]<- logratio_cor[i,j]

      }
    }
  } else {
    if (showtitle==TRUE){
          title<-"Spearman cor linear regression y~Integral(log(xi/xj))"
    }

    s=0
    for(i in (1:(k-1))){
      for(j in ((i+1):k)){
        s<-s+1
        lr_ij<-lrX[,s]
        if (is.null(covar)){
          model1<-lm(y_unique~lr_ij)
        } else {
          df<-data.frame(y_unique,lr_ij, covar)
          model1<-lm(y_unique~., data=df)
        }
        #logratio_cor[i,j]<- summary(model1)$adj.r.squared
        beta1<-model1$coefficients[[2]]
        z<-sign(beta1)*predict(model1)
        logratio_cor[i,j]<-cor(y_unique, z,method="spearman")
        logratio_cor[j,i]<- - logratio_cor[i,j]
      }
    }
  }

  o<-(1:k)

  if (decreasing ==T){
    o<-order(colSums(abs(logratio_cor)), decreasing = T)
  }


  M<-logratio_cor[o,o]
  colnames(M)<-o
  rownames(M)<-o
  if (shownames==TRUE) {
    rownames(M)<-colnames(x)[o]
  }


  if (y.binary==TRUE){
    if (requireNamespace("corrplot", quietly = TRUE)) {
      (top_lr_plot<-corrplot::corrplot(M[(1:maxrow),(1:maxcol)],tl.pos = "lt",title=title, mar=mar, cl.lim=c(min(M),max(M)), col = col2(200),is.corr = FALSE))
    } else {
      stop("corrplot package not loaded")
    }

  } else {
    if (requireNamespace("corrplot", quietly = TRUE)) {
      (top_lr_plot<-corrplot::corrplot(M[(1:maxrow),(1:maxcol)],tl.pos = "lt",title=title, mar=mar, cl.lim=c(min(M),max(M)),is.corr = FALSE))
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


#' coda_glmnet_longitudinal
#'
#' Microbial signatures in longitudinal studies.
#' Identification of a set of microbial taxa whose joint dynamics is associated with the phenotype of interest.
#' The algorithm performs variable selection through penalized regression over the summary of the log-ratio trajectories (AUC).
#' The result is expressed as the (weighted) balance between two groups of taxa.
#'
#' @param x abundance matrix or data frame in long format (several rows per individual)
#' @param y outcome (binary); data type: numeric, character or factor vector
#' @param x_time observation times
#' @param subject_id subject id
#' @param ini_time initial time to be analyzed
#' @param end_time end time to be analyzed
#' @param covar data frame with covariates (default = NULL)
#' @param lambda penalization parameter (default = "lambda.1se")
#' @param nvar number of variables to use in the glmnet.fit function (default = NULL)
#' @param alpha elastic net parameter (default = 0.9)
#' @param nfolds number of folds (default = 10)
#' @param showPlots if TRUE, shows the plots (default = FALSE)
#' @param coef_threshold coefficient threshold, minimum absolute value of the coefficient for a variable to be included in the model (default =0)
#'
#'
#' @return in case of binary outcome: list with "taxa.num","taxa.name","log-contrast coefficients","predictions","apparent AUC","mean cv-AUC","sd cv-AUC","predictions plot","signature plot","trajectories plot"
#'
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
#' data(ecam_filtered, package = "coda4microbiome")   # load the data
#'
#' ecam_results<-coda_glmnet_longitudinal (x=x_ecam[,(1:4)],y= metadata$diet,
#' x_time= metadata$day_of_life, subject_id = metadata$studyid, ini_time=0,
#' end_time=60,lambda="lambda.min",nfolds=4, showPlots=FALSE)
#'
#' ecam_results$taxa.num
#'
#-------------------------------------------------------------------------------
coda_glmnet_longitudinal <- function(x,y, x_time, subject_id, ini_time, end_time,
                                   covar=NULL, lambda="lambda.1se",nvar=NULL,alpha=0.9,nfolds=10, showPlots=TRUE, coef_threshold=0){
  #suppressWarnings()

  # library(glmnet)
  # library(pROC)

  yini<-y

  if (sum(x==0)>0){
    x<-impute_zeros(x)
  }
  logX1 = log(x);
  #subject_id<-as.numeric(as.factor(subject_id))
  #subject_id<-as.numeric(subject_id)

  nsubjects=length(unique(subject_id)); #num of subjects


  indexUser=seq_along(subject_id)[!duplicated(subject_id)];

  y_unique<-y[indexUser]

  if (!is.null(covar)){
    covar=covar[indexUser,]
  }


  y.binary<-ifelse(dim(table(y_unique))==2, TRUE, FALSE)


  if (y.binary==TRUE){
    y<-factor(y)
    y_unique<-as.factor(y_unique)
  }

  if (is.factor(y)){
    labelsy<-levels(y)
    y_unique<-factor(y_unique, labels=labelsy )

  }



  alpha0<-alpha

  kselect<-ncol(x)

  taxaselect<-(1:ncol(x))

  k<-ncol(x)
  # Compute all the column integrals
  intLogX <-NULL
  for (ki in (1:(ncol(logX1)))){
    print(paste('ind=', ki))
    yy=as.numeric(logX1[,ki]);
    integrals=integralFun(x_time, yy, subject_id, a=ini_time, b=end_time)
    intLogX<-cbind(intLogX, matrix(integrals));
  }
  print(dim(intLogX))


  m<-length(y_unique)
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
      lrX[,lloc] <- intLogX[,i]-intLogX[,j]
      lrcolnames<-c(lrcolnames,paste(paste("lr",i,sep=""),j,sep="."))
    }
  }
  colnames(lrX)<-lrcolnames


  idlrXsub<-idlrX
  lrXsub<-lrX

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)


  if (y.binary==TRUE){
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique, family = "binomial", alpha=alpha0, type.measure = "auc", nfolds=nfolds,keep=TRUE)
    } else {
      df0<-data.frame(y_unique,covar)
      model0<-glm(y_unique~., family = "binomial", data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique, family = "binomial" , offset=x0, alpha=alpha0, type.measure = "auc", nfolds=nfolds, keep=TRUE)
    }
  } else {
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique , alpha=alpha0, type.measure = "deviance", nfolds=nfolds, keep=TRUE)
    } else {
      df0<-data.frame(y_unique,covar)
      model0<-lm(y_unique~., data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique , offset=x0, alpha=alpha0, type.measure = "deviance", nfolds=nfolds, keep=TRUE)
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

  idlrXsub[lrselect,]

  coeflogcontrast<-rep(0,ncol(x))
  for (i in (1:length(coeflr))){
    coeflogcontrast[idlrXsub[i,1]]<-coeflogcontrast[idlrXsub[i,1]]+coeflr[i]
    coeflogcontrast[idlrXsub[i,2]]<-coeflogcontrast[idlrXsub[i,2]]-coeflr[i]
  }

  #coeflogcontrast<-2*coeflogcontrast/sum(abs(coeflogcontrast))
  #varlogcontrast<-which(abs(coeflogcontrast)>0)
  varlogcontrast<-which(abs(coeflogcontrast)>coef_threshold)
    coeflogcontrast<-coeflogcontrast[varlogcontrast]

  (names.select<-colnames(x)[varlogcontrast])

  (positive<-ifelse(coeflogcontrast>0,1,0))

  positive<-factor(positive, levels = c(0,1), labels = c("negative","positive"))


  logcontrast=as.matrix(lrXsub[,lrselect])%*%coeflr[lrselect]


  if (is.null(covar)){
    #predictions<-logcontrast
    predictions<-as.numeric(predict(lassocv,lrX,s=lambdavalue))

  } else {
    # if (y.binary==TRUE){
    #   df1<-data.frame(y_unique,logcontrast, covar)
    #   m1<-glm(y_unique~., family = "binomial", data=df1)
    #   predictions<-predict(m1)
    #
    # } else {
    #   df1<-data.frame(y_unique,logcontrast, covar)
    #   m1<-lm(y_unique~., data=df1)
    #   predictions<-predict(m1)
    # }

    #predictions<-x0+logcontrast
    predictions<-as.numeric(predict(lassocv,lrX,s=lambdavalue, newoffset=x0))


  }

  coeflogcontrast<-2*coeflogcontrast/sum(abs(coeflogcontrast))

  if (y.binary==TRUE){
    AUC_signature<-pROC::auc(pROC::roc(y_unique, as.numeric(predictions),quiet = TRUE))[[1]]
    if (length(varlogcontrast)==0) AUC_signature<- 0.5
    mcvAUC<-lassocv$cvm[idrow]
    sdcvAUC<-lassocv$cvsd[idrow]

  } else {
    mcvDev<-lassocv$cvm[idrow]
    sdcvDev<-lassocv$cvsd[idrow]
    Rsq <- 0
    if (length(varlogcontrast)>0){
      Rsq<-as.numeric(cor(predictions,y)^2)
    }
  }

  y<-y_unique

  plot1<-NULL
  plot2<-NULL

  if (length(lrselect>0)){

    plot1 <- plot_prediction(predictions,y, showPlots = showPlots)

    plot2<-plot_signature(names.select,coeflogcontrast, showPlots = showPlots)

  } else {
    print("No variables are selected. The prediction and the signature plots are not displayed.")
  }



  plot3<-NULL

  if (showPlots==TRUE){

  plot3<-plot_signature_curves(varlogcontrast,coeflogcontrast, x=x,y=yini, x_time, subject_id, ini_time, end_time)
  }

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
      "signature plot"=plot2,
      "trajectories plot"=plot3)
  } else {
    results <- list(
      "taxa.num" = varlogcontrast,
      "taxa.name" = names.select,
      "log-contrast coefficients" = coeflogcontrast,
      "predictions"=predictions,
      "apparent Rsq" = Rsq,
      "mean cv-Deviance"= mcvDev,
      "sd cv-Deviance"= sdcvDev,
      "predictions plot"=plot1,
      "signature plot"=plot2,
      "trajectories plot"=plot3)
  }
  return(results)
}


#' coda_glmnet_longitudinal0
#'
#' internal function
#'
#' @param x abundance matrix or data frame in long format (several rows per individual)
#' @param lrX log-ratio matrix
#' @param idlrX indices table in the log-ratio matrix
#' @param nameslrX colnames of the log-ratio matrix
#' @param y outcome (binary); data type: numeric, character or factor vector
#' @param x_time observation times
#' @param subject_id subject id
#' @param ini_time initial time to be analyzed
#' @param end_time end time to be analyzed
#' @param covar data frame with covariates (default = NULL)
#' @param ktop given number of selected taxa or compute the best number in case it is NULL (default = NULL)
#' @param lambda penalization parameter (default = "lambda.1se")
#' @param alpha elastic net parameter (default = 0.9)
#' @param nfolds number of folds
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
coda_glmnet_longitudinal0 <- function(x,lrX,idlrX,nameslrX,y, x_time, subject_id, ini_time, end_time,
                                    covar=NULL, ktop=NULL, lambda="lambda.1se",alpha=0.9,nfolds=10){
  #suppressWarnings()

  # library(glmnet)
  # library(pROC)

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  if (sum(x==0)>0){
    x<-impute_zeros(x)
  }
  logX1 = log(x);
  #subject_id<-as.numeric(as.factor(subject_id))

  nsubjects=length(unique(subject_id)); #num of subjects


  indexUser=seq_along(subject_id)[!duplicated(subject_id)];

  if (!is.null(covar)){
    covar=covar[indexUser,]
  }


  #y_unique<-as.numeric(y[indexUser])-1
  y_unique<-y[indexUser]


  if (y.binary == TRUE){
    y<-factor(y)
  y_unique<-as.factor(y_unique)
  }

  if (is.factor(y)){
    labelsy<-levels(y)
    y_unique<-factor(y_unique, labels=labelsy )

  }


  alpha0<-alpha

  kselect<-ncol(x)

  idlrXsub<-idlrX
  lrXsub<-lrX



  if (y.binary==TRUE){
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique, family = "binomial", alpha=alpha0, type.measure = "auc", nfolds=nfolds,keep=TRUE)
    } else {
      df0<-data.frame(y_unique,covar)
      model0<-glm(y_unique~., family = "binomial", data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique, family = "binomial" , offset=x0, alpha=alpha0, type.measure = "auc", nfolds=nfolds, keep=TRUE)
    }
  } else {
    if(is.null(covar)){
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique , alpha=alpha0, type.measure = "deviance", nfolds=nfolds, keep=TRUE)
    } else {
      df0<-data.frame(y_unique,covar)
      model0<-lm(y_unique~., data=df0)
      x0<-predict(model0)
      lassocv<-glmnet::cv.glmnet(lrXsub,y_unique , offset=x0, alpha=alpha0, type.measure = "deviance", nfolds=nfolds, keep=TRUE)
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

  #coeflogcontrast<-2*coeflogcontrast/sum(abs(coeflogcontrast))
  varlogcontrast<-which(abs(coeflogcontrast)>0)
  coeflogcontrast<-coeflogcontrast[varlogcontrast]

  (names.select<-colnames(x)[varlogcontrast])

  (positive<-ifelse(coeflogcontrast>0,1,0))

  positive<-factor(positive, levels = c(0,1), labels = c("negative","positive"))

  #df<-data.frame(taxa.name=names.select, taxa.num=varlogcontrast, coefficient=round(coeflogcontrast,digits = 2), positive)



  logcontrast=as.matrix(lrXsub[,lrselect])%*%coeflr[lrselect]
  # logcontrast<-logcontrast-mean(logcontrast)


  if (is.null(covar)){
    #predictions<-logcontrast
    predictions<-as.numeric(predict(lassocv,lrX,s=lambdavalue))

  } else {
    # if (y.binary==TRUE){
    #   df1<-data.frame(y_unique,logcontrast, covar)
    #   m1<-glm(y_unique~., family = "binomial", data=df1)
    #   predictions<-predict(m1)
    #
    # } else {
    #   df1<-data.frame(y_unique,logcontrast, covar)
    #   m1<-lm(y_unique~., data=df1)
    #   predictions<-predict(m1)
    # }

    #predictions<-x0+logcontrast
    predictions<-as.numeric(predict(lassocv,lrX,s=lambdavalue, newoffset=x0))

  }

  coeflogcontrast<-2*coeflogcontrast

  if (y.binary==TRUE){
    AUC_signature<-pROC::auc(pROC::roc(y_unique, as.numeric(predictions),quiet = TRUE))[[1]]
    if (length(varlogcontrast)==0) AUC_signature<- 0.5
    mcvAUC<-lassocv$cvm[idrow]
    sdcvAUC<-lassocv$cvsd[idrow]

  } else {
    mcvDev<-lassocv$cvm[idrow]
    sdcvDev<-lassocv$cvsd[idrow]
    Rsq <- 0
    if (length(varlogcontrast)>0){
      Rsq<-as.numeric(cor(predictions,y)^2)
    }
  }

  y<-y_unique

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
      "mean cv-Deviance"= mcvDev,
      "sd cv-Deviance"= sdcvDev)
  }
  return(results)
}




#' coda_glmnet_longitudinal_null
#'
#' Performs a permutational test for the coda_glmnet_longitudinal() algorithm:
#' It provides the distribution of results under the null hypothesis by
#' implementing the coda_glmnet_longitudinal() on different rearrangements of the response variable.

#'
#' @param x abundance matrix or data frame in long format (several rows per individual)
#' @param y outcome (binary); data type: numeric, character or factor vector
#' @param x_time observation times
#' @param subject_id subject id
#' @param ini_time initial time to be analyzed
#' @param end_time end time to be analyzed
#' @param niter number of sample iterations
#' @param covar data frame with covariates (default = NULL)
#' @param alpha elastic net parameter (default = 0.9)
#' @param lambda penalization parameter (default = "lambda.1se")
#' @param nfolds number of folds
#' @param sig significance value (default = 0.05)
#'
#'
#'
#' @return list with "accuracy" and "confidence interval"
#'
#' @export
#'
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' set.seed(123) # to reproduce the results
#'
#' data(ecam_filtered, package = "coda4microbiome")   # load the data
#'
#' x=x_ecam # microbiome abundance
#' x_time = metadata$day_of_life    # observation times
#' subject_id = metadata$studyid   # subject id
#' y= metadata$diet           # diet ("bd"= breast diet, "fd"=formula diet)
#' ini_time = 0
#' end_time = 90
#'
#'  coda_glmnet_longitudinal_null (x,y, x_time, subject_id, ini_time, end_time,
#'                                       lambda="lambda.min",nfolds=4, niter=3)
#'
#'
#-------------------------------------------------------------------------------
coda_glmnet_longitudinal_null<-function(x,y,x_time, subject_id, ini_time, end_time,
                                      niter=100,covar=NULL, alpha=0.9, lambda="lambda.1se", nfolds=10,
                                      sig=0.05){

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  if (sum(x==0)>0){
    x<-impute_zeros(x)
  }
  logX1 = log(x);

  #subject_id<-as.numeric(as.factor(subject_id))
  nsubjects=length(unique(subject_id)); #num of subjects


  indexUser=seq_along(subject_id)[!duplicated(subject_id)];

  if (!is.null(covar)){
    covar=covar[indexUser,]
  }

  #y_unique<-as.numeric(y[indexUser])-1

  y_unique<-y[indexUser]

  if (y.binary==TRUE){
    y<-factor(y)
  y_unique<-as.factor(y_unique)
  }

  if (is.factor(y)){
    labelsy<-levels(y)
    y_unique<-factor(y_unique, labels=labelsy )

  }

  k<-ncol(x)
  # Compute all the column integrals
  intLogX <-NULL
  for (ki in (1:(ncol(logX1)))){
    print(paste('ind=', ki))
    yy=as.numeric(logX1[,ki]);
    integrals=integralFun(x_time, yy, subject_id, a=ini_time, b=end_time)
    intLogX<-cbind(intLogX, matrix(integrals));
  }
  print(dim(intLogX))


  m<-length(y_unique)
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
      lrX[,lloc] <- intLogX[,i]-intLogX[,j]
      lrcolnames<-c(lrcolnames,paste(paste("lr",i,sep=""),j,sep="."))
    }
  }
  colnames(lrX)<-lrcolnames



  y1<-y

  accuracy<-rep(0,niter)
  for(i in (1:niter)){
    y1<-sample(y1)
    lr<-coda_glmnet_longitudinal0(x,lrX,idlrX,nameslrX,y=y1, x_time, subject_id, ini_time, end_time,
                                alpha=alpha,lambda=lambda,covar=covar,nfolds = nfolds)
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
    "confidence interval"=quantile(accuracy, c(sig/1,1-(sig/2)))
  )
  return(results)
}


#' plot_signature_curves
#'
#' Graphical representation of the signature trajectories
#'
#' @param varNum column number of the variables (taxa)
#' @param coeff coefficients (coefficients must sum-up zero)
#' @param x microbiome abundance matrix in long format
#' @param y binary outcome; data type: numeric, character or factor vector
#' @param x_time observation times
#' @param subject_id subject id
#' @param ini_time initial time to be analyzed
#' @param end_time end time to be analyzed
#' @param color color graphical parameter
#' @param showLabel graphical parameter (see help(label))
#' @param location graphical parameter (see help(label))
#' @param inset graphical parameter (see help(label))
#' @param cex graphical parameter (see help(label))
#' @param y.intersp graphical parameter (see help(label))
#' @param main_title title plot
#'
#' @return trajectories plot
#'
#' @export
#'
#'
#'
#' @importFrom plyr d_ply
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom grDevices recordPlot
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' x=matrix(c(2, 3, 4, 1, 2, 5, 10, 20, 15, 30, 40, 12), ncol=2)
#' x_time = c(0,10,20,1,15, 25)
#' subject_id = c(1,1,1,2,2,2)
#' y=c(0,0,0,1,1,1)
#' plot_signature_curves(varNum=c(1,2), coeff=c(1,-1), x, y,x_time, subject_id,
#'                        ini_time=0, end_time=25)
#'
#-------------------------------------------------------------------------------
plot_signature_curves<-function(varNum, coeff, x,y, x_time, subject_id, ini_time, end_time,
                                color=c("chocolate1","slateblue2"), showLabel=TRUE, location="bottomright", inset=c(0.01,0.02), cex=0.8, y.intersp=0.7, main_title=NULL){
  # library(plyr)
  # suppressWarnings()
  X1 <- impute_zeros(x)

  #subject_id<-as.numeric(as.factor(subject_id))

  logX1 = log(X1);
  xref=sort(unique(x_time));

  y.binary<-ifelse(dim(table(y))==2, TRUE, FALSE)

  logcontrast=as.matrix(logX1[,varNum])%*%coeff
  # logcontrast<-logcontrast-mean(logcontrast)
  yy<-logcontrast

  positive<-varNum[coeff>0]
  negative<-varNum[coeff<0]

  a=ini_time;
  b=end_time;
  indxPlot <- which((x_time>=a) & (x_time <=b))
  ymin<-min(yy[indxPlot])
  ymax<-max(yy[indxPlot])

  k1<-positive
  k2<-negative

  k1string<-paste(k1, collapse = ',')
  k2string<-paste(k2, collapse = ',')

  if (is.null(main_title)){
    main_title<-paste("Signature",k1string,"vs",k2string)
  }

  if (y.binary==TRUE){
    y<-factor(y)
  plot(x_time,yy, xlim=c(a, b), ylim=c(ymin, ymax), main = main_title, ylab = "", xlab = "time");
    if (showLabel==TRUE){
      graphics::legend(location, inset=inset,c(levels(as.factor(y))[1],levels(as.factor(y))[2]),
             col=color,pch=c(1,1), cex = cex, bty="n")

    }
  } else {
  plot(x_time,yy, xlim=c(a, b), ylim=c(ymin, ymax), main = main_title, ylab = "", xlab = "time");
  }
  df<-data.frame(x_time,yy, subject_id,y)
  colnames(df)=c("time", "yy", "id","Y")

  # colores <- c("orchid","gold1")
  #colores <- c("chocolate1","slateblue2")
  colores <- color
  #c("F3AE05","B09EC5")
  # d_ply(df,"id",function(x) lines(x$time, x$yy, lty=3, col=as.numeric(x$Y)+2))

  plyr::d_ply(df,"id",function(x) graphics::lines(x$time, x$yy, lty=3, col=ifelse(as.numeric(x$Y)==1,colores[1],colores[2])))
  n = length(unique(df$id)); #num of curves

  df0<-df

  if (y.binary == TRUE){
  df0 <- subset(df, as.factor(y)==levels(as.factor(y))[1])
  }
  xref0 = sort(unique(df0$time));
  id0 <- unique(df0$id)  # ***
  yref0 = rep(0,length(xref0));
  valors<-NULL
  for (i in (1:length(xref0))){
    valors<-NULL
    val=xref0[i];
    # nn=0;
    for (k in id0){  # ***
      curve <- subset(df0, df0$id==k)
      if (length(curve$time)>0){
        if ((min(curve$time)<=val) && (max(curve$time)>=val)){
          if (length(curve$time)>1){
            pt <- approx(curve$time,curve$yy,val);
            valors=c(valors,pt$y);
          }
        }
      }
    }
    if (length(valors)>0){
      yref0[i]=mean(valors); ###############
    }
  }
  # lines(xref0, yref0,col=3, pch=19, type="o")
  graphics::lines(xref0, yref0,col=colores[1], pch=19, type="o")

  if (y.binary ==TRUE){
  df1<-subset(df, as.factor(y)==levels(as.factor(y))[2])
  id1<-unique(df1$id)  # ***
  xref1=sort(unique(df1$time));
  yref1=rep(0,length(xref1));
  valors<-NULL
  for (i in (1:length(xref1))){
    valors<-NULL
    val=xref1[i];
    nn=0;
    for (k in id1){   #***
      curve <- subset(df1, df1$id==k)
      if (length(curve$time)>0){
        if ((min(curve$time)<=val) && (max(curve$time)>=val)){
          if (length(curve$time)>1){
            pt <- approx(curve$time,curve$yy,val);
            valors=c(valors, pt$y);
          }
        }
      }
    }
    if (length(valors)>0){
      yref1[i]=mean(valors); ###############
    }
  }
  # lines(xref1, yref1, col=4, pch= 19, type="o", lwd=3)
  graphics::lines(xref1, yref1, col=colores[2], pch= 19, type="o", lwd=3)
  }
  graphics::abline(h=0)

  if (y.binary ==TRUE){
  if (showLabel==TRUE){
    graphics::legend(location, inset=inset,c(levels(as.factor(y))[1],levels(as.factor(y))[2]),
           col=color,pch=c(1,1), cex = cex, bty="o", y.intersp = y.intersp)

  }
  }

  my_plot=grDevices::recordPlot();
  return(my_plot);

}


#' plotMedianCurve
#'
#' Internal plot function
#'
#'
#' @param iNum .
#' @param iDen .
#' @param X .
#' @param Y .
#' @param x_time .
#' @param subject_id .
#' @param ini_time .
#' @param end_time .
#'
#' @return .
#' @importFrom plyr d_ply
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom grDevices recordPlot
#'
#'
#'
#' @author M. Calle - T. Susin
#'
#-------------------------------------------------------------------------------
plotMedianCurve<-function(iNum, iDen, X,Y, x_time, subject_id, ini_time, end_time){
  # library(plyr)
  X1 <- impute_zeros(X)

  logX1 = log(X1);
  xref=sort(unique(x_time));

  #subject_id<-as.numeric(as.factor(subject_id))


  a=ini_time;
  b=end_time;

  #
  k1=iNum  # First variable
  kplus=length(k1);
  k2=iDen  # Second variable
  kminus=length(k2);


  # idnot00 <- which(abs(rowSums(X[,k1])+rowSums(X[,k2]))>0)
  colNum = logX1[ ,k1];
  if (kplus==1){
    sumColNum = colNum;
  }else{
    sumColNum = rowSums(colNum);
  }
  numerator_sums = matrix(sumColNum/kplus);
  colDen = logX1[ ,iDen];
  if (kminus==1){
    sumColDen = colDen;
  }else{
    sumColDen = rowSums(colDen);
  }
  denominator_sums = matrix(sumColDen/kminus);
  yy=numerator_sums - denominator_sums;

  k1string<-paste(k1, collapse = ',')
  k2string<-paste(k2, collapse = ',')

  plot(x_time,yy, xlim=c(a, b), main = paste("Balance",k1string,"vs",k2string), ylab = "", xlab = "time");
  df<-data.frame(x_time,yy, subject_id,Y)
  colnames(df)=c("time", "yy", "id","Y")
  d_ply(df,"id",function(x) graphics::lines(x$time, x$yy, lty=3, col=as.numeric(x$Y)+2))
  n = length(unique(df$id)); #num of curves

  df0 <- subset(df, as.numeric(Y)==min(as.numeric(Y)))
  xref0 = sort(unique(df0$time));
  id0 <- unique(df0$id)  # ***
  yref0 = rep(0,length(xref0));
  valors<-NULL
  for (i in (1:length(xref0))){
    valors<-NULL
    val=xref0[i];
    # nn=0;
    for (k in id0){  # ***
      curve <- subset(df0, df0$id==k)
      if (length(curve$time)>0){
        if ((min(curve$time)<=val) && (max(curve$time)>=val)){
          if (length(curve$time)>1){
            pt <- approx(curve$time,curve$yy,val);
            valors=c(valors,pt$y);
          }
        }
      }
    }
    if (length(valors)>0){
      yref0[i]=median(valors);
    }
  }
  graphics::lines(xref0, yref0,col=3, pch=19, type="o")


  df1<-subset(df, as.numeric(Y)==max(as.numeric(Y)))
  id1<-unique(df1$id)  # ***
  xref1=sort(unique(df1$time));
  yref1=rep(0,length(xref1));
  valors<-NULL
  for (i in (1:length(xref1))){
    valors<-NULL
    val=xref1[i];
    nn=0;
    for (k in id1){   #***
      curve <- subset(df1, df1$id==k)
      if (length(curve$time)>0){
        if ((min(curve$time)<=val) && (max(curve$time)>=val)){
          if (length(curve$time)>1){
            pt <- approx(curve$time,curve$yy,val);
            valors=c(valors, pt$y);
          }
        }
      }
    }
    if (length(valors)>0){
      yref1[i]=median(valors);
    }
  }
  #lines(xref, yref0,col=3, pch=19, type="o", lwd=3)
  graphics::lines(xref1, yref1, col=4, pch= 19, type="o", lwd=3)
  graphics::abline(h=0)

  my_plot=recordPlot();
  return(my_plot);

}

#' filter_longitudinal
#'
#' Filters those individuals and taxa with enough longitudinal information
#'
#' @param x abundance matrix or data frame in long format (several rows per individual)
#' @param taxanames names of different taxa
#' @param x_time observation times
#' @param subject_id subject id
#' @param metadata matrix sample data
#' @param ini_time initial time to be analyzed
#' @param end_time end time to be analyzed
#' @param percent_indv percentage of individuals with more than min_obs observations
#' @param min_obs required minimum number of observations per individual
#'
#' @return list with filtered abundance table, taxanames and metadata
#' @export
#'
#'
#' @author M. Calle - T. Susin
#'
#' @examples
#'
#' data(ecam_filtered, package = "coda4microbiome")   # load the data
#'
#' x=x_ecam # microbiome abundance
#' x_time = metadata$day_of_life    # observation times
#' subject_id = metadata$studyid   # subject id
#' ini_time = 0
#' end_time = 360
#'
#' data_filtered<-filter_longitudinal(x,taxanames,x_time, subject_id, metadata,
#'                                               ini_time, end_time, min_obs=4)
#'
#'
#-------------------------------------------------------------------------------
filter_longitudinal<-function(x,taxanames=NULL,x_time, subject_id, metadata,ini_time=min(x_time), end_time=max(x_time), percent_indv=0.7, min_obs=3){

  id1<-names(which(table(subject_id)>=min_obs))
  index1<-which(subject_id%in%id1)
  metadata<-metadata[index1,]
  x<-x[index1,]
  x_time<-x_time[index1]
  subject_id<-subject_id[index1]

  rownames(metadata)<-NULL
  indtime<-which((x_time>=ini_time)&(x_time<=end_time))
  metadata0<-metadata[indtime,]
  subject_id0<-subject_id[indtime]
  x0<-x[indtime,]

  presencex<-x0
  presencex[x0>0]<-1

  nindividuals<-length(table(subject_id0))
  numberobs<-matrix(0,nindividuals,ncol(x))
  for(i in (1:ncol(x))){
    numberobs[,i]<-tapply(presencex[,i],subject_id0,sum)
  }


  nobs<-apply(numberobs,2,function(x) sum(x>=min_obs))  # for each taxa, number of individuals with more than 3 observations
  indtaxa<-which(nobs>=percent_indv*nindividuals)

  indsubject<-NULL
  for (i in (1:length(subject_id0))){
    if (subject_id0[i] %in% subject_id0[-i]){
      indsubject<-c(indsubject,i)
    }
  }

  x<-x0[indsubject,indtaxa]
  if (!is.null(taxanames)){
  taxanames<-taxanames[indtaxa]
  } else {
    taxanames<-colnames(x)
  }
  metadata<-metadata0[indsubject,]

  #save(x,taxanames, metadata, file="data_filtered.RData")

  results <- list(
    "filtered abundance matrix" = x,

    "filtered taxa names" = taxanames,

    "filtered metadata" = metadata
  )

  return(results)


}



