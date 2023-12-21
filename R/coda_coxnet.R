#'coda_coxnet
#'
#' Microbial signatures in survival studies
#' The algorithm performs variable selection through an elastic-net penalized Cox regression conveniently adapted to CoDA.
#' The result is expressed as the (weighted) balance between two groups of taxa.
#' It allows the use of non-compositional covariates.
#'
#'
#' @param x abundance matrix or data frame (rows are samples, columns are variables (taxa))
#' @param time time to event or follow up time for right censored data. Must be a numericvector.
#' @param status event occurrence. Vector (type: numeric or logical) specifying 0, or FALSE, for no event occurrence, and 1, or TRUE, for event occurrence.
#' @param covar data frame with covariates (default = NULL)
#' @param lambda penalization parameter (default = "lambda.1se")
#' @param nvar number of variables to use in the glmnet.fit function (default = NULL)
#' @param alpha elastic net parameter (default = 0.9)
#' @param nfolds number of folds
#' @param showPlots if TRUE, shows the plots (default = TRUE)
#' @param coef_threshold coefficient threshold, minimum absolute value of the coefficient for a variable to be included in the model (default =0)
#'
#' @return list with "taxa.num","taxa.name","log-contrast coefficients","risk.score","apparent Cindex","mean cv-Cindex","sd cv-Cindex","risk score plot","signature plot".
#' @export
#' @import survival
#' @import glmnet
#' @import ggpubr
#'
#' @author M. Calle, M. Pujolassos, T. Susin
#'
#' @examples
#'
#' data(data_survival, package = "coda4microbiome")
#' time <- Event_time
#' status <- Event
#' set.seed(12345)
#' coda_coxnet(x = x,
#'            time = Event_time,
#'            status = Event,
#'            covar = NULL,
#'            lambda = "lambda.1se", nvar = NULL,
#'            alpha = 0.9, nfolds = 10, showPlots = TRUE, coef_threshold = 0)
#'
#'
#'
#-------------------------------------------------------------------------------
coda_coxnet <- function (x, time, status, covar = NULL, lambda = "lambda.1se", nvar = NULL,
                         alpha = 0.9, nfolds = 10, showPlots = TRUE, coef_threshold = 0)
{

#  library(survival)
#  library(glmnet)
#  library(ggpubr)

  # Abundance table
  x <- impute_zeros(x)
  kselect <- ncol(x)
  taxaselect <- (1:ncol(x))
  lrmatrix <- logratios_matrix(x)
  lrX <- lrmatrix[[1]]
  idlrX <- lrmatrix[[2]]
  nameslrX <- lrmatrix[[3]]

  # Response variable

  y = Surv(time, status)       # Cox response variable (time & status)


  if (is.null(covar)) { # cox glmnet WITHOUT covars

    cvfit <- glmnet::cv.glmnet(lrX, y, family = "cox", type.measure = "C",
                               alpha = alpha, nfolds = nfolds,
                               keep = TRUE)
  } else {                # cox glmnet for BINARY Y, WITH covars

    df0 <- data.frame(as.matrix(y), covar)
    model0 <- coxph(Surv(time, status) ~ ., data = df0) # coxph to integrate all covars as one into the cox model
    x0 <- predict(model0)

    cvfit <- glmnet::cv.glmnet(lrX, y, family = "cox",
                               type.measure = "C",
                               nfolds = nfolds, alpha = alpha,
                               keep = TRUE,
                               offset = x0)
  }
  if (showPlots == TRUE) {
    plot(cvfit)
  }
  if (!is.null(nvar)) {
    rowlasso <- max(which(cvfit$glmnet.fit$df <= nvar))
    lambda <- cvfit$glmnet.fit$lambda[rowlasso]
  }
  lambdavalue <- lambda
  if (is.character(lambda)) {
    if (lambda == "lambda.1se")
      lambdavalue <- cvfit$lambda.1se
    if (lambda == "lambda.min")
      lambdavalue <- cvfit$lambda.min
  }
  idrow <- max(which(cvfit$glmnet.fit$lambda >= lambdavalue))
  coeflr <- as.vector(coef(cvfit, s = lambda)) #[-1]
  lrselect <- which(coeflr != 0)
  coeflogcontrast <- rep(0, ncol(x))
  for (i in (1:length(coeflr))) {
    coeflogcontrast[idlrX[i, 1]] <- coeflogcontrast[idlrX[i,
                                                          1]] + coeflr[i]
    coeflogcontrast[idlrX[i, 2]] <- coeflogcontrast[idlrX[i,
                                                          2]] - coeflr[i]
  }
  varlogcontrast <- which(abs(coeflogcontrast) > coef_threshold)
  coeflogcontrast <- coeflogcontrast[varlogcontrast]

  (names.select <- colnames(x)[varlogcontrast])
  (sign <- ifelse(coeflogcontrast > 0, 1, 0))
  sign <- factor(sign, levels = c(0, 1), labels = c("negative",
                                                    "positive"))
  logcontrast = as.matrix(log(x)[, varlogcontrast]) %*% coeflogcontrast # Bacterial signature

  if (is.null(covar)) {
    predictions <- logcontrast
  } else {
    predictions<-x0+logcontrast
  }

  coeflogcontrast<-2*coeflogcontrast/sum(abs(coeflogcontrast))

  if (length(varlogcontrast) == 0){
    Cindex_signature <- 0.5
  } else {
    Cindex_signature <- glmnet::Cindex(pred=predictions, y)     # Apparent C-Index
  }
  mcvCindex <- cvfit$cvm[idrow]
  sdcvCindex <- cvfit$cvsd[idrow]

  #plot1 <- plot_prediction(predictions, y, showPlots = showPlots)
  plot1 <- plot_riskscore(predictions, x, time, status, showPlots = showPlots)
  plot2 <- plot_signature(names.select, coeflogcontrast, showPlots = showPlots)

  #if (y.binary == TRUE) {
  results <- list(taxa.num = varlogcontrast, taxa.name = names.select,
                  `log-contrast coefficients` = coeflogcontrast, risk.score = predictions,
                  `apparent Cindex` = Cindex_signature, `mean cv-Cindex` = mcvCindex,
                  `sd cv-Cindex` = sdcvCindex,
                  `risk score plot` = plot1,
                  `signature plot` = plot2)
  return(results)
}

###############################################################################

################################################################################
###################### Plot1: Signature Survival curves ########################
################################################################################
#' plot_survcurves
#'
#' Plots survival curves stratifying samples based on their signature value.
#' Signature value for stratification can be set by the user.
#'
#' @param risk.score microbial risk score values resulting from the coda_coxnet model
#' @param time time to event or follow up time for right censored data. Must be a vector (type:numeric) specifying time to event for each sample for right censored data (what about interval data?).
#' @param status event occurrence. Vector (type: numeric or logical) specifying 0 or FALSE for no event occurrence, and 1 or TRUE for event occurrence.
#' @param strata.quantile cut-off quantile (values 0, 1) for sample stratification based on signature value. Default is set to 0.5 quantile of the signature.
#'
#' @return return an object of class ggsurvplot which is list containing the following:
#' plot: the survival plot (ggplot object).
#' table: the number of subjects at risk table per time (ggplot object).
#' data.survplot: data used to plot the survival curves (data.frame).
#' data.survtable: data used to plot the tables under the main survival curves (data.frame).
#'
#' @export
#' @import survival
#' @import survminer
#'
#' @author M. Calle, M. Pujolassos, T. Susin
#'
#' @examples
#'
#' set.seed(12345)
#'
#' data(data_survival, package = "coda4microbiome")
#' time <- Event_time
#' status <- Event
#' crohn_cox <- coda_coxnet(x = x,
#'                          time = Event_time,
#'                          status = Event,
#'                          covar = NULL,
#'                          lambda = "lambda.1se", nvar = NULL,
#'                          alpha = 0.9, nfolds = 10, showPlots = TRUE, coef_threshold = 0)
#' plot_survcurves(risk.score = crohn_cox$risk.score,
#'                  time,
#'                  status,
#'                  strata.quantile = 0.5)
#'
#'
#' #-------------------------------------------------------------------------------

plot_survcurves <- function(risk.score, time, status, strata.quantile = 0.5)
{
#  library(survminer)
#  library(survival)
  c_sign = c("#F8766DFF", "#00BFC4FF")

  # Take predictions (risk.score) from the coda_coxnet model
  df <- data.frame(risk.score, status, time)

  # Set survival classification based on signature predictions (risk.score)
  if (is.null(strata.quantile)){
    strata.quantile = 0.5
    stratalabs = c(paste("below", strata.quantile, "quantile", sep = " "),
                 paste("above", strata.quantile, "quantile", sep = " "))

    cutoff = quantile(df$risk.score, strata.quantile)[[1]]

    df$class <- ifelse (df$risk.score <= cutoff, stratalabs[1], stratalabs[2])
    #df$class <- factor(df$class, levels = c(0, 1), labels = stratalabs)
  } else {
    if (0 < strata.quantile | strata.quantile < 1) {
      stratalabs = c(paste("below", strata.quantile, "quantile", sep = " "),
                   paste("above", strata.quantile, "quantile", sep = " "))

      cutoff = quantile(df$risk.score, strata.quantile)[[1]]

      df$class <- ifelse (df$risk.score <= cutoff, stratalabs[1], stratalabs[2])
      #df$class <- factor(df$class, levels = c(0, 1), labels = stratalabs)
    } else {
      print("strata.quantile must be between 0 and 1" )
    }
  }
  fit1 = survfit(Surv(time, status) ~ class, data=df)
  plotlabs = gsub("class=", "", names(fit1[["strata"]]))
  survplot <- ggsurvplot(fit1, data=df,
                         conf.int=TRUE, pval=TRUE, risk.table=TRUE,
                         risk.table.y.text = FALSE,
                         legend.labs = plotlabs,
                         palette = c_sign)
  return(survplot)
}

###############################################################################

################################################################################
###################### Plot2: Heatmap of Microbial risk scores #########################
################################################################################
#' plot_riskscore
#'
#' Plots samples ordered by microbial risk score values along time to event.
#'
#' @param risk.score microbial risk score values resulting from the coda_coxnet model
#' @param x original survival data
#' @param time time to event or follow up time for right censored data. Must be a vector (type:numeric) specifying time to event for each sample for right censored data.
#' @param status event occurrence. Vector (numeric or logical) specifying 0 (or FALSE) for no event occurrence, and 1 (or TRUE) for event occurrence.
#' @param showPlots (default: TRUE)
#'
#' @return returns an object of class HeatmapList.
#'
#' @export
#' @import ComplexHeatmap
#' @import circlize
#'
#' @author M. Calle, M. Pujolassos, T. Susin
#'
#' @examples
#'
#' set.seed(12345)
#'
#' data(data_survival, package = "coda4microbiome")
#' time <- Event_time
#' status <- Event
#' crohn_cox <- coda_coxnet(x = x,
#'                          time = Event_time,
#'                          status = Event,
#'                          covar = NULL,
#'                          lambda = "lambda.1se", nvar = NULL,
#'                          alpha = 0.9, nfolds = 10, showPlots = TRUE, coef_threshold = 0)
#' plot_riskscore(risk.score = crohn_cox$risk.score,
#'                     x = x,
#'                     time = Event_time,
#'                     status =  Event,
#'                     showPlots = TRUE)
#'
#'
#' #-------------------------------------------------------------------------------
plot_riskscore <- function(risk.score, x, time, status, showPlots = TRUE)
{
#  library(ComplexHeatmap)
#  library(circlize)

  # Set colors:
  c_sign = c("#F8766DFF", "#00BFC4FF")
  #c_sign = c("deepskyblue3", "coral1")
  c_status = c("gray95", "chocolate1")
  c_time = c("gray80", "chocolate1")

  # Arrange features data
    # Signature
  signature <- matrix(risk.score, dimnames = list(rownames(x), "signature"))
    # Status
  if (is.numeric(status)){
    e <- as.factor(status)
    event <- matrix(e, dimnames = list(rownames(x), "event"))
  }
  if (is.logical(status)){
    e <- factor(status, levels = c("FALSE", "TRUE"), labels = c("0", "1"))
    event <- matrix(e, dimnames = list(rownames(x), "event"))
  }
    # Time
  eventime <- matrix(time, dimnames = list(rownames(x), "time"))

  # Get the right order according to signature value
  samorder <- rownames(x)[order(signature[,1], decreasing = F)]
  signature <-signature[samorder,]
  event <- event[samorder,]
  eventime <- eventime[samorder,]

  # Heatmaps
  tsignature <- as.matrix(t(signature))
  rownames(tsignature) <- "Microbial risk score"
  hm_sign <- Heatmap(tsignature,
                     name = "Microbial risk score",
                     # cluster and dendogram
                     cluster_rows = F,
                     cluster_columns = F,
                     # Rows arrangements (samples)
                     #row_order = samorder,
                     row_names_side = "left",
                     show_row_dend = FALSE,
                     show_row_names = TRUE,
                     # Columns arrangements (signature)
                     column_title_side = "bottom",
                     show_column_names = FALSE,
                     show_column_dend = FALSE,

                     col =  colorRamp2(breaks = c(quantile(signature, 0)[[1]],
                                                  quantile(signature, 1)[[1]]),
                                       colors = c_sign),
                     border = FALSE
                     #width = unit(1, "cm")
  )

  tevent <- as.matrix(t(event))
  rownames(tevent) <- "Status"
  hm_event <- Heatmap(tevent,
                      name = "Status",
                      # cluster and dendogram
                      cluster_rows = F,
                      cluster_columns = F,
                      # Rows arrangements (samples)
                      ##row_order = samorder,
                      row_names_side = "left",
                      show_row_dend = FALSE,
                      show_row_names = TRUE,
                      # Columns arrangements (signature)
                      column_title_side = "bottom",
                      show_column_names = FALSE,
                      show_column_dend = FALSE,

                      col =  c("0" = c_status[[1]], "1" = c_status[[2]]),
                      border = FALSE
                      #width = unit(1, "cm")
  )
  point_time = HeatmapAnnotation(name = "Time",
                                 Time = anno_points(eventime,  #which = "row",
                                                    gp = grid::gpar(col = ifelse(event == "1", c_time[[2]], c_time[[1]])),
                                                    axis_param = list(direction = "reverse")),
                                 # Annotation labels
                                 height = unit(5, "cm"),
                                 annotation_name_side = "left",
                                 annotation_name_rot = 0,
                                 gap = unit(2, "mm"),
                                 border = TRUE
  )
  ht <- point_time %v% hm_event %v% hm_sign

  if (showPlots == TRUE) {
    draw(ht)
  }
  return(ht)
}

