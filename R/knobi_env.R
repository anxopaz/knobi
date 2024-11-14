#' @title KBPM environmental analysis
#'
#' @name knobi_env
#'
#' @description Analyze and model the relationships between surplus production (SP) and environmental covariable(s) to test whether productivity changes in response to environmental fluctuations in three steps: (1) correlation analysis between the environmental variable(s) at different delays (lags) and the KBPM residuals through Pearson's correlation or autoregressive models; (2) selection of which lagged environmental variable(s) is included in the environmental KBPM models fit; (3) KBPM environmental fitting, where environmental effects are included as additive and multiplicative effects in the KBPM formulation (see details).
#'
#' @param knobi_results The output object of \code{\link{knobi_fit}} function (main package function).
#' @param data A list containing the following data: \itemize{
#' \item env: data frame containing the values of every environmental variable(s) in each column. Rows represent years.
#' \item years: years in which the environmental variable(s) are reported.}
#' @param control Optional. List containing the following settings: \itemize{
#' \item nlag: this argument is used to test, in the correlation analysis, the lags smaller or equal to 'nlag' (natural number). This means that correlation between KBPM residuals_{t} and X_{t-lag}, being X the environmental variable and lag the selected value from sequence {0,1,...,nlag}, is computed. The lag corresponding to the highest correlation among the KBPM residuals and the corresponding time lagged environmental covariable is considered in the environmental model (unless otherwise specified in 'selected_var'). By default, 'nlag=3'. See details.
#' \item lag: optional numerical vector providing the lag value(s) to consider in the relation between the KBPM surplus production and the environmental variable(s). The length of this argument must be equal to the number of environmental variables included. Applies only if 'nlag' argument is not provided.
#' \item start_c: optional. Numerical vector providing the starting values of the environmental parameter 'c' for the additive and multiplicative models, respectively. By default, start_c=c(1,1). See details.
#' \item ar_cor: optional. Logical. By default this argument is FALSE, meaning that the correlation between the KBPM residuals and the environmental variable(s) is analyzed through Pearson correlation test, as described above. If this argument is "TRUE", the relationship  between the KBPM residuals and the environmental variables is analyzed by fitting autoregressive (AR) models including each one of the environmental time lagged variables as explanatory covariables. The environmental variable whose model reports the lowest Akaike information criterion (AIC) is selected to be included in the environmental KBPM fit. See details.
#' \item plot3d: optional. Logical. If this argument is TRUE, 3D plots reporting the surplus production curve conditioned to a grid of environmental values are provided. FALSE by default.
#' \item selected_var: optional. Character. By default, the fit is carried out  including the highest correlated environmental variable at the corresponding time lag. However, if this argument is equal to the name of one of the environmental variables, this variable is used in the environmental fit considering the lag derived from the correlation analysis.
#' \item multicovar: optional. Logical. TRUE  means that the environmental model includes all the input environmental covariables at the same time, up to a maximum of 5. By default this argument is FALSE, which means that only the environmental variable reporting the highest correlation is included (after lagging it if corresponds).}
#' @param plot_out Logical. TRUE means that a file with the plot of the environmental fits is created. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_dir Optional. Directory to create the folder and save the plots. Required when 'plot_out=TRUE'. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_filename Optional. Name of the folder that will contain the plots. Required when 'plot_out=TRUE'. The default value is the input in the \code{\link{knobi_fit}} function.
#'
#' @details
#' Additive environmental model adds the following term on the right hand of equation (1) or (2) described in \code{\link{knobi_fit}} function: \eqn{cX_{t-lag}B_{t}}, being \eqn{X_{t-lag}} the environmental variable and \eqn{B_{t}} the biomass or SSB at time \eqn{t-lag}.
#' Multiplicative environmental model multiplies the right hand of equation (1) or (2) by \eqn{exp(cX_{t-lag})}.
#'
#' If ar_cor argument is "TRUE", the correlation analysis between the KBPM residuals and the environmental variable(s) is carried out as follows.
#' First, an AR model is fitted to the residuals.
#' \deqn{r_t=\sum_{i=1}^{p}\beta_{i}r_{t-i}+\epsilon_{t}}
#' being \eqn{r_t} the KBPM base residual for year \eqn{t} and \eqn{p} the AR model order, estimated as the maximum time lag at which the absolute value of the residuals partial autocorrelation is large than \eqn{qnorm(0.975)/\sqrt(N_r)}, being \eqn{N_r} the length of the residuals series.
#' Then, AR models are fitted considering each one of the lagged environmental variable(s),
#' \deqn{r_{t,lag}=\sum_{i=1}^{p}\beta_{i}r_{t-i}+X_{t-lag}+\epsilon_{t}}, for \eqn{lag=0,1,...,nlag}.
#' being \eqn{X_{t,lag}} the lagged environmental variable at year {t-lag}. Then, we have an autoregressive model for each of the lagged environmental variables.
#' The AIC values of the above models are compared, and the lagged environmental variable whose model reports the lowest AIC is used in the KBPM fit, except if the argument 'lag' is used.
#'
#' @return A list containing the results of the three-step environmental analysis is provided. \itemize{
#' \item add: estimates of the additive model parameters.
#' \item mult: estimates of the multiplicative model parameters.
#' \item BRPs: reference points (RPs) estimates for each model assuming X_t=0 (see vignettes for RPs equations depending on X_t).
#' \item df: data frame with the information used in the fit.
#' \item selected_var: environmental variable(s) used in the fit.
#' \item selected_lag: data frame providing the time lag of the environmental variable(s) in the KBPM fit.
#' \item lag_cor: correlation between the environmental variable(s) and the KBPM residuals for each one of the time lags.
#' \item env_aic: if 'ar_cor=TRUE', AIC values of each one of the autoregressive models (see details).
#' \item scaled_var: standardized environmental variable(s) used in the fit, i.e. the environmental variable(s) minus its mean and divided by its standard deviation (sd).
#' \item plots3D: list with the 3D plots objects (if 'plot3d=TRUE').
#' \item residuals: Pearson residuals from each model fit (base KBPM, additive model and multiplicative model).
#' \item error_table: array of performance and accuracy measures for each model: Standard error of the regression (SER), coefficient of determination (R-squared), adjusted coefficient of determination (adj-R-squared), Akaike information criterion (AIC), root-mean-squared error (RMSE) and mean absolute percentage error (MAPE) of observed and estimated surplus production values; and the statistic and p-value of the F-test. p-values lower than the considered significance level lead to the rejection of the null hypothesis, which implies that the KBPM environmental model is statistically better than the KBPM model that does not consider environmental information (The null hypothesis states that the environmental model fits the data as well as the base one).}
#' Results are also reported through plots which are displayed in the plot window and also saved (if 'plot_out=TRUE') in the provided directory or in the same directory as knobi_fit's plots.
#' The first plot reports the correlation analysis between the environmental variable(s) and the KBPM residuals. The second one reports the fitted SP values derived from the base model (no environmental information) and from the environmental ones.
#' If 'multicovar=FALSE' and 'plot3d=TRUE', 3D plots reporting the surplus production curve conditioned to a grid of environmental values are also reported.
#' If 'multicovar=TRUE', a plot with the Pearson correlation between the environmental variables is reported too.
#'
#' @author
#' \itemize{
#' \item{Anxo Paz}
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{M. Grazia Pennino}
#' }
#'
#' @examples
#'
#' \dontrun{
#'
#' # First, run the example of knobi_fit function
#'
#' # Then, provide environmental data series
#'
#' Env <- knobi_data$Env
#'
#' # The environmental data series must start in the first year of the KBPM fit data
#' # minus the provided nlag or lag
#' years <- knobi_results$df$Year # See knobi_fit example to obtain the knobi_results object
#' ind <- which(Env[,1]==years[1])
#' ind1 <- which(Env[,1]==years[length(years)])
#' nlag <- 5
#' Env <- Env[(ind-nlag):ind1,]
#'
#' # Now we create the environmental list
#' data <- list(env=data.frame(AMO=Env$AMO,Tmax=Env$TMax),
#'            years=Env$years)
#' control <- list(nlag=nlag)
#'
#' knobi_environmental <- knobi_env(knobi_results,data,control)
#' knobi_environmental
#' }
#'
#' @export


utils::globalVariables(c("correlation","AIC","SP"))

knobi_env <- function( knobi_results, data, control=NULL, plot_out=FALSE, plot_filename=NULL, plot_dir=NULL){

  if( is.null( control)) { control <- list()}

  if( is.null( control$ar_cor)) { control$ar_cor <- FALSE}

  input <- knobi_results$input
  bdata <- knobi_results$df
  pella <- knobi_results$control$pella

  env0 <- as.data.frame(data$env)
  env_names <- names(env0)
  y_env <- data$years
  res <- knobi_results$residuals

  if( is.null(control$start_c)){ start_c <- c(1,1)} else { start_c <- control$start_c}

  if(plot_out==T){

    old_dir <- getwd()

    if (is.null(plot_dir)) plot_dir <- knobi_results$control$plot_settings$plot_dir
    if (is.null(plot_filename)) plot_filename <- knobi_results$control$plot_settings$plot_filename

    setwd(paste0(plot_dir,"/",plot_filename))

  }

  df <- cbind( res, bdata)

  f_year <- df$Year[1]; l_year <- df$Year[length(df$Year)]

  env <- list()
  res_env <- list()

  df_env <- df

  res_env$selected_lag <- array( NA, dim=c(length(env_names),2))
  colnames(res_env$selected_lag) <- c("lag","correlation")
  rownames(res_env$selected_lag) <- env_names

  if ( is.null(control$lag)){ lag <- ifelse( is.null(control$nlag), 3, control$nlag)
  } else { lag <- max(control$lag); res_env$selected_lag[,1] <- control$lag}

  # Env data ----------------------------

  data_env <- list()

  res_env$lag_cor <- array( NA, dim=c(length(env_names), lag+1))

  vec_env <- "lag_0"

  for (i in 1:lag) vec_env <- c(vec_env,paste0("lag_",i))

  colnames(res_env$lag_cor) <- vec_env
  rownames(res_env$lag_cor) <- env_names


  for(j in env_names){

    data_env[[j]] <- df

    ind <- which( y_env == f_year)
    ind1 <- which( y_env == l_year)

    if( length(ind) > 0){ data_env[[j]]$env0 <- env0[ind:ind1,j]} else {
      warning('The length of the environmental variable is not enough to use the number of input lags')}

    for ( i in 1:lag){
      ind <- which( y_env == (f_year-i))
      ind1 <- which( y_env == (l_year-i))
      if(length(ind) > 0){ data_env[[j]][,5+i] <- env0[[j]][ind:ind1]} else {
          warning('The length of the environmental variable is not enough to use the number of input lags')}
    }

    bname <- ifelse( knobi_results$control$method=="Biomass", 'B', 'SSB')

    colnames(data_env[[j]]) <- c("res","SP",bname,"years",vec_env)

    if(!is.na(input$Recruitment[1])) data_env[[j]]$R <- input$Recruitment

    env[[j]] <- data_env[[j]][,c("res",vec_env)]

    cor <- round(cor(as.matrix(env[[j]]),use="na.or.complete"),4)
    cor <- cor[,1]; cor <- cor[-1]

    if (is.null(control$lag)){

      res_env$selected_lag[j,2] <- cor[[which(max(abs(cor))==abs(cor))[1]]]
      res_env$selected_lag[j,1] <- which(max(abs(cor))==abs(cor))[1]-1

    } else {

      res_env$selected_lag[j,2] <- cor[res_env$selected_lag[j,1]+1]

    }

    df_env[,j] <- scale(env[[j]][,res_env$selected_lag[j,1]+2])

    res_env$lag_cor[j,] <- cor

  }


  # AR cor ------------------

  if( control$ar_cor==TRUE){

    pacf_res <- stats::pacf( res, plot=FALSE)$acf[,1,1]

    ref <- stats::qnorm(0.975)/sqrt(length(res))

    auto <- max(which(abs(pacf_res)>=ref),0)
    if(auto==0){ warning("KBPM base residuals are not autocorrelated. An AR(0) model is fitted for SP residuals.")}

    fit_base  <-  stats::arima0( res, order = c(auto,0,0))

    env_aic <- array( NA, dim=c(length(env_names),lag+1))
    colnames(env_aic) <- vec_env; rownames(env_aic) <- env_names

    for(j in env_names) for(i in vec_env) env_aic[j,i] <- stats::arima0( res, order = c(auto,0,0), xreg = env[[j]][,i])$aic

    env_aic <- cbind(base=rep(fit_base$aic,length(env_names)),env_aic)
    res_env$env_aic <- env_aic

    env_aic_c <- env_aic-fit_base$aic+2
    min_aic <- env_aic_c[which.min(env_aic_c)]
    if(min_aic>=0) warning("AR models considering environmental variable(s) do not really improve AR model considering only the residuals")

    colnames(res_env$selected_lag)[2] <- "aic"

    if (is.null(control$lag)){

      for(j in env_names){
        res_env$selected_lag[j,1] <- which(env_aic == min(env_aic[j,-1]), arr.ind=TRUE)[2]-2
        res_env$selected_lag[j,2] <- min(env_aic[j,-1])
      }

    } else {

      for(j in env_names){
        res_env$selected_lag[j,2] <- env_aic[j,(res_env$selected_lag[j,1]+2)]
      }
    }

    for(j in env_names){
      df_env[,j] <- scale(env[[j]][,res_env$selected_lag[j,1]+2])
    }


  }


  # corr plot ----------------

  lagf <- NULL
  corlist <- NULL

  if( control$ar_cor == FALSE){

    for(i in vec_env){
      lagf <- c(lagf,rep(i,length(env_names)))
      corlist <- c(corlist,res_env$lag_cor[,i])
    }

    envcorplot_df <- data.frame(correlation=corlist,lag=lagf,factor=rep(env_names,length(vec_env)))

    envcorplot <- ggplot2::ggplot(data=envcorplot_df,ggplot2::aes(x=lag,y=correlation,group=factor,color=factor)) +
      ggplot2::theme_bw() + ggplot2::geom_point() + ggplot2::geom_line(linetype = "dashed") + ggplot2::ylim(-1,1) +
      ggplot2::labs(title="Environmental correlation with base KBPM SP residuals", subtitle=knobi_results$input$Stock,
                    y="Correlation",x="") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

    print(envcorplot)

  } else {

    for(i in vec_env){
      lagf <- c(lagf,rep(i,length(env_names)))
      corlist <- c(corlist,res_env$env_aic[,i])
    }

    maximo <- max(env_aic)+2
    minimo <- min(env_aic)-2

    envcorplot_df <- data.frame(AIC=corlist,lag=lagf,factor=rep(env_names,length(vec_env)))

    envcorplot <- ggplot2::ggplot(data=envcorplot_df,ggplot2::aes(x=lag,y=AIC,group=factor,color=factor)) +
      ggplot2::theme_bw() + ggplot2::geom_point() + ggplot2::geom_line(linetype = "dashed") + ggplot2::ylim(minimo,maximo) +
      ggplot2::labs(title="AIC comparison", subtitle=knobi_results$input$Stock,
                    y="AIC",x="") +
      ggplot2::geom_hline(yintercept = fit_base$aic) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                     legend.background = ggplot2::element_rect(fill = "transparent"))

    print(envcorplot)

  }


  if (plot_out==TRUE){
    p  <-  grDevices::recordPlot()
    grDevices::jpeg("corplot.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  multicovar <- ifelse( is.null(control$multicovar), FALSE, control$multicovar)

  # Fit  ---------------

  if(multicovar == FALSE){

    if( is.null( control$selected_var)){

      if( control$ar_cor == FALSE){ selected_var <- env_names[which.max(abs(res_env$selected_lag[,2]))]
      } else { selected_var <- env_names[which.min(abs(res_env$selected_lag[,2]))]}

    } else { selected_var <- control$selected_var}

  } else {

    selected_var <- env_names

    if(!is.null(control$plot3d)) control$plot3d <- NULL

  }


  envsel <- df_env[,selected_var]
  res_env$selected_var <- colnames( envsel) <- selected_var

  if( multicovar){

    cor <- round(cor(as.matrix(envsel),use="na.or.complete"),4)
    for(i in 1:ncol(cor)) cor[i,i] <- 0
    cor <- max(abs(cor))
    if(cor >= 0.35){
      p.mat  <-  corrplot::cor.mtest(envsel)$p
      corrplot::corrplot(round(cor(as.matrix(envsel),use="na.or.complete"),2),method='circle',
                         title="Covariables correlation",mar=c(0,0,2,0),addCoef.col = "black",
                         type="lower",diag=T,p.mat = p.mat, sig.level = 0.05)
      print(round(cor(as.matrix(envsel),use="na.or.complete"),4))
      cat("\n")
      warning("The covariables are highly correlated (see corrplot). To avoid estimation issues, consider removing variables or setting multicovar=FALSE.")
    }

  }

  Data <- list( data=df, start_r=as.numeric(knobi_results$params[['r']]),
    start_K=as.numeric(knobi_results$params[['K']]), start_c=start_c[1], start_p=1)

  kbpm_mult  <- kbpm_fit( Data, envsel, pella, 'Mult')
  res_env$mult$params  <-  kbpm_mult
  res_env$mult$BRPs  <-  BRP( res_env$mult, pella)

  kbpm_add  <- kbpm_fit( Data, envsel, pella, 'Add')
  res_env$add$params  <-  kbpm_add
  res_env$add$BRPs  <-  BRP( res_env$add, pella)

  envdf <- df
  envdf <- cbind( envdf, envsel)

  bv  <- envdf$base
  bv1  <-  envdf$mult <- predict_model( res_env$mult, df$x, pella, 'Mult', envsel)
  bv2  <-  envdf$add <- predict_model( res_env$add, df$x, pella, 'Add', envsel)

  res_env$scaled_var <- envsel

  if( multicovar == F) rownames(res_env$scaled_var) <- df_env$Year-as.numeric(res_env$selected_lag[selected_var,1])

  if(is.null(control$plot3d)){plots3d <- FALSE} else {plots3d <- control$plot3d}


  # 3D plots ---------------------------

  if(plots3d==TRUE){

    x  <-  envdf$x
    y  <-  envdf[[selected_var]]
    ysc <- envsel
    z  <-  envdf$y

    r_a <- kbpm_add['r']; K_a <- kbpm_add['K']; c_a <- kbpm_add['c']
    p_a <- ifelse( pella, kbpm_add['p'], 1)

    r_m <- kbpm_mult['r']; K_m <- kbpm_mult['K']; c_m <- kbpm_mult['c']
    p_m <- ifelse( pella, kbpm_mult['p'], 1)


    grid.lines  <-  400
    cut_a <- max( K_a+c_a*K_a*max(y), K_a+c_a*K_a*min(y))
    x.pred_a  <-  seq(0, cut_a, length.out = grid.lines)
    x.pred_m  <-  seq(0, K_m, length.out = grid.lines)

    y.pred  <-  c( seq(min(y), 0, length.out = grid.lines/2), seq( 0, max(y), length.out = grid.lines/2))
    xy_a  <-  expand.grid(x = x.pred_a, y = y.pred)
    xy_m  <-  expand.grid(x = x.pred_m, y = y.pred)

    z.pred_a  <-  (r_a/p_a) * xy_a$x * (1-(xy_a$x/K_a)^(p_a)) + c_a * xy_a$y * xy_a$x
    z.pred_m  <-  exp(c_m*xy_m$y) * ((r_m/p_m) * xy_m$x * (1-(xy_m$x/K_m)^p_m))
    z.pred_a  <-  matrix(z.pred_a,nrow=grid.lines,ncol=grid.lines)
    z.pred_m  <-  matrix(z.pred_m,nrow=grid.lines,ncol=grid.lines)

    y2  <-  y * attr(ysc, 'scaled:scale') + attr(ysc, 'scaled:center')
    y.pred2  <-  y.pred * attr(ysc, 'scaled:scale') + attr(ysc, 'scaled:center')

    cat("\n This can take a while... \n")

    plot3D::scatter3D(x, y2, z, bty="b2", pch = 19, cex = 1, cex.axis = 0.6,
                      colkey = list(length = 0.5, width = 0.5, cex.clab = 0.8, cex.axis = 0.8),
                      xlim = c(0,max(x.pred_a)), zlim = c(0,max(z.pred_a)*1.5),
                      theta = 28, phi = 20, ticktype = "detailed", clim = c(0,max(z.pred_a,z)),
                      xlab = "SSB", ylab = selected_var, zlab = "SP", clab = "Surplus production",
                      surf = list(x = x.pred_a, y = y.pred2, z = z.pred_a, facets = NA, fit = z),
                      main = "Additive model: Production curve", sub = knobi_results$input$Stock)

    plot3d_add <- grDevices::recordPlot()

    if (plot_out == TRUE){
      grDevices::jpeg("additive_model.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(plot3d_add)
      grDevices::dev.off()
    }

    cat("\n ... Just a little more. May the 4th be with you ... \n")

    plot3D::scatter3D(x, y2, z, bty="b2", pch = 19, cex = 1, cex.axis = 0.6,
                      colkey = list(length = 0.5, width = 0.5, cex.clab = 0.8, cex.axis = 0.8),
                      xlim = c(0, max(x.pred_m)), zlim = c(0,max(z.pred_m)*1.5),
                      theta = 28, phi = 18, ticktype = "detailed", clim=c(0,max(z.pred_m,z)),
                      xlab = "SSB", ylab = selected_var, zlab = "SP", clab = "Surplus production",
                      surf = list(x = x.pred_m, y = y.pred2, z = z.pred_m, facets = NA, fit = z),
                      main = "Multiplicative model: Production curve", sub = knobi_results$input$Stock)
    plot3d_mult <- grDevices::recordPlot()

    if (plot_out==TRUE){
      grDevices::jpeg("multiplicative_model.jpeg",width=2500, height=2500,res=300)
      grDevices::replayPlot(plot3d_mult)
      grDevices::dev.off()
    }

    cat(paste0("\n ... Done! :) \n"))


    res_env$plots3D <- list(additive_plot=plot3d_add,multiplicative_plot=plot3d_mult)

  }

  envplot_df <- data.frame(SP=c(df$y,bv,bv2,bv1),Year=rep(df$Year,4),
                         factor=c(rep("Observed",length(bv)),rep("Base model",length(bv)),
                                  rep("Environmental Additive",length(bv)),
                                  rep("Environmental Multiplicative",length(bv))))

  env_plot <- ggplot2::ggplot(data=envplot_df,ggplot2::aes(x=Year,y=SP,color=factor)) + ggplot2::theme_bw() +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min(envplot_df$SP),max(envplot_df$SP)) +
    ggplot2::labs(title="Environmental fits", subtitle=knobi_results$input$Stock,
                  y="Surplus Production") +
    ggplot2::theme(legend.position = c(0.15,0.85), plot.title = ggplot2::element_text(hjust = 0.5),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),legend.title=ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent"))

  print(env_plot)

  if (plot_out==TRUE){
    p  <-  grDevices::recordPlot()
    grDevices::jpeg("fits_env.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  errors  <-  kbpm_error( knobi_results, res_env, plot_out)

  res_env$error_table <- errors$error_table
  res_env$residuals <- errors$residuals

  if(multicovar){
    mscale <- mcenter <- 1:ncol(res_env$scaled_var)
    for (i in 1:ncol(res_env$scaled_var)){
      mscale[i] <- attr( df_env[,selected_var[i]], 'scaled:scale')
      mcenter[i] <- attr( df_env[,selected_var[i]], 'scaled:center')}
    attr( res_env$scaled_var, 'scaled:scale') <- mscale
    attr( res_env$scaled_var, 'scaled:center') <- mcenter}

  res_env$input <- list( control = control, data = data, basecontrol = knobi_results$control)
  res_env$base <- list( params = knobi_results$params, BRPs = knobi_results$BRPs)

  res_env$BRPs <- rbind( res_env$base$BRPs, res_env$add$BRPs, res_env$mult$BRPs)
  rownames( res_env$BRPs) <- c( 'Base', 'Add', 'Mult')

  envdf <- envdf

  res_env$df <- envdf

  if (plot_out==TRUE){
    cat(paste0("\n Plots successfully saved in '",getwd(),"'"),". \n")
    setwd(old_dir)
  }

  class( res_env) <- 'knobi'

  return( res_env)

}
