#' @title Known Biomass Production Model (KBPM) fit
#'
#' @name knobi_fit
#'
#' @description This function, that is the main function of the \code{knobi} package, fits a type of surplus production models named known-biomass production models (KBPM) (MacCall, 2002). The surplus production curve is fitted using the observed catch time series and the biomass or SSB (Spawning Stock Biomass) estimates derived from the fit of a data-rich stock assessment model.
#'
#' @param data A list containing the data. \itemize{
#' \item Catch: catch time series observations.
#' \item Biomass: biomass time series estimated through a data-rich stock assessment model. If available, otherwise enter SSB in the next argument.
#' \item SSB: spawning stock biomass time series estimated through a data-rich stock assessment model. If available, otherwise introduce biomass in the previous argument.
#' \item years: optional. Years in which catch observations have been measured. By default, an increasing natural sequence from 1 to the total number of catch observations.
#' \item Stock: optional. Character string with the stock name for using in the plot subtitles.
#' \item Recruitment: optional. Recruitment time series estimated through a data-rich stock assessment model. See details.
#' \item F_input: optional. Fishing mortality time series estimated through a data-rich stock assessment model. See details.
#' \item RP: optional. Estimates of the following biological reference points derived from the fit of a data-rich stock assessment model (see details): \itemize{
#' \item MSY: Maximum Sustainable Yield.
#' \item F_MSY: fishing mortality at MSY.
#' \item B_MSY: biomass at MSY (or SSB at MSY depending on the value of the 'method' argument, see control list).
#' \item K: virgin biomass.}}
#' @param control A list containing the following control parameters. \itemize{
#' \item pella: Logical. TRUE means that Pella-Tomlinson model is used whereas FALSE means that Schaefer model is fitted (by default). See details.
#' \item start_r: optional. Start value of the growth rate model parameter 'r' (i.e. intrinsic rate of natural increase). See details.
#' \item start_K: optional. Start value of the carrying capacity model parameter 'K' (the maximum population size for growth to be positive). See details.
#' \item start_p: optional. Start value of the shape model parameter 'p'. Used for Pella-Tomlinson model. See details.
#' \item method: sets whether the fit is carried using 'SSB' or 'Biomass'. This argument is only required if both time series, 'SSB' and 'Biomass', are provided.}
#' @param plot_out Logical. TRUE means that files are created with the input time series plots and the fitting results. FALSE by default.
#' @param plot_dir Optional. Directory for saving the plot files. Required when plot_out=TRUE. Current directory by default.
#' @param plot_filename Optional. Name of the folder that will contain the plot files. By default, "knobi_results". Required when plot_out=TRUE.
#'
#'
#' @return The list of the results of the KBPM fit contains: \itemize{
#' \item params: model parameters estimates.
#' \item df: the data used for the model.
#' \item BRPs: biological reference points estimates: \itemize{
#' \item K: virgin biomass.
#' \item MSY: maximum sustainable yield (MSY).
#' \item B_MSY: biomass at MSY.
#' \item F_MSY: fishing mortality at MSY.
#' \item MSYoverK: ratio of MSY and K.}
#' \item residuals: Pearson's residuals calculated as (observations-estimates)/sd(observations).
#' \item error_table: data frame of accuracy and model performance measures: \itemize{
#' \item SER: standard error of the regression, calculated as the root of the rate between the sum of residual squares and the degrees of freedom of the regression.
#' \item R-squared: coefficient of determination.
#' \item adj-R-squared: adjusted coefficient of determination.
#' \item AIC: Akaike information criterion.
#' \item RMSE: root mean squared error (observed vs. estimated values).
#' \item MAPE: mean absolute percentage error (observed vs. estimated values).}
#' \item input: the input list updated with the annual average biomass ($input$Average_Biomass), the surplus production ($input$SP), and the F estimates derived from the fit, in ($input$F_output).
#' \item control: the updated control settings of the fit.
#' \item optimx: list of some results provided by \code{\link[optimx]{optimx}}: \itemize{
#' \item value: value of the objective function in the minimization process.
#' \item convergence: integer code. ‘0’ indicates successful completion in the optimization.
#' \item message: character string giving any additional information returned by the optimizer, or NULL.}
#' }
#' The plots are displayed in the plot window and also saved (if 'plot_out=TRUE') in the provided directory or in the current one (if 'plot_dir' is not provided). The following input quantities are plotted: fishing mortality time series, SSB or biomass, surplus production and catch time series. Plots of catch over fishing mortality, fishing mortality over SSB (or biomass), and catch over SSB (or biomass) time series with a smooth line from a "loess" regression are also available. Plot of input-output time series of fishing mortality includes horizontal lines for both, input and output, fishing mortalities at MSY (one line if input F_MSY is NULL). Fishing mortality relative to F_MSY is also plotted including an horizontal line at one (note that values greater than 1 indicate an F greater than F_MSY). The analogous SSB (or biomass) plots are also reported. On the other hand, the fitted surplus production curve is plotted twice with the SSB (or biomass) and SP observations (first) and with the catch and SP observations (second). Finally, a plot with the KBPM residuals is reported.
#'
#'
#' @author
#' \itemize{
#' \item{Anxo Paz}
#' \item{Marta Cousido-Rocha}
#' \item{Santiago Cerviño López}
#' \item{Maria Grazia Pennino}
#' }
#'
#' @details The KBPMs implemented in the current package are explained below.
#' Schaefer model (1954) (1):
#' \deqn{SP_{t} = r B_{t} (1-(B_{t}/K))}
#' where \eqn{SP_{t}} is the surplus production, \eqn{B_{t}} is the biomass or SSB averaged (mean of two consecutive years), \eqn{r} is the population growth rate parameter, and \eqn{K} is the virgin biomass. The subscript \eqn{t} denotes the time (years).
#' Pella and Tomlinson model (1969) (2):
#' \deqn{SP_{t} = (r/p) B_{t} (1-(B_{t}/K)^{p})}
#' where \eqn{SP_{t}} is the surplus production, \eqn{B_{t}} is the biomass or SSB averaged,  \eqn{r} is the population growth rate parameter, \eqn{K} is the virgin biomass and \eqn{p} is the asymmetry parameter. The subscript \eqn{t} denotes the time (years).
#'
#' The recruitment values, F_input estimates and F_type argument are included to get an overview of the stock status, but are not needed for the KBPM fitting. Similarly, the input reference points are only used for comparison purposes.
#'
#' KBPMs have also proven their usefulness for the multispecies management objectives. More precisely, KBPM approach can be applied to analyze the dynamic of the total aggregated biomass and catch of all targeted fish species in a community, defining in this way a simple data-limited ecosystem model to assess the ecosystem status (see example).
#'
#' @examples
#'
#' library(knobi)
#'
#' data( knobi_data)
#' hake_n <- knobi_data$hake_n
#'
#'
#' # Then, create the data object.
#'
#' data<-list()
#' data$SSB<-hake_n$SSB # We take the SSB in our data.
#' data$Catch<-hake_n$catches # We take the catch in our data.
#' data$F_input<-hake_n$F # We take the F in our data.
#' # Reference points estimates from ICES stock assessment model:
#' # ICES. 2021. Working Group for the Bay of Biscay and the Iberian Waters Ecoregion
#' # (WGBIE). ICES Scientific Reports. 3:48.1101 pp.
#' data$RP<-list(F_MSY=0.259, B_MSY=207398, MSY=75052, K=NA)
#' # In this case, B_MSY is equal to SSB_MSY, since control$method<-"SSB"
#' # (see control list below).
#' data$years<-hake_n$Year    # Years corresponding to the catch values
#' # (can be different than the years corresponding to SSB or biomass).
#'
#'
#' # Now we define the control.
#'
#' control<-list(
#'   pella = "TRUE")   # Logical. T means that the Pella-Tomlinson model is used.
#'                     #          F means that the Schaefer model is employed.
#'
#' # Finally, we can fit the model
#'
#' knobi_results<-knobi_fit(data,control,plot_out=TRUE,plot_filename="results")
#' knobi_results
#'
#'
#' #### Fitting multispecific KBPM
#'
#' # Below, a multistock approximation aggregating the
#' # northern and southern stocks of sardine is performed.
#'
#' # Firstly, read southern stock data
#' sardine1 <- knobi_data$sardine_s
#'
#' # Secondly, read northern stock data
#' sardine2 <- knobi_data$sardine_n
#'
#' # Extract common years of data in both stocks
#' index <- which(sardine1$Year %in% sardine2$Year)
#' sardine1 <- sardine1[index,]
#'
#' # Create a data.frame where the SSB and the catch are
#' # the sum of such data in the two stocks
#' years<-sardine1$Year
#' sardine <- data.frame(years=years,SSB=sardine1$SSB+sardine2$SSB,
#'                       catch=sardine1$catches+sardine2$catches)
#'
#' # Once the total SSB and catch are available
#' # we follow previous KBPM illustration
#' data<-list()
#' data$SSB<-sardine$SSB
#' data$Catch<-sardine$catch
#' data$years<-sardine$years
#'
#' knobi_results2<-knobi_fit(data)
#' knobi_results2
#'
#'
#' @references
#' Schaefer, M.B. (1954). Some Aspects of the Dynamics of Populations Important to the Management of the Commercial Marine Fisheries. Bulletin of the Inter-American Tropical Tuna Commission. 1:26-56.
#' Pella, J.J., Tomlinson, P.K. (1969). A generalized stock-production model. Bulletin of the Inter-American Tropical Tuna Commission. 13:421–58.
#' MacCall, A. (2002). Use of Known-Biomass Production Models to Determine Productivity of West Coast Groundfish Stocks. North American Journal of Fisheries Management, 22, 272-279.
#'
#' @export

utils::globalVariables(c("y", "f", "f_factor","fr","b","br","b_factor"))

knobi_fit <- function( data, control=NULL, plot_out=FALSE, plot_filename=NULL, plot_dir=NULL){

  if( plot_out == TRUE){

    old_dir <- getwd()

    if ( is.null(plot_dir)) plot_dir <- old_dir
    if ( is.null(plot_filename)) plot_filename <- "knobi_results"

    setwd( plot_dir)

    if ( !plot_filename %in% list.dirs( full.names = FALSE)) dir.create( plot_filename)

    setwd( paste0( plot_dir, "/", plot_filename))

    plot_settings <- list( plot_filename=plot_filename, plot_dir=plot_dir)
    control$plot_settings <- plot_settings

  }

  if( is.null(control)) control<-list()

  set.seed( 123)

  # Check input data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if( is.null(data$Biomass) & is.null(data$SSB)){ stop("Biomass or SSB time series should be provided.")}

  if( !is.null(data$Biomass) & !is.null(data$SSB) & is.null(control$method)){
    stop( "No information about 'control$method', that sets whether the fit is carried using 'SSB' or 'Biomass'.")}

  if( is.null(data$Catch)){ stop("Catch time series should be provided.")}

  C <- data$Catch

  if( is.null(data$years)){ years <- data$years <- 1:length(C)} else { years <- data$years}

  if( is.null(data$Biomass)){ data$Biomass <- NA}
  if( is.null(data$SSB)){ data$SSB<-NA}

  if( is.na(data$Biomass[1])){ control$method<-"SSB"}
  if( is.na(data$SSB[1])){ control$method<-"Biomass"}

  if( length(years) != (length(C))){ stop( 'Length of catch time series is different than the length of years vector.')}

  if( control$method == "Biomass"){ B <- data$Biomass} else { B <- data$SSB}

  if( length(B)!=(length(C)+1)){ warning('The length of the catch time series is reduced according to biomass or SSB time series length.')}
  if( length(B)!=(length(C)) & length(B)!=(length(C)+1)){ stop('The biomass MUST be provided for the same years as the catch or for such years and the next one.')}


  # Compute SP values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SP <- B[-1]; l <- length(B); B_aver <- B[-1]

  for ( i in 1:(l-1)){
    SP[i] <- as.numeric(B[i+1]-B[i]+C[i])
    B_aver[i] <- (B[i+1]+B[i])/2
  }

  B <- B_aver


  # Update correctly dimension ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if ( length(C) != length(B)){
    lb <- length(B)
    C <- C[1:lb]
    years <- years[1:lb]
    data$Catch <- C
    data$years <- years
  }

  if( is.null(data$F_input)){ data$F_input <- NA} else {
    data$F_input <- data$F_input[1:length(data$Catch)]}

  if( is.null(data$Recruitment)){ data$Recruitment <- NA} else {
    data$Recruitment <- data$Recruitment[1:length(data$Catch)]}

  # Save SP, F and average biomass -------------

  data$Average_Biomass <- B_aver
  data$SP <- SP
  F_out <- C/B_aver

  # Care with 0 values
  ind <- which(B_aver==0)
  F_out[ind] <- NA
  data$F_output <- F_out


  # Plots input data --------------

  # Plot trends

  plotInput( data, control$method, plot_out)


  # Fit the model  -----------

  # Define new data.frame with the required info

  df <- data.frame( x = data$Average_Biomass, y = data$SP, Year=years)
  Year <- years

  if( is.null(control$pella)){ control$pella <- FALSE}
  pella <- ifelse( control$pella, TRUE, FALSE)

  start_r <- ifelse( is.null( control$start_r), 0.5, control$start_r)
  start_K <- ifelse( is.null( control$start_K), max(df$x), control$start_K)
  start_p <- ifelse( is.null( control$start_p), 1, control$start_p)

  Data <- list( data=df, start_r=start_r, start_K=start_K, start_p=start_p)

  model <- kbpm_fit( Data, pella = control$pella, model = 'Base')

  attr( model$par,'status') <- NULL

  fit <- list(
    params = model$par,
    optimx = list( value = as.numeric(model$value), convergence = model$convergence, message = model$message))


  # Output plots -----------------------

  x <- df$x
  r <- fit$params['r']
  K <- fit$params['K']

  av <- seq( 0, K, length.out = 3*length(x))
  bv <- predict_model(fit, av, control$pella, 'Base')

  if(!is.null(data$Stock)){ subtitle <- data$Stock} else { subtitle <- NULL}

  df_aux<-data.frame(av,bv)
  if(control$method=="Biomass"){
    xtit<-"SP curve and observed Biomass and SP"
    xaxis<-"Biomass"
    xleg<-"Observed biomass"
  } else {
    xtit<-"SP curve and observed SSB and SP"
    xaxis<-"Spawning Biomass (SSB)"
    xleg<-"Observed SSB"
  }

  vec<-min(df_aux$bv,df$y)
  vec1<-max(df_aux$bv,df$y)

  fit_plot <- ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_line( data=df_aux, ggplot2::aes(x=av,y=bv), linewidth=1.2) + ggplot2::ylim(vec,vec1) +
    ggplot2::geom_point( data = df[c(1,nrow(df)),], ggplot2::aes(x=x,y=y,color=Year)) +
    ggplot2::geom_text( data = df[c(1,nrow(df)),], ggplot2::aes(x=x,y=y,color=Year,
        label=Year), vjust=-1, show.legend = FALSE) +
    ggplot2::geom_point(data=df,ggplot2::aes(x=x,y=y,color=Year)) +
    ggplot2::geom_path(data=df,ggplot2::aes(x=x,y=y,color=Year)) +
    ggplot2::labs(title=xtit,x =xaxis, y = "Surplus Production (SP)") +
    ggplot2::guides(size="none",col=ggplot2::guide_legend(title=xleg)) +
    ggplot2::scale_color_gradient(breaks=c(Year[1],Year[length(Year)])) +
    ggplot2::theme(legend.position = c(.9,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   axis.line=ggplot2::element_line())

  if(!is.null(subtitle)) fit_plot<-fit_plot+ggplot2::labs(subtitle=subtitle)

  print(fit_plot)

  if(plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("fit.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }



  df$C <- C
  df$FM <- data$F_output
  df$base <- predict_model( fit, df$x, pella, 'Base')

  fit$df <- df


  vec<-min(df_aux$bv,df$C)
  vec1<-max(df_aux$bv,df$C)

  fitc_plot<-ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_line(data=df_aux,ggplot2::aes(x=av,y=bv), linewidth=1.2) + ggplot2::ylim(vec,vec1) +
    ggplot2::geom_point(data=df[c(1,nrow(df)),],ggplot2::aes(x=x,y=C,color=Year)) +
    ggplot2::geom_text(data=df[c(1,nrow(df)),],ggplot2::aes(x=x,y=C,color=Year,
         label=Year,vjust=-1), show.legend = FALSE) +
    ggplot2::geom_point(data=df,ggplot2::aes(x=x,y=C,color=Year)) +
    ggplot2::geom_path(data=df,ggplot2::aes(x=x,y=C,color=Year)) +
    ggplot2::labs(title="SP curve and observed catch and SP",
                  x ="Spawning Biomass (SSB)", y = "Surplus Production (SP)") +
    ggplot2::guides(size="none",col=ggplot2::guide_legend(title="Observed catch")) +
    ggplot2::scale_color_gradient(breaks=c(Year[1],Year[length(Year)])) +
    ggplot2::theme(legend.position = c(.9,.85), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  if(!is.null(subtitle)) fitc_plot<-fitc_plot+ggplot2::labs(subtitle=subtitle)

  print(fitc_plot)

  if(plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("fit_catch.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }

  # Compute reference points ---------------
  # (Jacobson et al. 2002) and Winker et al. (2018)

  BRPs <- fit$BRPs <- BRP( fit, control$pella)

  # Extract

  Bmsy <- BRPs['B_MSY']
  Fmsy <- BRPs['F_MSY']

  # Plots RF ---------------------
  # F over years

  F_out <- data$F_output
  if(is.na( data$F_input[1])) {F_inp<-rep(NA,length(F_out))} else {F_inp<-data$F_input}
  if(is.null(data$RP$F_MSY)) {Fmsy_inp<-NA} else {Fmsy_inp<-data$RP$F_MSY}
  Frel_out<-F_out/Fmsy
  Frel_inp<-F_inp/Fmsy_inp

  if(is.null(data$RP$B_MSY)) {Bmsy_inp<-NA} else {Bmsy_inp<-data$RP$B_MSY}
  Brel_out<-B_aver/Bmsy
  Brel_inp<-B_aver/Bmsy_inp

  max_f<-max(c(F_inp,F_out,Fmsy,Fmsy_inp)*1.1,na.rm = TRUE)
  min_f<-min(c(F_inp,F_out,Fmsy,Fmsy_inp),na.rm = TRUE)
  max_b<-max(c(B_aver,Bmsy,Bmsy_inp)*1.1,na.rm = TRUE)
  min_b<-min(c(B_aver,Bmsy,Bmsy_inp),na.rm = TRUE)
  max_fr<-max(c(Frel_inp,Frel_out,1.1),na.rm = TRUE)
  min_fr<-min(c(Frel_inp,Frel_out,0.9),na.rm = TRUE)
  max_br<-max(c(Brel_inp,Brel_out,1.1),na.rm = TRUE)
  min_br<-min(c(Brel_inp,Brel_out,0.9),na.rm = TRUE)

  if(control$method=="SSB"){bfac<-c(rep("SSB from KBPM",length(Year)),rep("original SSB",length(Year)))
  } else {bfac<-c(rep("B from KBPM",length(Year)),rep("original B",length(Year)))}

  if(control$method=="SSB"){
    btit<-"Spawning Biomass (SSB) over time"
    baxis<-"Spawning biomass (SSB)"
    Brtit<-"Relative SSB over time"} else {
      btit<-"Biomass over time"
      baxis<-"Biomass"
      Brtit<-"Relative Biomass over time"}


  plot_df<-data.frame(Year=rep(Year,2),f=c(F_out,F_inp),fr=c(Frel_out,Frel_inp),
                      f_factor=c(rep("F from KBPM",length(Year)),rep("original F",length(Year))),
                      b=rep(B_aver,2),br=c(Brel_out,Brel_inp),b_factor=bfac)


  f_plot<-ggplot2::ggplot(data=subset(plot_df,!is.na(plot_df$f)),ggplot2::aes(x=Year,y=f,color=f_factor)) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min_f,max_f) +
    ggplot2::labs(title="Fishing Mortality (F) over time",y = "Fishing mortality (F)") +
    ggplot2::geom_hline(yintercept=c(Fmsy,Fmsy_inp),linetype="dashed",
                        color = c("#F8766D","#00BFC4"), na.rm=T) +
    ggplot2::annotate("text",x=Year[length(Year)]-1,y=Fmsy,label="Fmsy (KBPM)",
                      color = "#F8766D",size=3,vjust=-1) +
    ggplot2::guides(col=ggplot2::guide_legend(title="")) +
    ggplot2::theme(legend.position = c(.9,.95), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  if(!is.na(Fmsy_inp)){
    f_plot <- f_plot +
      ggplot2::annotate("text",x=Year[1]+1,y=Fmsy_inp,label="original Fmsy",color = "#00BFC4",na.rm=T,size=3,vjust=-1)
  }

  print(f_plot)

  if(plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("F_absolute.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }

  # F over years (relative)

  fr_plot<-ggplot2::ggplot(data=subset(plot_df,!is.na(plot_df$fr)),ggplot2::aes(x=Year,y=fr,color=f_factor)) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min_fr,max_fr) + ggplot2::geom_hline(yintercept=1) +
    ggplot2::labs(title="Relative Fishing Mortality (F) over time",y = "Fishing mortality (F)") +
    ggplot2::guides(col=ggplot2::guide_legend(title="")) +
    ggplot2::theme(legend.position = c(.9,.95), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  print(fr_plot)

  if(plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("F_relative.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }

  # Biomass over years

  b_plot<-ggplot2::ggplot(data=plot_df,ggplot2::aes(x=Year,y=b,color="#F8766D")) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::ylim(min_b,max_b) + ggplot2::labs(title=btit,y=baxis) +
    ggplot2::geom_hline(yintercept=c(Bmsy,Bmsy_inp),linetype="dashed",
                        color = c("#F8766D","#00BFC4"), na.rm=T) +
    ggplot2::annotate("text",x=Year[length(Year)]-1,y=Bmsy,label="Bmsy (KBPM)",
                      color = "#F8766D",size=3,vjust=-1) +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  if(!is.na(Bmsy_inp)){
    b_plot <- b_plot +
      ggplot2::annotate("text",x=Year[1]+1,y=Bmsy_inp,label="original Bmsy",color = "#00BFC4",na.rm=T,size=3,vjust=-1)
  }

  print(b_plot)

  if(plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("B_absolute.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }


  br_plot<-ggplot2::ggplot(data=subset(plot_df,!is.na(plot_df$br)),ggplot2::aes(x=Year,y=br,color=b_factor)) +
    ggplot2::theme_bw() + ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::ylim(min_br,max_br) +
    ggplot2::labs(title=Brtit ,y =baxis) + ggplot2::geom_hline(yintercept=1) +
    ggplot2::guides(col=ggplot2::guide_legend(title="")) +
    ggplot2::theme(legend.position = c(.9,.95), legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5))

  print(br_plot)

  if(plot_out==TRUE){
    p <- grDevices::recordPlot()
    grDevices::jpeg("B_relative.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()
  }

  fit$input <- data
  fit$control <- control

  errors <- kbpm_error( fit, plot_out=plot_out)

  fit$error_table <- errors$error_table
  fit$residuals <- errors$residuals

  if( plot_out==TRUE) setwd(old_dir)

  class(fit)="knobi"

  return(fit)

}
