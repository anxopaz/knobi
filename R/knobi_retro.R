#' @title KBPM retrospective analysis
#'
#' @name knobi_retro
#'
#' @description This function performs a retrospective analysis that evaluates the robustness of the KBPM fit to the systematic deletion of recent data.
#'
#' @param knobi_results A list containing the results of the KBPM fit. Object provided by \code{\link{knobi_fit}} function (main package function).
#' @param nR Number of retrospective patterns. 5 by default. See details.
#' @param yR Optional. Vector of catch time series final years in each one of the retrospective models. See details.
#' @param yR0 Optional. Vector of catch time series starting years in each one of the retrospective models. The same length of 'yR' vector is required. By default, the catch time series is assumed to start in the same year as in the original fit.
#' @param env_results Optional. A list containing the results of the environmental KBPM fit. Object provided by \code{\link{knobi_env}} function.
#' @param plot_out Logical. TRUE means that a file with the plot of the retrospective fits is created. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_dir Optional. Directory to create the folder and save the plots. Required when 'plot_out=TRUE'. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_filename Optional. Name of the folder that will contain the plots. Required when 'plot_out=TRUE'. The default value is the input in the \code{\link{knobi_fit}} function.
#'
#' @details There are different options for defining retrospective fits:
#' (1) Usage of 'nR' argument. This argument specifies the number of retrospective patterns. Furthermore, it is implicit in the use of this argument that the retrospective patterns will consist of the systematic deletion of the last year of data up to the number of years determined by 'nR'.
#' (2) Usage of 'yR' argument. This argument specifies the catch time series final years in each one of the retrospective models, which allows greater flexibility by being able to choose the year from which we delete information. Then, the number of retrospective patterns will be the length of the vector 'yR'. In this option it is possible to set different starting years for the respective time series final years through the argument 'yR0'.
#' If both arguments are provided, the package will use 'yR'.
#'
#' @return A list containing the retrospective analysis: parameter estimates and reference points for each one of the models.
#' The estimated surplus production curves from the retrospective analysis are plotted. The plot is displayed in the plot window and saved (if plot_out=TRUE) in the provided directory or in the current directory.
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
#' library(knobi)
#'
#' # See knobi_fit example to obtain the knobi_results object
#' knobi_retrospectives<-knobi_retro(knobi_results,plot_out=T)  # default nR=5
#' knobi_retrospectives
#'
#' knobi_retro(knobi_results,nR=3)
#' knobi_retro(knobi_results,yR=c(2010,2015))
#' knobi_retro(knobi_results,yR=c(2010,2015),yR0=c(1995,2000))
#'
#' # See knobi_env example to obtain the env_results object
#' knobi_retro(knobi_results,env_results=knobi_environmental,yR=c(2010,2015),yR0=c(1995,2000))
#'
#' }
#'
#' @export


utils::globalVariables(c("y","B","SP"))

knobi_retro <- function( knobi_results, env_results=NULL, nR=5, yR=NULL, yR0=NULL, plot_out=F, plot_filename=NULL, plot_dir=NULL){

  df <- knobi_results$df

  Year <- df$Year
  lastyear <- max(Year)

  x <- df$x

  pars <- knobi_results$params

  pella <- knobi_results$control$pella

  if(plot_out==T){

    old_dir <- getwd()

    if (is.null(plot_dir)) plot_dir <- knobi_results$control$plot_settings$plot_dir
    if (is.null(plot_filename)) plot_filename <- knobi_results$control$plot_settings$plot_filename

    setwd(paste0(plot_dir,"/",plot_filename))

  }

  r <- pars[['r']]
  K <- pars[['K']]
  if ( pella){ p <- pars[['p']]}

  av <- seq( 0, K, length.out = 3*length(x))
  bv <- predict_model( knobi_results, av, pella, 'Base')

  df_aux <- data.frame( av, bv)


  if(is.null(yR)){

    nR <- nR
    yR <- lastyear-1:nR
    yR0 <- rep(Year[1],length(yR))

  } else {

    yR <- yR
    nR <- length(yR)
    if( is.null(yR0)) yR0<-rep( Year[1], length(yR))
    if( length(yR)!=length(yR0)) stop("yR and yR0 must have the same length")

  }

  modelretro <- list()
  brps <- NULL
  params <- NULL
  names_retro <- NULL

  df_plot <- data.frame( B=df_aux$av, SP=df_aux$bv, factor=paste(min(df$Year), "-", lastyear))

  for (i in 1:nR){

    newdf <- subset( df, Year<=yR[i] & Year>=yR0[i])
    iname <- paste( yR0[i], "-", yR[i])
    names_retro <- c( names_retro, iname)

    Data <- list(data=newdf, start_r=0.5, start_K=max(newdf$x), start_p=1)

    ipars <- kbpm_fit( Data, pella = pella, model = "Base")$par
    attr( ipars,'status') <- NULL
    params <- rbind( params, ipars)
    fit <- list( params = ipars)

    ibrps <- BRP( fit, pella)
    fit$brps <- ibrps

    brps <- rbind( brps, ibrps)
    rownames(brps)[i] <- rownames(params)[i] <- iname

    modelretro[[iname]] <- fit

    iK <- ipars['K']
    aretro <- seq(0, iK, length.out = 3*length(x))
    bretro <- predict_model( fit, aretro, pella, 'Base')

    i_df_plot <- data.frame( B=aretro, SP=bretro, factor=iname)
    df_plot <- rbind( df_plot, i_df_plot)

  }

  df_plot$factor <- as.factor( df_plot$factor)

  if( knobi_results$control$method == "SSB"){
    btit<-"SP curve and observed SSB and SP"; baxis<-"Spawning biomass (SSB)"; bleg<-"observed SSB"} else {
      btit<-"SP curve and observed Biomass and SP"; baxis<-"Biomass"; bleg<-"observed biomass"}


  max_y <- max( df_plot$SP, df$y)
  min_y <- min( df_plot$SP, df$y)

  retro_plot <- ggplot2::ggplot() + ggplot2::theme_bw() + ggplot2::ylim(min_y,max_y) +
    ggplot2::geom_point( data=df[c(1,nrow(df)),], ggplot2::aes( x=x, y=y), color='darkblue', size=3, show.legend=FALSE) +
    ggplot2::geom_text( data=df[c(1,nrow(df)),], ggplot2::aes( x=x, y=y, label=Year), color='darkblue', vjust=-1, size=4, show.legend = FALSE) +
    ggplot2::geom_point( data=df, ggplot2::aes( x=x, y=y), color='darkblue', show.legend=FALSE) +
    ggplot2::geom_path( data=df, ggplot2::aes( x=x, y=y), color='darkblue') +
    ggplot2::labs( title=btit, subtitle=knobi_results$input$Stock, x=baxis, y="Surplus Production (SP)") +
    ggplot2::geom_line( data=df_plot, ggplot2::aes( x=B, y=SP, color=factor), linewidth=1) +
    ggplot2::scale_color_manual( values=c(2:(nR+1),1)) +
    ggplot2::guides( size="none",col=ggplot2::guide_legend(title="")) +
    ggplot2::theme( legend.position = c(.89,0.75), legend.background = ggplot2::element_rect(fill = "transparent"),
                    plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5),
                    axis.line=ggplot2::element_line())

  print(retro_plot)

  if(plot_out==TRUE){

    p <- grDevices::recordPlot()
    grDevices::jpeg("fits_retro.jpeg",width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()

    cat(paste0("Plot successfully saved in '",getwd(),"'"),"\n")
    setwd(old_dir)

  }

  brps <- rbind( knobi_results$BRPs, brps)
  params <- rbind( knobi_results$params, params)

  rownames(brps)[1] <- rownames(params)[1] <- paste(min(df$Year), "-", lastyear)

  Retros <- list( BRPs = brps, params = params)


  if( !is.null(env_results)){

    Retros <- list( base = list( BRPs = brps, params = params))

    df <- env_results$df
    evar <- env_results$selected_var

    addres <- env_results$add
    multres <- env_results$mult

    add_brps <- mult_brps <- NULL
    add_params <- mult_params <- NULL

    for (i in 1:nR){

      newdf <- subset( df, Year<=yR[i] & Year>=yR0[i])
      iname <- names_retro[i]

      ienv <- as.matrix( newdf[,evar])

      Data <- list( data=newdf, start_r=0.5, start_K=max(newdf$x), start_p=1, start_c=1)

      iadd_pars <- kbpm_fit( Data, ienv, pella, 'Add')
      imult_pars <- kbpm_fit( Data, ienv, pella, 'Mult')

      add_params <- rbind( add_params, iadd_pars)
      mult_params <- rbind( mult_params, imult_pars)

      iadd_fit <- list( params = iadd_pars)
      imult_fit <- list( params = imult_pars)

      iadd_brps <- BRP( iadd_fit, pella)
      imult_brps <- BRP( imult_fit, pella)

      add_brps <- rbind( add_brps, iadd_brps)
      mult_brps <- rbind( mult_brps, imult_brps)

      rownames(add_brps)[i] <- rownames(add_params)[i] <-
        rownames(mult_brps)[i] <- rownames(mult_params)[i] <- iname

    }

    add_brps <- rbind( addres$BRPs, add_brps)
    mult_brps <- rbind( multres$BRPs, mult_brps)

    add_params <- rbind( addres$params, add_params)
    mult_params <- rbind( multres$params, mult_params)

    rownames(add_brps)[1] <- rownames(add_params)[1] <-
      rownames(mult_brps)[1] <- rownames(mult_params)[1] <- paste(min(df$Year), "-", lastyear)

    Retros$add <- list( BRPs = add_brps, params = add_params)
    Retros$mult <- list( BRPs = mult_brps, params = mult_params)

  }

  return( Retros)

}

