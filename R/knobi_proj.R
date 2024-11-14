#' @title KBPM projections
#'
#' @name knobi_proj
#'
#' @description Projection of future population and fishery dynamics is carried out for a given set of management targets. More precisely, the function projects the time series of biomass (or spawning biomass) and then the surplus production for a set of future catch or fishing mortality values.
#'
#' @param knobi_results The output object of \code{\link{knobi_fit}} function (main package function).
#' @param env_results Optional. The output object of \code{\link{knobi_env}} function. The environmental KBPM fit is required to forecast considering the environmental variable(s) values in the future.
#' @param Ct Optional. Vector, data frame or matrix that establishes the catch values for each one of the projected years. Different catch scenarios are allowed and can be defined in each of the data frame or matrix columns. Therefore, the length vector (in the case of one scenario) or the number of rows of the data frame or matrix (in the case of multiple scenarios) must be equal to the number of projected years. Projections can be carried out based on a set of future catch values or fishing mortality values, so only one of two arguments, 'Ct' or 'f', is required.
#' @param f Optional. Vector, data frame or matrix that establishes the fishing mortality values for each one of the projected years. Different fishing mortality scenarios are allowed and can be defined in each of the data.frame or matrix columns. Therefore, the length vector (in the case of one scenario) or the number of rows of the data frame or matrix (in the case of multiple scenarios) must be equal to the number of projected years. Projections can be carried out based on a set of future catch values or fishing mortality values, so only one of two arguments, 'Ct' or 'f', is required.
#' @param env Optional. Environmental variable(s) projections required if the environmental fit is considered to forecast the population and fishery dynamics. This fit considers the variable(s) selected in the \code{\link{knobi_env}} function. If the 'multicovar' argument of \code{\link{knobi_env}} is FALSE, tihs argument is a vector, data frame or matrix containing the values of the environmental covariates (unstandardized) for the projection years (rows) and the different catch or fishing mortality settings (columns).  On the other hand, if the 'multicovar' argument of \code{\link{knobi_env}} is TRUE, the current argument must be a list, and each entry must be a data frame or matrix corresponding to each catch or fishing mortality setting containing the values of the environmental covariates for that scenario.
#' @param plot_out Logical. TRUE means that a file with the plot of the forecasts is created. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_dir Optional. Directory to create the folder and save the plots. Required when 'plot_out=TRUE'. The default value is the input in the \code{\link{knobi_fit}} function.
#' @param plot_filename Optional. Name of the folder that will contain the plots. Required when 'plot_out=TRUE'. The default value is the input in the \code{\link{knobi_fit}} function.
#'
#' @return A list containing the projection results. \itemize{
#' \item base: three-dimensional matrix containing the projected time series derived from the KBPM base model (see details in \code{link{knobi_fit}}) for each of the catch or fishing mortality scenarios. The first dimension (rows) corresponds to the years, the second dimension (columns) represent the stock quantities (biomass or SSB, surplus production, F and catches) and the third dimension refers to the projection scenarios.
#' \item add: four-dimensional matrix containing the projected time series derived from the KBPM additive model (see \code{\link{knobi_env}} details) for each of the catch (or fishing mortality) scenarios combined with each of the different environment settings. The first dimension (rows) corresponds to the years, the second dimension (columns) concern the stock quantities (biomass or SSB, surplus production, F and catches), the third dimension refers to the environmental projection scenarios and the fourth dimension corresponds to the catch (or fishing mortality) scenarios. It is only returned if the arguments 'env_results' and 'env' are provided.
#' \item mult: four-dimensional matrix containing the projected time series derived from the KBPM multiplicative model (see \code{\link{knobi_env}} details) for each of the catch (or fishing mortality) scenarios combined with each of the different environment settings. The first dimension (rows) corresponds to the years, the second dimension (columns) concern the stock quantities (biomass or SSB, surplus production, F and catches), the third dimension refers to the environmental projection scenarios and the fourth dimension corresponds to the catch (or fishing mortality) scenarios. It is only returned if the arguments 'env_results' and 'env' are provided.
#' \item df: data frame containing historical and projected biomass (or SSB), catch, SP and fishing mortality values for the different scenarios.
#' \item plots: data frame containing the plots with the projections, for each scenario or combination of scenarios and for each model.}
#' The resulting plots are displayed in the plot window and are also saved (if plot_out="TRUE") in the  provided directory or in the same directory as \code{link{knobi_fit}}.
#' If only the base KBPM is considered in the projections, a plot is presented with a panel reporting the catch, fishing mortality, SP and biomass (or SSB) for each catch or fishing mortality scenario.
#' Furthermore, if the environmental fit is also considered, the plots are displayed with panels for each catch or fishing mortality scenario and for each environmental scenario.
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
#' ### Projecting through catch with no environmental information
#'
#' # Then, create the data frame containing the selected catch for the projected
#' # years. In this illustration, within each scenario, the catch values are
#' # constant through the projected years. Three scenarios are considered:
#' # (i) catch value equal to the last historical catch multiplied by 1,
#' # (ii) last historical catch multiplied by 1.2 and
#' # (iii) last historical catch multiplied by 0.8.
#'
#' catch<-rep(knobi_results$df$C[length(knobi_results$df$C)],5)
#'
#' Ct<-data.frame(catch=catch,
#'               catch08=0.8*catch,
#'               catch12=1.2*catch)
#'
#' # Then, knobi_proj function can be applied
#'
#' knobi_proj(knobi_results, Ct=Ct)
#'
#'
#' ### With environmental information
#'
#' # In this case, in addition to the previous example, the 'knobi_env' example
#' # has to be run at first
#'
#' # We include the future values of the environmental variable(s) in a data
#' # frame containing the environmental covariable values for the projected
#' # years. Three scenarios are considered:
#' # (i) Constant maximum year temperature equal to 19,
#' # (ii) Temperature equal to 19 and then constant raise of 0.5 degrees
#' # (iii) Temperature equal to 19 and then constant raise of 1 degrees
#'
#' env<-data.frame(Tmax_cte=c(19,19,19,19,19),
#'                Tmax_05=c(19,19.5,20,20.5,21),
#'                Tmax_1=c(19,20,21,22,23))
#'
#' # Based on the previous objects we can apply the projection function.
#'
#' knobi_proj(knobi_results, knobi_environmental, Ct=Ct, env=env)
#'
#'
#' ### Through fishing mortality without environmental information
#'
#' # Alternatively, projections can be based on fishing mortality.
#' # The scenarios presented below have been created from the estimated F_msy of
#' # knobi_fit analysis.
#'
#' fmsy<-knobi_results$BRPs['F_MSY']
#' ff<-rep(fmsy,8)
#' f<-data.frame(f=ff,f12=ff*1.2,f08=ff*0.8)
#'
#' knobi_proj(knobi_results, f=f)
#'
#'
#' ### Through fishing mortality with environmental information
#'
#' knobi_proj(knobi_results, f=f[1:5,], env_results=env_results, env=env)
#'
#'
#' # In case of multicovar<-TRUE in knobi_env, a list is required in which
#' # each item is a data frame for each environmental scenario
#'
#' env<-list(climate_1=data.frame(AMO=c(0.2,0.2,0.3,0.3,0.4),
#'                               Tmax=c(19,19,20,20,21)),
#'           climate_2=data.frame(AMO=c(0.2,0.3,0.4,0.5,0.6),
#'                               Tmax=c(19,20,21,22,23)))
#'
#' knobi_proj(knobi_results, knobi_environmental2, Ct=Ct, env=env)
#' }
#'
#' @export


utils::globalVariables(c("B", "Sc", "SP","C","FM","Value","Variable","EnvSc"))

knobi_proj <- function( knobi_results, env_results=NULL, Ct=NULL, f=NULL, env=NULL,
                         plot_out=FALSE, plot_filename=NULL, plot_dir=NULL){


  years <- knobi_results$df$Year
  lastyear <- years[length(years)]

  pella <- knobi_results$control$pella

  if( is.null(Ct) & is.null(f)) stop('You must provide catch or f time series')

  # models ---------------

  if( is.null(f)){ byc <- T
    model <- function(Bt1,Bt,Xt,K,r,p) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))-Xt-Bt1
    model_a <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))+Et%*%c*Bt-Xt-Bt1
    model_m <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))*exp(Et%*%c)-Xt-Bt1
  } else { byc <- F
    model <- function(Bt1,Bt,Xt,K,r,p) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))-Xt*((Bt1+Bt)/2)-Bt1
    model_a <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))+Et%*%c*Bt-Xt*((Bt1+Bt)/2)-Bt1
    model_m <- function(Bt1,Bt,Xt,K,r,p,c,Et) Bt+(r/p)*((Bt1+Bt)/2)*(1-((Bt1+Bt)^p)/(K^p*2^p))*exp(Et%*%c)-Xt*((Bt1+Bt)/2)-Bt1
  }

  if(plot_out==T){

    old_dir <- getwd()

    if (is.null(plot_dir)) plot_dir <- knobi_results$control$plot_settings$plot_dir
    if (is.null(plot_filename)) plot_filename <- knobi_results$control$plot_settings$plot_filename

    setwd(paste0(plot_dir,"/",plot_filename))

  }


  # scenarios & years ------------------

  if( byc) Xt <- as.matrix(Ct) else Xt <- as.matrix(f)

  n_esc <- ncol(Xt)
  if(is.null(colnames(Xt))) colnames(Xt) <- if( n_esc==1){ 'Projection'} else { paste0("Projection_",1:n_esc)}

  sc_names <-colnames(Xt)

  ny <- nrow(Xt)
  nyears <- c( lastyear:c(lastyear+ny))
  pyears <- nyears[-1]
  fyears <- c(years,pyears)
  ly <- length(nyears)


  # init ---------------------

  Bt_ini <- if( knobi_results$control$method == 'SSB') knobi_results$input$SSB else knobi_results$input$Biomass
  B0 <- Bt_ini[length(Bt_ini)]

  params<-knobi_results$params
  r<-params['r']; K<-params['K']; p<-ifelse(pella,params['p'],1)

  base_Bt <- array( c( B0, rep(0,ly-1)), c( ly, n_esc), dimnames = list( nyears+1, sc_names))

  base_Baver <- base_SP <- base_Ct <- base_Ft <-
    array(rep(0,ly-1),c(ly-1,n_esc), dimnames = list( pyears,sc_names))

  df <- pdf0 <- knobi_results$df
  df$Sc <- 'input'
  df$Model <- 'input'
  df$base <- NULL

  df <- rbind( df, data.frame( x = df$x, y = pdf0$base, Year = df$Year, C = df$C, FM = df$FM,
                               Sc = 'input', Model = 'base'))

  # base loop ------------------

  for(j in sc_names){

    for(i in 1:ny){

      Bi <- base_Bt[i,j]

      if( Bi == 0){ Bt1 <- 0} else {

        Xi <- Xt[i,j]

        v <- NULL

        try( v <- stats::uniroot( model, c(0,2*K), Bt=Bi, Xt=Xi, K=K, r=r, p=p)$root, silent=TRUE)
        Bt1 <- ifelse( is.null(v), 0, v)

        if( byc){ base_Ct[i,j] <- Xi; base_Ft[i,j] <- Xi/((Bi+Bt1)/2)} else {
          base_Ft[i,j] <- Xi; base_Ct[i,j] <- Xi*((Bi+Bt1)/2)}}

      if(Bt1<=1e-10){

        base_Bt[c(i+1),j] <- 0
        base_Baver[i,j] <- Bi/2
        base_SP[i,j] <- (r/p)*((Bi)/2)*(1-((Bi)^p)/(K^p*2^p))
        base_Ct[i,j] <- Bi+base_SP[i,j]
        base_Ft[i,j] <- base_Ct[i,j]/(ifelse(Bi==0,1,0.5*Bi))

      } else {

        base_Bt[c(i+1),j]<-Bt1
        base_SP[i,j]<-as.numeric(Bt1-Bi+base_Ct[i,j])
        base_Baver[i,j]<-(Bt1+Bi)/2

      }
    }

    df <- rbind( df, data.frame( x = base_Baver[,j], y = base_SP[,j], Year = pyears, C = base_Ct[,j],
                       FM = base_Ft[,j], Sc = j, Model = 'base'))

  }

  if(any(base_Bt==0)) for(i in which(base_Bt[nrow(base_Bt),]==0))
      warning(paste0('Introduced catch or F in "',sc_names[i],'" scenario lead to stock collapse'))


  # environmental -------------

  if(is.null(env_results)==FALSE & is.null(env)==TRUE) {stop('Environmental data is required')}
  if(is.null(env_results)==TRUE & is.null(env)==FALSE) {stop('Environmental fit results are required')}

  if(is.null(env_results)==FALSE){

    df$EnvSc <- NA

    edf <- env_results$df

    for(mi in c('mult','add'))
      df <- rbind( df, data.frame( x = edf$x, y = edf[[mi]], Year = edf$Year, C = edf$C,
                                      FM = edf$FM, Sc = 'input', Model = mi, EnvSc = NA))


    env <- env

    selvar <- env_results$selected_var
    basevar <- env_results$scaled_var

    nvar <- length(selvar)
    csn <- if(nvar>1) paste0('c',1:nvar) else 'c'

    add <- env_results$add$params
    r_a <- add['r']; K_a <- add['K']; p_a <- ifelse(pella,add['p'],1)
    c_a <- add[csn]

    mult <- env_results$mult$params
    r_m <- mult['r']; K_m <- mult['K']; p_m <- ifelse(pella,mult['p'],1)
    c_m <- mult[csn]

    multicovar <- ifelse( length(env_results$selected_var)>1, T, F)

    if( multicovar){

      nsc <- length(env)
      scnames <- if( is.null(names(env))) paste0( 'Environmental_',1:nsc) else names(env)
      varnames <- colnames( env[[1]]); nvars <- ncol(env[[1]])
      if( any(selvar != varnames)) stop( 'Different environmental variables are provided')
      Et <-  array( NA, dim=c(dim(env[[1]]),nsc), dimnames=(list(NULL,varnames,scnames)))
      for(i in 1:nsc) for(j in 1:nvars)
        Et[,j,i] <- (env[[i]][,j] - attr(basevar,"scaled:center")[j])/attr(basevar,"scaled:scale")[j]

    } else {

      env <- as.matrix(env)
      nsc <- ncol(env)
      scnames <- if( is.null(colnames(env))) paste0( 'Environmental_',1:nsc) else colnames(env)
      Et <-  array( NA, dim=c( nrow(env),1,nsc), dimnames=(list(NULL,NULL,scnames)))
      for(i in 1:nsc) Et[,,i] <- (env[,i] - attr(basevar,"scaled:center"))/attr(basevar,"scaled:scale")

    }

    add_Bt <- mult_Bt <- array( c( B0, rep(0,ly-1)), c( ly, n_esc, nsc), dimnames = list( nyears+1, sc_names, scnames))

    add_Baver <- add_SP <- add_Ct <- add_Ft <- mult_Baver <- mult_SP <- mult_Ct <- mult_Ft <-
      array( 0, c( ly-1, n_esc, nsc), dimnames = list( nyears[-1], sc_names, scnames))


    for(j in sc_names){

      for(n in scnames){

        for(i in 1:ny){

          Bi_a <- add_Bt[i,j,n]
          Bi_m <- mult_Bt[i,j,n]
          Xi <- Xt[i,j]
          Ei <- Et[i,,n]

          if( Bi_a == 0){ Bt1 <- 0} else {

            va <- NULL

            try( va <- stats::uniroot( model_a, c(0,2*(K_a+Ei%*%c_a)),
                              Bt=Bi_a, Xt=Xi, K=K_a, r=r_a, p=p_a, c=c_a, Et=Ei)$root, silent=TRUE)

            Bt1 <- if(is.null(va)) 0 else va

            if( byc){ add_Ct[i,j,n] <- Xi; add_Ft[i,j,n] <- Xi/((Bi_a+Bt1)/2)} else {
              add_Ft[i,j,n] <- Xi; add_Ct[i,j,n] <- Xi*((Bi_a+Bt1)/2)}
            }

          if( Bt1 <= 1e-10){

            add_Bt[c(i+1),j,n] <- 0
            add_SP[i,j,n] <- (r_a/p_a)*((Bi_a)/2)*(1-((Bi_a)^p_a)/(K_a^p_a*2^p_a))
            add_Baver[i,j,n] <- Bi_a/2
            add_Ct[i,j,n] <- Bi_a+add_SP[i,j,n]+Ei%*%c_a*Bi_a
            add_Ft[i,j,n] <- add_Ct[i,j,n]/(ifelse(Bi_a==0,1,0.5*Bi_a))

          } else {

            add_Bt[c(i+1),j,n] <- Bt1
            add_SP[i,j,n] <- as.numeric(Bt1-Bi_a+add_Ct[i,j,n])
            add_Baver[i,j,n] <- (Bt1+Bi_a)/2

          }


          if( Bi_m == 0){ Bt1 <- 0} else {

            vm <- NULL

            try( vm <- stats::uniroot( model_m, c(0,2*(K_a+Ei%*%c_a)),
                              Bt=Bi_m, Xt=Xi, K=K_m, r=r_m, p=p_m, c=c_m, Et=Ei)$root, silent=TRUE)

            Bt1 <- if(is.null(vm)) 0 else vm

            if( byc){ mult_Ct[i,j,n] <- Xi; mult_Ft[i,j,n] <- Xi/((Bi_m+Bt1)/2)} else {
              mult_Ft[i,j,n] <- Xi; mult_Ct[i,j,n] <- Xi*((Bi_m+Bt1)/2)}
            }

          if( Bt1 <= 1e-10){

            mult_Bt[c(i+1),j,n] <- 0
            mult_SP[i,j,n] <- (r_m/p_m)*((Bi_m)/2)*(1-((Bi_m)^p_m)/(K_m^p_m*2^p_m))
            mult_Baver[i,j,n] <- Bi_m/2
            mult_Ct[i,j,n] <- Bi_m+mult_SP[i,j,n]*exp(Ei%*%c_m)
            mult_Ft[i,j,n] <- mult_Ct[i,j,n]/(ifelse(Bi_m==0,1,0.5*Bi_m))

          } else {

            mult_Bt[c(i+1),j,n] <- Bt1
            mult_SP[i,j,n] <- as.numeric(Bt1-Bi_m+mult_Ct[i,j,n])
            mult_Baver[i,j,n] <- (Bt1+Bi_m)/2

          }



        }

        df <- rbind( df, data.frame( x = add_Baver[,j,n], y = add_SP[,j,n], Year = pyears, C = add_Ct[,j,n],
                           FM = add_Ft[,j,n], Sc = j, Model = 'add', EnvSc = n),
                     data.frame( x = mult_Baver[,j,n], y = mult_SP[,j,n], Year = pyears, C = mult_Ct[,j,n],
                           FM = mult_Ft[,j,n], Sc = j, Model = 'mult', EnvSc = n))

      }

  }

  }

  ndf <- df

  colnames(df)[which(colnames(df) == 'x')] <- ifelse( knobi_results$control$method=='SSB','SSB','Biomass')
  colnames(ndf)[which(colnames(ndf) == 'x')] <- 'B'
  colnames(df)[which(colnames(df) == 'y')] <- colnames(ndf)[which(colnames(ndf) == 'y')] <- 'SP'

  dummyb <- subset( ndf, Year == nyears[1] & Model =='base')

  if(is.null(env_results)==FALSE){
    dummya <- subset( ndf, Year == nyears[1] & Model =='add')
    dummym <- subset( ndf, Year == nyears[1] & Model =='mult')}

  for(i in sc_names){
    idummyb <- dummyb
    idummyb$Sc <- i
    ndf <- rbind( ndf, idummyb)
  }

  if(is.null(env_results)==FALSE){ for(i in sc_names){ for(j in scnames){

      idummya <- dummya; idummym <- dummym
      idummya$Sc <- idummym$Sc <-i; idummya$EnvSc <- idummym$EnvSc <-j
      ndf <- rbind( ndf, idummya, idummym)

    }}}

  ndf <- subset( ndf, Year >= nyears[1]-5)

  baseplots <- list()

  bdf <- subset( ndf, Model != 'add' & Model != 'mult')

  baseplots[['B']] <- ggplot2::ggplot( bdf, ggplot2::aes( x = Year, y = B, linetype = Model, color = Sc)) +
    ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
    ggplot2::labs(title = ifelse( knobi_results$control$method=='SSB','SSB projections','Biomass projections'),
                  x = "Year", y = ifelse( knobi_results$control$method=='SSB','SSB (t)','Biomass (t)'), color = "Scenario") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  baseplots[['SP']] <- ggplot2::ggplot( bdf, ggplot2::aes( x = Year, y = SP, linetype = Model, color = Sc)) +
    ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
    ggplot2::labs(title = "SP projections", x = "Year", y = "Surplus Production (t)", color = "Scenario") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  baseplots[['C']] <- ggplot2::ggplot( bdf, ggplot2::aes( x = Year, y = C, linetype = Model, color = Sc)) +
    ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
    ggplot2::labs(title = "Catch projections", x = "Year", y = "Catch (t)", color = "Scenario") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  baseplots[['FM']] <- ggplot2::ggplot( bdf, ggplot2::aes( x = Year, y = FM, linetype = Model, color = Sc)) +
    ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
    ggplot2::labs(title = "F projections", x = "Year", y = "Fishing Mortality", color = "Scenario") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  bdfl <- tidyr::pivot_longer( bdf, cols = c( B, SP, C, FM),
                               names_to = "Variable",
                               values_to = "Value")

  bdfl$Variable[which(bdfl$Variable=='B')] <- if(knobi_results$control$method=='SSB') 'SSB (t)' else 'B (t)'
  bdfl$Variable[which(bdfl$Variable=='C')] <- 'Catch (t)'
  bdfl$Variable[which(bdfl$Variable=='SP')] <- 'SP (t)'
  bdfl$Variable[which(bdfl$Variable=='FM')] <- 'F'

  baseplots[['all']] <- ggplot2::ggplot( bdfl, ggplot2::aes( x = Year, y = Value, linetype = Model, color = Sc)) +
    ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
    ggplot2::facet_grid(rows=dplyr::vars(Variable), scales='free_y') +
    ggplot2::labs(title = "Base KBPM projections", y='', color = "Scenario") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                   plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

  print(baseplots[['all']])

  if (plot_out==TRUE){
    p  <-  grDevices::recordPlot()
    grDevices::jpeg(paste0('base_proj.jpeg'),width=2500, height=2500,res=300)
    grDevices::replayPlot(p)
    grDevices::dev.off()}

  forecast <- list( df=df, base = list( B_aver=base_Baver, Catch=base_Ct, FM=base_Ft, SP=base_SP),
                    plots = list( base_KBPM= baseplots))

  rquants <- c('B_aver','Catch','FM','SP')

  baseres <- array( NA, c( ly-1, n_esc, 4), dimnames = list( nyears[-1], sc_names, rquants))

  for(q in rquants) baseres[,,q] <- forecast$base[[q]]

  forecast$base <- baseres

  if(is.null(env_results)==FALSE){

    plotsbyenv <- list()

    for( j in scnames){

      jdf <- subset( ndf, EnvSc == j | Sc == 'input' & Model != 'input')

      plotsbyenv[[j]][['B']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = B, linetype = Model, color = Sc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = ifelse( knobi_results$control$method=='SSB','SSB projections','Biomass projections'),
                      x = "Year", y = ifelse( knobi_results$control$method=='SSB','SSB (t)','Biomass (t)'), color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      plotsbyenv[[j]][['SP']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = SP, linetype = Model, color = Sc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = "SP projections", x = "Year", y = "Surplus Production (t)", color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      plotsbyenv[[j]][['C']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = C, linetype = Model, color = Sc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = "Catch projections", x = "Year", y = "Catch (t)", color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      plotsbyenv[[j]][['FM']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = FM, linetype = Model, color = Sc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = "F projections", x = "Year", y = "Fishing Mortality", color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      jdfl <- tidyr::pivot_longer( jdf, cols = c( B, SP, C, FM),
                                   names_to = "Variable",
                                   values_to = "Value")

      jdfl$Variable[which(jdfl$Variable=='B')] <- if(knobi_results$control$method=='SSB') 'SSB (t)' else 'B (t)'
      jdfl$Variable[which(jdfl$Variable=='C')] <- 'Catch (t)'
      jdfl$Variable[which(jdfl$Variable=='SP')] <- 'SP (t)'
      jdfl$Variable[which(jdfl$Variable=='FM')] <- 'F'

      plotsbyenv[[j]][['all']] <- ggplot2::ggplot( jdfl, ggplot2::aes( x = Year, y = Value, linetype = Model, color = Sc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = paste0( '"', j, '" scenario projections'), x = "Year", y = "", color = "Scenario") +
        ggplot2::facet_grid(rows=dplyr::vars(Variable), scales='free_y') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      print(plotsbyenv[[j]][['all']])

      if (plot_out==TRUE){
        p  <-  grDevices::recordPlot()
        grDevices::jpeg(paste0(j,'_proj.jpeg'),width=2500, height=2500,res=300)
        grDevices::replayPlot(p)
        grDevices::dev.off()}

    }


    plotsbyC <- list()

    for( j in sc_names){

      jdf <- subset( ndf, Sc == j | Sc == 'input' & Model != 'input')

      jdf$EnvSc[which(is.na(jdf$EnvSc))] <- 'input'

      plotsbyC[[j]][['B']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = B, linetype = Model, color = EnvSc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = ifelse( knobi_results$control$method=='SSB','SSB projections','Biomass projections'),
                      x = "Year", y = ifelse( knobi_results$control$method=='SSB','SSB (t)','Biomass (t)'), color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      plotsbyC[[j]][['SP']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = SP, linetype = Model, color = EnvSc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = "SP projections", x = "Year", y = "Surplus Production (t)", color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      plotsbyC[[j]][['C']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = C, linetype = Model, color = EnvSc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = "Catch projections", x = "Year", y = "Catch (t)", color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      plotsbyC[[j]][['FM']] <- ggplot2::ggplot( jdf, ggplot2::aes( x = Year, y = FM, linetype = Model, color = EnvSc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = "F projections", x = "Year", y = "Fishing Mortality", color = "Scenario") +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      jdfl <- tidyr::pivot_longer( jdf, cols = c( B, SP, C, FM),
                                   names_to = "Variable",
                                   values_to = "Value")

      jdfl$Variable[which(jdfl$Variable=='B')] <- if(knobi_results$control$method=='SSB') 'SSB (t)' else 'B (t)'
      jdfl$Variable[which(jdfl$Variable=='C')] <- 'Catch (t)'
      jdfl$Variable[which(jdfl$Variable=='SP')] <- 'SP (t)'
      jdfl$Variable[which(jdfl$Variable=='FM')] <- 'F'

      plotsbyenv[[j]][['all']] <- ggplot2::ggplot( jdfl, ggplot2::aes( x = Year, y = Value, linetype = Model, color = EnvSc)) +
        ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = nyears[1], linetype = "longdash") +
        ggplot2::labs(title = paste0( '"', j, '" scenario projections'), x = "Year", y = "", color = "Scenario") +
        ggplot2::facet_grid(rows=dplyr::vars(Variable), scales='free_y') +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::theme(legend.position = 'bottom', legend.background = ggplot2::element_rect(fill = "transparent"),
                       plot.title = ggplot2::element_text(hjust = 0.5), axis.line=ggplot2::element_line())

      print(plotsbyenv[[j]][['all']])

      if (plot_out==TRUE){
        p  <-  grDevices::recordPlot()
        grDevices::jpeg(paste0(j,'_proj.jpeg'),width=2500, height=2500,res=300)
        grDevices::replayPlot(p)
        grDevices::dev.off()}

    }

    forecast[['add']] <- list( B_aver=add_Baver, Catch=add_Ct, FM=add_Ft, SP=add_SP)
    forecast[['mult']] <- list( B_aver=mult_Baver, Catch=mult_Ct, FM=mult_Ft, SP=mult_SP)
    forecast$plots[['env']] <- list( byC= plotsbyC, byEnv=plotsbyenv)

    addres <- multres <- array( NA, c( ly-1, n_esc, nsc, 4),
          dimnames = list( nyears[-1], sc_names, scnames, rquants))

    for(q in rquants) addres[,,,q] <- forecast$add[[q]]; multres[,,,q] <- forecast$mult[[q]]

    forecast$add <- addres; forecast$mult <- multres


  }


  if(plot_out==TRUE){
    cat(paste0("\n Plots successfully saved in '",getwd(),"'"),". \n")
    setwd(old_dir)
  }

  class(forecast) <- 'knobi'

  return(forecast)

}

