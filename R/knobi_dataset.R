#' @title Input data (sardine and hake stocks and environmental variables)
#'
#' @description
#' 
#' Assessment summary table for the hake stock in ICES subareas 4, 6, and 7, and ICES divisions 3.a, 8.a, 8.b, and 8.d (northern hake stock), as well as for the sardine stocks: northern stock (ICES divisions 8.a, 8.b, and 8.d) and southern stock (ICES divisions 8.c and 9.a).  
#' 
#' The data corresponds to the assessment conducted in 2021 and was derived using the \pkg{icesSAG} package. Environmental data is also reported, including the Atlantic Multi-decadal Oscillation (AMO) and the North Atlantic Oscillation (NAO) indices.
#' 
#' @details The stock's data can be obtained using the following R code:
#' 
#' \preformatted{
#' # Hake northern stock data:  
#' # hake_n <- icesSAG::getSAG(stock = "hke.27.3a46-8abd", year = 2021)
#' # hake_n <- hake_n[-nrow(hake_n),]  
#' 
#' # Sardine southern stock data:  
#' # sardine_s <- icesSAG::getSAG(stock = "pil.27.8c9a", year = 2021)[1:43,]  
#' 
#' # Sardine northern stock data:  
#' # sardine_n <- icesSAG::getSAG(stock = "pil.27.8abd", year = 2021)[1:21,]  
#' }
#'
#' The AMO data, i.e. the Atlantic Multi-decadal Oscillation index, and the NAO data, i.e. the North Atlantic Oscillation index — both large-scale climate variability patterns that affect the entire North Atlantic Ocean — were downloaded from the UCAR website: \url{https://climatedataguide.ucar.edu/}.
#'
#' @usage data(knobi_dataset)
#'
#' @format A list where items correspond to the assessment summary tables for the ICES northern hake stock, the ICES southern and northern sardine stocks, and the environmental data frame.
#'
#' @source \code{getSAG} function of the \pkg{icesSAG} package and the UCAR website: \url{https://climatedataguide.ucar.edu/}.
#' 
#'
"knobi_dataset"
