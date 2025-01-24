#' @title Input data (sardine and hake stocks and environmental variables)
#'
#' @description
#' 
#' Assessment summary table for the hake stock in ICES subareas 4, 6, and 7, and ICES divisions 3.a, 8.a, 8.b, and 8.d (northern hake stock), as well as for the sardine stocks: northern stock (ICES divisions 8.a, 8.b, and 8.d) and southern stock (ICES divisions 8.c and 9.a).  
#' 
#' The data corresponds to the assessment conducted in 2021 and was derived using the icesSAG package. Environmental data is also reported, including the maximum temperature in Vigo, a city in Galicia (Spain), and the AMO (Atlantic Multi-decadal Oscillation).
#' 
#' @details
#' 
#' The stock's data can be obtained using the following R code:
#' 
#' ```r
#' # Hake northern stock data  
#' 
#' hake_n <- icesSAG::getSAG( stock = "Hake", year = 2021)  
#' 
#' hake_n <- subset( hake_n, summary_data[,17] == "hke.27.3a46-8abd")  
#' 
#' hake_n <- nhake[-nrow(hake_n),]  
#'
#' # Sardine southern stock data  
#' 
#' sardine_s <- icesSAG::getSAG( stock = "pil.27.8c9a", year = 2021)  
#' 
#' sardine_s <- sardine_s[1:43,]  
#' 
#'
#' # Sardine northern stock data  
#' 
#' sardine_n <- icesSAG::getSAG( stock = "pil.27.8abd", year = 2021)  
#' 
#' sardine_n <- sardine_n[1:21,]  
#' ```
#'
#' The 'Tmax' data, i.e. the maximum temperature in Vigo, a city in Galicia (Spain), came from the Meteogalicia website (https://www.meteogalicia.gal), downloaded by hand.  
#' 
#' The 'AMO' data, i.e. The Atlantic Multi-decadal Oscillation is a large-scale climate variability pattern that affects the entire North Atlantic Ocean, was also downloaded by hand from the UCAR website (https://climatedataguide.ucar.edu/).
#'
#' @usage data(knobi_dataset)
#'
#' @format A list where items corresponds to the assessment summary tables for the ICES northern hake stock and the ICES southern and northern sardine stocks, and the environmental data frame.
#'
#' @source `getSAG` function of the "icesSAG" package and the Meteogalicia (https://www.meteogalicia.gal) and UCAR website (https://climatedataguide.ucar.edu/).
#'
"knobi_dataset"
