#' knobi_data
#'
#' The ICES European hake's Northern stock and the Southern and Northern stock of the sardine (Sardina pilchardus) come from the ICES reports in year 2021, that can be downloaded with the icesSAG package as follows:
#'
#' # Hake northern stock data
#' hake_n <- icesSAG::getSAG( stock = "Hake", year = 2021)
#' hake_n <- subset( hake_n, summary_data[,17] == "hke.27.3a46-8abd")
#' hake_n <- nhake[-nrow(hake_n),]
#'
#' # Sardine southern stock data
#' sardine_s <- icesSAG::getSAG( stock = "pil.27.8c9a", year = 2021)
#' sardine_s <- sardine_s[1:43,]
#'
#' # Sardine northern stock data
#' sardine_n <- icesSAG::getSAG( stock = "pil.27.8abd", year = 2021)
#' sardine_n <- sardine_n[1:21,]
#'
#' The Tmax data came from the Meteogalicia website (https://www.meteogalicia.gal), downloaded by hand.
#' The AMO data was also downloaded by hand from the UCAR website (https://climatedataguide.ucar.edu/).
#'
#' @usage data(knobi_data)
#'
#' @format A list where the items corresponds with the northern hake information, the northern and southern sardine and the environmental information.
#'
#' @source `getSAG` function of the "icesSAG" package and the Meteogalicia and UCAR website.
#'
"knobi_data"
