
plotInput <- function( data, method = 'SSB', plot_out){

  years=as.numeric(data$years)
  SP=as.numeric(data$SP)
  C=as.numeric(data$Catch)
  SSB=as.numeric(data$Average_Biomass)

  if(!is.null(data$Stock)){subtitle=data$Stock} else {subtitle=NULL}
  if(is.na(data$F_input[1])==T){
    data=data.frame(Year=as.numeric(data$years),SP=as.numeric(data$SP),
                    C=as.numeric(data$Catch), SSB=as.numeric(data$Average_Biomass))
  } else {
    F=as.numeric(data$F_input)
    data=data.frame(F=as.numeric(data$F_input),Year=as.numeric(data$years),SP=as.numeric(data$SP),
                    C=as.numeric(data$Catch), SSB=as.numeric(data$Average_Biomass))}

  # Panel 1
  # First plot
  if(!is.null(data$F)){
    p1=ggplot2::ggplot(data,ggplot2::aes(x=Year, y=F)) +
      ggplot2::geom_line(data=data[!is.na(data$F),]) +
      ggplot2::geom_point(shape=21, color="black", fill="#56B4E9", size=3) +
      ggplot2::ggtitle("F")+ggplot2::theme_bw()
    if(!is.null(subtitle)){
      p1=p1+ggplot2::labs(subtitle=subtitle)
    }}
  # Second plot

  p2tit <- ifelse( method == 'SSB', "Average SSB", 'Average Biomass')

  p2=ggplot2::ggplot(data,ggplot2::aes(x=Year, y=SSB)) +
    ggplot2::geom_line(data=data[!is.na(data$SSB),]) +
    ggplot2::geom_point(shape=21, color="black", fill="#56B4E9", size=3) +
    ggplot2::ggtitle(p2tit)+ggplot2::theme_bw()
  if(!is.null(subtitle)){
    p2=p2+ggplot2::labs(subtitle=subtitle)
  }


  # Third plot
  p3=ggplot2::ggplot(data,ggplot2::aes(x=Year, y=SP)) +
    ggplot2::geom_line(data=data[!is.na(data$SP),]) +
    ggplot2::geom_point(shape=21, color="black", fill="#56B4E9", size=3) +
    ggplot2::ggtitle("Surplus Production")+ggplot2::theme_bw()
  if(!is.null(subtitle)){
    p3=p3+ggplot2::labs(subtitle=subtitle)
  }
  # Fourth plot
  p4=ggplot2::ggplot(data,ggplot2::aes(x=Year, y=C)) +
    ggplot2::geom_line(data=data[!is.na(data$C),]) +
    ggplot2::geom_point(shape=21, color="black", fill="#56B4E9", size=3) +
    ggplot2::ggtitle("Catch")+ggplot2::theme_bw()
  if(!is.null(subtitle)){
    p4=p4+ggplot2::labs(subtitle=subtitle)
  }
  if(!is.null(data$F)){
    if(plot_out==T){
      grDevices::jpeg("plotInput1.jpeg",width=2500, height=2000,res=300)
      gridExtra::grid.arrange(p1, p2,p3,p4, nrow = 2)
      grDevices::dev.off()}
    gridExtra::grid.arrange(p1, p2,p3,p4, nrow = 2)}

  # Panel 2

  # First plot
  if(!is.null(data$F)){
    ind=is.na(data$C)
    ind=which(ind==TRUE)
    if(length(ind)>0){data1=data[-ind,]} else {data1=data}
    ind1=is.na(data1$F);ind1=which(ind1==TRUE)
    if(length(ind1)>0){data1=data1[-ind1,]} else {data1=data}

    p5=ggplot2::ggplot(data1,ggplot2::aes(x=F, y=C)) +
      ggplot2::geom_smooth(formula = y ~ x, method="loess") +
      ggplot2::geom_point(shape=21, color="black", fill="#56B4E9", size=3) +
      ggplot2::ggtitle("Catch over F")+ggplot2::theme_bw()
    if(!is.null(subtitle)){
      p5=p5+ggplot2::labs(subtitle=subtitle)
    }


    # Second plot

    ind=is.na(data$SSB)
    ind=which(ind==TRUE)
    if(length(ind)>0){data1=data[-ind,]}else {data1=data}
    ind1=is.na(data1$F);ind1=which(ind1==TRUE)
    if(length(ind1)>0){data1=data1[-ind1,]}else {data1=data}

    p6tit <- ifelse( method == 'SSB', "F over average SSB", 'F over average B')

    p6=ggplot2::ggplot(data1,ggplot2::aes(x=SSB, y=F)) +
      ggplot2::geom_smooth(formula = y ~ x, method="loess") +
      ggplot2::geom_point(shape=21, color="black", fill="#56B4E9", size=3) +
      ggplot2::ggtitle(p6tit)+ggplot2::theme_bw()
    if(!is.null(subtitle)){
      p6=p6+ggplot2::labs(subtitle=subtitle)
    }
  }

  # Third plot

  ind=is.na(data$SSB)
  ind=which(ind==TRUE)
  if(length(ind)>0){ data1=data[-ind,]}else {data1=data}
  ind1=is.na(data1$C);ind1=which(ind1==TRUE)
  if(length(ind1)>0){data1=data1[-ind1,]}else {data1=data}

  p7tit <- ifelse( method == 'SSB', "Catch over average SSB", 'Catch over average B')

  p7=ggplot2::ggplot(data1,ggplot2::aes(x=SSB, y=C)) +
    ggplot2::geom_smooth(formula = y ~ x, method="loess") +
    ggplot2::geom_point(shape=21, color="black", fill="#56B4E9", size=3) +
    ggplot2::ggtitle(p7tit)+ggplot2::theme_bw()
  if(!is.null(subtitle)){
    p7=p7+ggplot2::labs(subtitle=subtitle)
  }

  if(is.null(data$F)){
    if(plot_out==T){
      grDevices::jpeg("plotInput.jpeg",width=2500, height=2000,res=300)
      gridExtra::grid.arrange(p2,p3,p4,p7, nrow = 2)
      grDevices::dev.off()}
    gridExtra::grid.arrange(p2,p3,p4,p7, nrow = 2)
  } else {
    if(plot_out==T){
      grDevices::jpeg("plotInput2.jpeg",width=2500, height=1000,res=300)
      gridExtra::grid.arrange(p5,p6,p7, layout_matrix = rbind(c(1,1),c(2,3)))
      grDevices::dev.off()
      gridExtra::grid.arrange(p5,p6,p7, layout_matrix = rbind(c(1,1),c(2,3)))
    } else {
      gridExtra::grid.arrange(p5,p6,p7, layout_matrix = rbind(c(1,1),c(2,3)))}}
}

